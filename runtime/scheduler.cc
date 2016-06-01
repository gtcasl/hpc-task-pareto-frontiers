#include <task.h>
#include <scheduler.h>
#include <cassert>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <mpi.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h> /* For mode constants */
#include <fcntl.h> /* For O_* constants */
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <sys/types.h>
#include <signal.h>


#define __GNU_SOURCE
#include <sched.h>

#define error(...) fprintf(stderr, __VA_ARGS__)

#define DEFAULT_NUM_THREADS 228

char taskBuffer[1024];

Scheduler* Scheduler::global = 0;

void 
Scheduler::init(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc_);
  //assert(nproc_ > 1);
  nworkers_ = nproc_ - 1;

  int num_threads = DEFAULT_NUM_THREADS;
  auto threads_env = getenv("NUMTHREADS");
  if(threads_env != NULL){
    num_threads = std::atoi(threads_env);
  }
  auto power_limit_env = getenv("POWERLIMIT");
  if(power_limit_env != NULL){
    power_limit_ = std::atoi(power_limit_env);
    available_power_ = power_limit_;
  }
  for(int i = 0; i < num_threads; ++i){
    available_cores_.insert(i);
  }
}

static const char* mmap_fname = "/mmap_heap";

void
Scheduler::addNeededBuffer(BufferBase* buf, size_t size){
  buf->mmap_offset = total_buffer_size_;
  if (size % 4096){
    size_t extra = 4096 - size % 4096;
    size += extra;
  }
  total_buffer_size_ += size;
  buffers_.push_back(buf);
}

void
Scheduler::assignBuffer(BufferBase* buf){
  if (next_copy_ % ncopies_ == 0){
    buf->mmap_offset = buf->mmap_offset % total_buffer_size_;
  }
  else {
    buf->mmap_offset += total_buffer_size_;
  }
  buf->buffer = relocatePointer(buf->mmap_offset);
}

void
Scheduler::resetIter()
{
  next_copy_ = 0;
  for (auto& buf : buffers_){
    assignBuffer(buf);
  }
}

void
Scheduler::nextIter()
{
  next_copy_++;
  for (auto& buf : buffers_){
    assignBuffer(buf);
  }
}

void
Scheduler::allocateHeap(int ncopies)
{

  ncopies_ = ncopies;
  mmap_size_ = total_buffer_size_ * ncopies;

  int fd;
  if(rank_ == 0 || rank_ == 1){
    int res = shm_unlink(mmap_fname);
    fd = shm_open(mmap_fname, O_RDWR | O_CREAT | O_EXCL, S_IRWXU);
    if(fd < 0){
      error("invalid fd %d shm_open on %s: error=%d: rank %d\n",
        fd, mmap_fname, errno, rank_);
      perror(NULL);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if(rank_ != 0 && rank_ != 1){
    fd = shm_open(mmap_fname, O_RDWR, S_IRWXU);
    if (fd < 0){
      error("invalid fd %d shm_open on %s: error=%d: rank %d\n",
        fd, mmap_fname, errno, rank_);
      perror(NULL);
    }
  }

  ftruncate(fd, mmap_size_);

  assert(mmap_size_ % 4096 == 0 && "BAD MMAP SIZE");
  mmap_buffer_ = mmap(NULL, mmap_size_, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
  if (mmap_buffer_ == ((void*)-1)){
    error("bad mmap on shm_open %s:%d: error=%d\n",
      mmap_fname, fd, errno);
    perror(NULL);
  }

  for (auto& buf : buffers_){
    assignBuffer(buf);
  }

}

void
Scheduler::run(Task* root)
{
  if (rank_ == 0){
    runMaster(root);
    terminateWorkers();
  } else {
    runWorker();
  }
}

void
Scheduler::finalize()
{
  MPI_Finalize();
}

void
Scheduler::deallocateHeap()
{
  munmap(mmap_buffer_, mmap_size_);
  shm_unlink(mmap_fname);
}

void
Scheduler::terminateWorkers()
{
  int nproc; MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  for (int dst=1; dst < nproc; ++dst){
    int dummy = 42;
    MPI_Send(&dummy, 1, MPI_INT, dst, terminate_tag, MPI_COMM_WORLD);
  }
}

void
Scheduler::runWorker()
{
  int parent = 0;
  MPI_Status stat;
  while (1){
    MPI_Recv(taskBuffer, sizeof(taskBuffer), MPI_BYTE, parent, 
      MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    int size;
    MPI_Get_count(&stat, MPI_BYTE, &size);
    if (stat.MPI_TAG == terminate_tag){
      return;
    } else {
      auto t = reinterpret_cast<Task*>(taskBuffer);
      auto runner = TaskRunner::get(t->typeID());
      if (!runner){
        fprintf(stderr, "No runner registered for type ID %d\n", t->typeID());
        abort();
      }
      double start = getTime();
      runner->run(t, size);
      double stop = getTime();
      double elapsed_seconds = stop - start;

      // Send elapsed time and num cores to master
      MPI_Request rqst1, rqst2;
      MPI_Isend(&elapsed_seconds, 1, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &rqst1);
      int ncores = t->getNumThreads();
      MPI_Isend(&ncores, 1, MPI_INT, parent, 0, MPI_COMM_WORLD, &rqst2);
    }
  }
}

double
Scheduler::getTime() const 
{
#ifdef no_timespec
  struct timeval start_time;
  gettimeofday(&start_time, NULL);
  return start_time.tv_sec + 1e-6*start_time.tv_usec;
#else
  struct timespec start_time;
  clock_gettime(CLOCK_MONOTONIC, &start_time);
  return start_time.tv_sec + 1e-9*start_time.tv_nsec;
#endif
}

void Scheduler::returnCpu(int cpu)
{
  auto res = available_cores_.insert(cpu);
  assert(res.second && "Error: Trying to double free a cpu");
}

int Scheduler::claimCpu()
{
  assert(numAvailableCores() > 0 && "Error: should never try to claim a cpu if none are available");
  int cpu = *available_cores_.begin();
  available_cores_.erase(available_cores_.begin());
  return cpu;
}

uint32_t Scheduler::readMICPoweruW() const
{
  uint32_t power = 0;
#ifndef no_miclib
  struct mic_power_util_info* pinfo;
  if(mic_get_power_utilization_info(mic_device_, &pinfo) != E_MIC_SUCCESS){
    std::cerr << "Error: Unable to read power utilization info" << std::endl;
    return power;
  }
  if(mic_get_inst_power_readings(pinfo, &power) != E_MIC_SUCCESS){
    std::cerr << "Error: Unable to read power utilization info" << std::endl;
  }
  mic_free_power_utilization_info(pinfo);
#endif
  return power;
}

void Scheduler::overflow(int signum, siginfo_t*, void*){ 
  if(signum == SIGALRM){
    auto power_uW = global->readMICPoweruW();
    global->cumulative_power_ += power_uW;
    ++global->num_power_samples_;
    if(power_uW > global->max_power_){
      global->max_power_ = power_uW;
    }
  }
}

class Logger{
    std::ofstream logfile_;
  public:
    Logger(const char* fname) : logfile_{fname} {
      logfile_ << std::fixed;
    }
    
    void log(){
      logfile_ << "\n";
    }

    template<class T, class U, class... Ts>
    void log(const T& t, const U& u, const Ts&... args){
      logfile_ << t << "=" << u << ",";
      log(args...);
    }
};

void
AdvancedScheduler::runMaster(Task* root)
{
  /* TODO: Add logging to document when decisions are made,
   * ie, out of power, out of cores, could use more cores, etc
   */
  Logger logger{"scheduler.log"};
  logger.log("message", "setting maximum number of threads",
             "nthreads", numAvailableCores());

  std::list<int> availableWorkers;
  for (int i=1; i < nworkers(); ++i){ //leave off 0
    availableWorkers.push_back(i);
  }

  std::list<Task*> runningTasks;
  std::list<Task*> pendingTasks;
  std::map<Task*, std::unordered_set<int> > taskCpuAssignments;

  // initialize power measurement on this rank
  struct sigaction sa;
  sa.sa_sigaction = overflow;
  sa.sa_flags = SA_SIGINFO;
  if(sigaction(SIGALRM, &sa, nullptr) != 0){
    logger.log("error", "Unable to set up signal handler");
    abort();
  }
  struct itimerval work_time;
  work_time.it_value.tv_sec = 0;
  work_time.it_value.tv_usec = 50000; // sleep for 50 ms
  work_time.it_interval.tv_sec = 0;
  work_time.it_interval.tv_usec = 50000; // sleep for 50 ms
  setitimer(ITIMER_REAL, &work_time, NULL);

  double start_time = getTime();

  int tick_number = 0;

  pendingTasks.push_back(root);
  do{
    std::list<Task*>::iterator tmp,
      it = runningTasks.begin(),
      end = runningTasks.end();
    while (it != end){
      tmp = it++;
      Task* t = *tmp;
      if (t->checkDone()){
        // read back the elapsed time from the runner and log to file
        double elapsed_seconds;
        int nthreads;
        MPI_Recv(&elapsed_seconds, 1, MPI_DOUBLE, t->worker() + 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&nthreads, 1, MPI_INT, t->worker() + 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        double start_time = t->getStartTime();
        logger.log("end_task", "",
                   "name", TaskRunner::get_name(t->typeID()),
                   "start_time", start_time,
                   "elapsed_seconds", elapsed_seconds,
                   "nthreads", nthreads,
                   "released_listeners", t->getNumListeners());
        availableWorkers.push_front(t->worker());
        runningTasks.erase(tmp);
        t->clearListeners(pendingTasks);
        for(const auto& cpu : taskCpuAssignments[t]){
          returnCpu(cpu);
        }
        available_power_ += TaskRunner::get_powers(t->typeID())[nthreads - 1];
        taskCpuAssignments[t].clear();
      }
    }

    bool increment_tick = true;
    if(!pendingTasks.empty()){

      std::map<Task*,int> task_thread_assignments;
      int sum_s = 0;
      for(const auto& t : pendingTasks){
        int minthreads = TaskRunner::get_min_threads(t->typeID());
        task_thread_assignments[t] = minthreads;
        sum_s += minthreads;
      }
      if(sum_s > numAvailableCores()){
        logger.log("message", "Scaling desired threads to those available");
        for(auto& s : task_thread_assignments){
          s.second = std::floor((double)numAvailableCores() / sum_s * s.second);
        }
      }else{
        logger.log("message", "Enough threads to go around");
      }

      // Shuffle threads around until we minimize makespan
      std::map<Task*,double> est_makespans;
      for(const auto& t : pendingTasks){
        // Estimate the longest time from this task to completion of the DAG.
        // estimateTime uses the best possible exec time for each task on longest
        // path, including this task, so we have to subtract it out and in its
        // place use the actual predicted execution time for this task with the 
        // number of cores it has been assigned.
        double est_makespan = std::numeric_limits<double>::infinity();
        if(task_thread_assignments[t] != 0){
          est_makespan = t->estimateTime() -
            TaskRunner::get_min_time(t->typeID()) +
            TaskRunner::get_times(t->typeID())[task_thread_assignments[t] - 1];
        }
        est_makespans[t] = est_makespan;
      }
      Task* min = std::min_element(std::begin(est_makespans), std::end(est_makespans),
                                   [](const std::pair<Task*,double>& a,
                                      const std::pair<Task*,double>& b){
                                   return a.second < b.second;
                                   })->first;
      Task* max = std::max_element(std::begin(est_makespans), std::end(est_makespans),
                                   [](const std::pair<Task*,double>& a,
                                      const std::pair<Task*,double>& b){
                                   return a.second < b.second;
                                   })->first;
      double min_makespan = est_makespans[max];
      while(1){
        if(min == max){
          break;
        }
        if(task_thread_assignments[min] == 0){
          break;
        }
        --task_thread_assignments[min];
        ++task_thread_assignments[max];
        est_makespans[min] = std::numeric_limits<double>::infinity();
        est_makespans[max] = std::numeric_limits<double>::infinity();
        if(task_thread_assignments[min] != 0){
          est_makespans[min] = min->estimateTime() -
            TaskRunner::get_min_time(min->typeID()) +
            TaskRunner::get_times(min->typeID())[task_thread_assignments[min] - 1];
        }
        if(task_thread_assignments[max] != 0){
          est_makespans[max] = max->estimateTime() -
            TaskRunner::get_min_time(max->typeID()) +
            TaskRunner::get_times(max->typeID())[task_thread_assignments[max] - 1];
        }
        Task* new_min = std::min_element(std::begin(est_makespans), std::end(est_makespans),
                                         [](const std::pair<Task*,double>& a,
                                            const std::pair<Task*,double>& b){
                                             return a.second < b.second;
                                         })->first;
        Task* new_max = std::max_element(std::begin(est_makespans), std::end(est_makespans),
                                         [](const std::pair<Task*,double>& a,
                                            const std::pair<Task*,double>& b){
                                             return a.second < b.second;
                                         })->first;
        // new config is not an improvement
        if(est_makespans[new_max] >= min_makespan){
          ++task_thread_assignments[min];
          --task_thread_assignments[max];
          break;
        }

        // replace min/max/min_makespan
        min = new_min;
        max = new_max;
        min_makespan = est_makespans[max];
      }

      int sum_p = 0;
      for(const auto& t : pendingTasks){
        if(task_thread_assignments[t] != 0){
          sum_p += TaskRunner::get_powers(t->typeID())[task_thread_assignments[t] - 1];
        }
      }
      if(sum_p > available_power_){
        logger.log("message", "Need to reduce workload to fit within power budget",
                   "available", available_power_,
                   "sum_p", sum_p);
        while(sum_p > available_power_){
          // the loser is the task that slows the least with a delta in power
          // losing task: Task*, #threads, delta t, delta p
          // time goes up, power goes down
          Task* losing_task = nullptr;
          int new_threads;
          double delta_t{0.0};
          double delta_p;
          for(const auto& t : pendingTasks){
            int old_t = task_thread_assignments[t];
            int new_t = TaskRunner::get_next_least_powerful_num_threads(t->typeID(),
                                                                        old_t);
            double old_time = TaskRunner::get_times(t->typeID())[old_t - 1];
            double new_time = TaskRunner::get_times(t->typeID())[new_t - 1];
            if(new_time - old_time < delta_t){
              double old_p = TaskRunner::get_powers(t->typeID())[old_t - 1];
              double new_p = TaskRunner::get_powers(t->typeID())[new_t - 1];
              losing_task = t;
              new_threads = new_t;
              delta_t = new_time - old_time;
              delta_p = new_p - old_p;
            }
          }
          std::cout << "losing: " << losing_task << std::endl;
          for(const auto& t : pendingTasks){
            std::cout << "pending: " << t << std::endl;
          }
          assert(losing_task != nullptr && "Error: Unable to find a losing task");
          assert((new_threads == 0 ||
                 new_threads < task_thread_assignments[losing_task]) &&
                 "Error: by lowering the power we're increasing the number of threads. Uh oh.");
          logger.log("message", "Lowering power",
                     "name", TaskRunner::get_name(losing_task->typeID()),
                     "old_threads", task_thread_assignments[losing_task],
                     "new_threads", new_threads,
                     "delta_t", delta_t,
                     "delta_p", delta_p);
          task_thread_assignments[losing_task] = new_threads;
          sum_p += delta_p;
        }
      } else {
        logger.log("message", "Enough power to go around");
      }

      available_power_ -= sum_p;

      if(increment_tick){
        ++tick_number;
        increment_tick = false;
      }
      /*
       * For each pending task,
       * if task_thread_assignments[task] == 0, continue
       * else
       * remove task from pendingTasks
       * assign it task_thread_assignments[task] threads
       * pop a worker
       * run task on worker
       * add task to runningTasks
       */
      auto task_it = std::begin(pendingTasks);
      while(task_it != std::end(pendingTasks)){
        Task* task = *task_it;
        auto removed = task_it++;
        if(task_thread_assignments[task] == 0){
          logger.log("message", "Attempting to schedule task with zero threads",
                     "task", TaskRunner::get_name(task->typeID()),
                     "tick", tick_number);
          continue;
        }
        assert(!availableWorkers.empty() && 
          "ERROR: Attempting to schedule without enough workers. Use more.");
        pendingTasks.erase(removed);
        // assign task_thread_assignments[task] threads
        for(int i = 0; i < task_thread_assignments[task]; ++i){
          // get an available cpu
          int cpu = claimCpu();
          // add it to the task
          task->addCpu(cpu);
          taskCpuAssignments[task].insert(cpu);
        }
        logger.log("start_task", "",
                   "name", TaskRunner::get_name(task->typeID()),
                   "tick", tick_number,
                   "start_time", getTime(),
                   "nthreads", task_thread_assignments[task]);
        int worker = availableWorkers.front(); 
        availableWorkers.pop_front();
        task->run(worker, tick_number);
        runningTasks.push_back(task);
      }
    }
  } while (!runningTasks.empty());

  double end_time = getTime();
  double avg_power_W = (double) cumulative_power_ / num_power_samples_ / 1000000.0;

  // finalize power measurement on this rank
  struct sigaction sa2;
  sa2.sa_handler = SIG_IGN;
  sigaction(SIGALRM, &sa2, nullptr);
  work_time.it_value.tv_sec = 0;
  work_time.it_value.tv_usec = 0;
  setitimer(ITIMER_REAL, &work_time, NULL);

  logger.log("average_power_W", avg_power_W);
  logger.log("max_power_W", max_power_ / 1000000.0);
  logger.log("time_s", end_time - start_time);
}

void
BasicScheduler::runMaster(Task* root)
{
  std::list<int> availableWorkers;
  for (int i=1; i < nworkers(); ++i){ //leave off 0
    availableWorkers.push_back(i);
  }

  std::list<Task*> runningTasks;
  std::list<Task*> pendingTasks;

  root->run(0,0); //run the first task on worker 0
  runningTasks.push_back(root);
  while (!runningTasks.empty()){
    std::list<Task*>::iterator tmp,
      it = runningTasks.begin(),
      end = runningTasks.end();
    while (it != end){
      tmp = it++;
      Task* t = *tmp;
      if (t->checkDone()){
        availableWorkers.push_front(t->worker());
        runningTasks.erase(tmp);
        t->clearListeners(pendingTasks);
      }
    }

    while (!pendingTasks.empty()){
      if (availableWorkers.empty()){
        break;
      }

      int worker = availableWorkers.front(); 
      availableWorkers.pop_front();

      Task* t = pendingTasks.front();
      pendingTasks.pop_front();

      t->run(worker,0);
      runningTasks.push_back(t);
    }
  }
}
