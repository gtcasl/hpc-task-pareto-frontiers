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

#define __GNU_SOURCE
#include <sched.h>

#define error(...) fprintf(stderr, __VA_ARGS__)

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
}

static const char* mmap_fname = "mmap_heap";

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

  if(rank_ == 0){
    shm_unlink(mmap_fname);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  ncopies_ = ncopies;

  int fd = shm_open(mmap_fname, O_RDWR | O_CREAT | O_EXCL, S_IRWXU);
  if (fd < 0){ //oops, someone else already created it
    fd = shm_open(mmap_fname, O_RDWR, S_IRWXU);
  }
  if (fd < 0){
    error("invalid fd %d shm_open on %s: error=%d",
      fd, mmap_fname, errno);
  }

  mmap_size_ = total_buffer_size_ * ncopies;

  ftruncate(fd, mmap_size_);

  assert(mmap_size_ % 4096 == 0 && "BAD MMAP SIZE");
  mmap_buffer_ = mmap(NULL, mmap_size_, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
  if (mmap_buffer_ == ((void*)-1)){
    error("bad mmap on shm_open %s:%d: error=%d",
      mmap_fname, fd, errno);
  }

  //printf("allocated heap pointer %p of size %lu on rank %d\n", 
  //  mmap_buffer_, mmap_size_, rank_);
  
  for (auto& buf : buffers_){
    assignBuffer(buf);
  }

}

void
Scheduler::run(Task* root)
{
  if (rank_ == 0){
    runMaster(root);
  } else {
    runWorker();
  }
}

void
Scheduler::stop()
{
  if (rank_ == 0){
    terminateWorkers();
  } else {
    //do nothing
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

void
BasicScheduler::runMaster(Task* root)
{
  /* TODO: Add logging to document when decisions are made,
   * ie, out of power, out of cores, could use more cores, etc
   */
  std::ofstream logfile{"scheduler.log"};
  logfile << "task\ttick\tstart\telapsed\tthreads\n";

  std::list<int> availableWorkers;
  for (int i=1; i < nworkers(); ++i){ //leave off 0
    availableWorkers.push_back(i);
  }

  std::list<Task*> runningTasks;
  std::list<Task*> pendingTasks;
  std::map<Task*, std::unordered_set<int> > taskCpuAssignments;

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
        logfile << std::fixed
                << TaskRunner::get_name(t->typeID()) << "\t"
                << t->getStartTick() << "\t"
                << start_time << "\t"
                << elapsed_seconds << "\t"
                << nthreads << "\n";
        availableWorkers.push_front(t->worker());
        runningTasks.erase(tmp);
        t->clearListeners(pendingTasks);
        for(const auto& cpu : taskCpuAssignments[t]){
          returnCpu(cpu);
        }
        taskCpuAssignments[t].clear();
      }
    }

    bool increment_tick = true;
    if(!pendingTasks.empty()){
      assert(availableWorkers.size() >= pendingTasks.size() && 
        "ERROR: Attempting to schedule without enough workers. Use more.");

      std::map<Task*,int> s_star;
      int sum_s = 0;
      for(const auto& t : pendingTasks){
        int minthreads = TaskRunner::get_min_threads(t->typeID());
        s_star[t] = minthreads;
        sum_s += minthreads;
      }
      if(sum_s > numAvailableCores()){
        std::cout << "Scaling desired threads to those available\n";
        for(auto& s : s_star){
          s.second = std::floor((double)numAvailableCores() / sum_s * s.second);
        }
      }else{
        //std::cout << "Enough threads to go around\n";
      }

      // Shuffle threads around until we minimize makespan
      std::map<Task*,double> blevels;
      for(const auto& t : pendingTasks){
        // Estimate the longest time from this task to completion of the DAG.
        // estimateTime uses the best possible exec time for each task on longest
        // path, including this task, so we have to subtract it out and in its
        // place use the actual predicted execution time for this task with the 
        // number of cores it has been assigned.
        double blevel = t->estimateTime() -
          TaskRunner::get_min_time(t->typeID()) +
          TaskRunner::get_times(t->typeID())[s_star[t]];
        blevels[t] = blevel;
      }
      while(1){
        Task* min = std::min_element(std::begin(blevels), std::end(blevels),
                                    [](const std::pair<Task*,double>& a,
                                       const std::pair<Task*,double>& b){
                                        return a.second < b.second;
                                    })->first;
        Task* max = std::max_element(std::begin(blevels), std::end(blevels),
                                    [](const std::pair<Task*,double>& a,
                                       const std::pair<Task*,double>& b){
                                        return a.second < b.second;
                                    })->first;
        if(min == max){
          break;
        }
        if(s_star[min] == 0){
          break;
        }
        double new_bmin = min->estimateTime() -
          TaskRunner::get_min_time(min->typeID()) +
          TaskRunner::get_times(min->typeID())[s_star[min]-1];
        double new_bmax = max->estimateTime() -
          TaskRunner::get_min_time(max->typeID()) +
          TaskRunner::get_times(max->typeID())[s_star[max]+1];
        if(blevels[max] <= new_bmax){
          break;
        }
        --s_star[min];
        ++s_star[max];
        blevels[min] = new_bmin;
        blevels[max] = new_bmax;

      }

      if(increment_tick){
        ++tick_number;
        increment_tick = false;
      }
      /*
       * For each pending task,
       * if s_star[task] == 0, continue
       * else
       * remove task from pendingTasks
       * assign it s_star[task] threads
       * pop a worker
       * run task on worker
       * add task to runningTasks
       */
      auto task_it = std::begin(pendingTasks);
      while(task_it != std::end(pendingTasks)){
        Task* task = *task_it;
        auto removed = task_it++;
        if(s_star[task] == 0){
          continue;
        }
        pendingTasks.erase(removed);
        // assign s_star[task] threads
        for(int i = 0; i < s_star[task]; ++i){
          // get an available cpu
          int cpu = claimCpu();
          // add it to the task
          task->addCpu(cpu);
          taskCpuAssignments[task].insert(cpu);
        }
        int worker = availableWorkers.front(); 
        availableWorkers.pop_front();
        task->run(worker, tick_number);
        runningTasks.push_back(task);
      }
    }
  } while (!runningTasks.empty());
}

