#include <task.h>
#include <scheduler.h>
#include <cassert>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
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
#include <time.h>
#include <miclib.h>
#include <vector>
#include <cstdlib>


#define __GNU_SOURCE
#include <sched.h>

#define error(...) fprintf(stderr, __VA_ARGS__); perror(NULL); abort()

#define DEFAULT_NUM_THREADS 228

Scheduler* Scheduler::global = 0;

Scheduler::~Scheduler(){
  global = 0;
  mic_close_device(mic_device_);
}

Scheduler::Scheduler() :
  cumulative_power_(0),
  num_power_samples_(0),
  max_power_(0),
  power_limit_(1000.0),
  available_power_(power_limit_),
  do_profiling_(false)
{
  if (global){
    fprintf(stderr, "only allowed one instance of a scheduler at a time\n");
    abort();
  }
  global = this;
  // initialize power measurement
  if(mic_open_device(&mic_device_, 0) != E_MIC_SUCCESS){
    std::cerr << "Error: Unable to open MIC" << std::endl;
    abort();
  }
}

void 
Scheduler::init(int argc, char** argv)
{
  int num_threads = DEFAULT_NUM_THREADS;
  auto threads_env = getenv("NUMTHREADS");
  if(threads_env != NULL){
    num_threads = std::atoi(threads_env);
  }
  auto power_limit_env = getenv("POWERLIMIT");
  if(power_limit_env != NULL){
    power_limit_ = std::atol(power_limit_env);
    available_power_ = power_limit_;
  }
  available_cores_ = CPUList(num_threads);
  putenv("MKL_DYNAMIC=FALSE");
}

void
Scheduler::run(Task* root)
{
  runMaster(root);
}

void
Scheduler::finalize()
{
}

double
getTime()
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
  //auto l = [&](const std::pair<int,bool>& a){return a.first == cpu;};
  auto l = std::make_pair(cpu, false);
  auto found = std::find(std::begin(available_cores_.cpus),
                         std::end(available_cores_.cpus),
                         l);
  if(found == std::end(available_cores_.cpus)){
    assert(true && "Error: Trying to free a cpu we don't have access to");
  }
  assert(!found->second && "Error: Trying to double free a cpu");
  found->second = true;
  available_cores_.num_available++;
}

int Scheduler::claimCpu()
{
  assert(available_cores_.num_available > 0 && "Error: should never try to claim a cpu if none are available");
  for(auto& e : available_cores_.cpus){
    if(e.second){
      e.second = false;
      available_cores_.num_available--;
      return e.first;
    }
  }
  assert(false && "Shouldn't get here");
  return -1;
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
      logfile_ << std::endl;
    }

    template<class T, class U, class... Ts>
    void log(const T& t, const U& u, const Ts&... args){
      logfile_ << t << "=" << u << ",";
      log(args...);
    }
};

struct ThreadStats{
  ThreadStats(Task* t, double e, int n) : task{t}, elapsed_seconds{e}, nthreads{n} {}
  Task* task;
  double elapsed_seconds;
  int nthreads;
};

void* run_task(void* params){
  Task* task = (Task*) params;
  double start = getTime();
  task->run();
  double stop = getTime();
  double elapsed_seconds = stop - start;
  int ncores = task->getNumThreads();
  return new ThreadStats(task, elapsed_seconds, ncores);
}

void
SequentialScheduler::runMaster(Task* root)
{
  int iter_num_threads = numAvailableCores();
  Logger logger{"scheduler.log"};
  logger.log("message", "cataloging configuration",
             "max_threads", numAvailableCores());

  double start_time = getTime();
  std::list<Task*> pendingTasks;

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
  if(setitimer(ITIMER_REAL, &work_time, NULL) == -1){
    std::cout << "Unable to set timer" << std::endl;
  }

  cumulative_power_ = 0;
  num_power_samples_ = 0;

  pendingTasks.push_back(root);
  while(!pendingTasks.empty()){
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
    Task* task = pendingTasks.front();
    pendingTasks.pop_front();
    for(int i = 0; i < numAvailableCores(); ++i){
      // get an available cpu
      //int cpu = claimCpu();
      // add it to the task
      task->addCpu(i);
    }
    logger.log("start_task", "",
               "name", Names[task->typeID()],
               "nthreads", iter_num_threads);
    max_power_ = 0;
    cumulative_power_ = 0;
    num_power_samples_ = 0;
    pthread_t thread;
    if(do_profiling_ && task->canProfile()){
      task->enableProfiling();
    }
    pthread_create(&thread, NULL, run_task, (void*)task);
    void* retval_v;
    pthread_join(thread, &retval_v);
    ThreadStats* retval = (ThreadStats*) retval_v;
    // read back the elapsed time from the runner and log to file
    double elapsed_seconds = retval->elapsed_seconds;
    int nthreads = retval->nthreads;
    double avg_power_W = (double) cumulative_power_ / num_power_samples_ / 1000000.0;
    logger.log("end_task", "",
               "name", Names[task->typeID()],
               "elapsed_seconds", elapsed_seconds,
               "nthreads", nthreads,
               "end_time", getTime(),
               "avg_power_W", avg_power_W,
               "max_power_W", max_power_ / 1000000.0,
               "avg_time_s", elapsed_seconds / task->getIters(),
               "released_listeners", task->getNumListeners());
    task->clearListeners(pendingTasks);
    delete retval;
  }

  double end_time = getTime();

  // finalize power measurement on this rank
  struct sigaction sa2;
  sa2.sa_handler = SIG_IGN;
  sigaction(SIGALRM, &sa2, nullptr);
  work_time.it_value.tv_sec = 0;
  work_time.it_value.tv_usec = 0;
  setitimer(ITIMER_REAL, &work_time, NULL);

  logger.log("time_s", end_time - start_time);
  logger.log("message", "ALL DONE");
}
  
void
SimpleScheduler::runMaster(Task* root)
{
  Logger logger{"scheduler.log"};
  logger.log("message", "cataloging configuration",
             "max_threads", numAvailableCores());

  std::list<pthread_t> runningTasks;
  std::list<Task*> pendingTasks;
  std::map<Task*, std::unordered_set<int> > taskCpuAssignments;

  int availableWorkers = 3;

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
    for(auto thread = begin(runningTasks);
        thread != end(runningTasks);){
      void* retval_v;
      if(pthread_tryjoin_np(*thread, &retval_v) == 0){
        // read back the elapsed time from the runner and log to file
        availableWorkers++;
        ThreadStats* retval = (ThreadStats*) retval_v;
        double elapsed_seconds = retval->elapsed_seconds;
        int nthreads = retval->nthreads;
        Task* task = retval->task;
        logger.log("end_task", "",
                   "name", Names[task->typeID()],
                   "elapsed_seconds", elapsed_seconds,
                   "nthreads", nthreads,
                   "end_time", getTime(),
                   "released_listeners", task->getNumListeners());
        task->clearListeners(pendingTasks);
        for(const auto& cpu : taskCpuAssignments[task]){
          returnCpu(cpu);
        }
        available_power_ += Powers[task->typeID()][nthreads];
        taskCpuAssignments[task].clear();
        runningTasks.erase(thread++);
      } else {
        thread++;
      }
    }

    bool increment_tick = true;
    if(!pendingTasks.empty()){
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
      int launches = pendingTasks.size() > availableWorkers ? availableWorkers : pendingTasks.size();
      int threads_per_worker = availableWorkers == 0 ? 0 : numAvailableCores() / launches;
      for(int i = 0; i < launches && threads_per_worker > 0; i++){
        auto task = pendingTasks.front();
        pendingTasks.pop_front();
        // assign task_thread_assignments[task] threads
        for(int i = 0; i < threads_per_worker; ++i){
          // get an available cpu
          int cpu = claimCpu();
          // add it to the task
          task->addCpu(cpu);
          taskCpuAssignments[task].insert(cpu);
        }
        if(increment_tick){
          ++tick_number;
          increment_tick = false;
        }
        logger.log("start_task", "",
                   "name", Names[task->typeID()],
                   "tick", tick_number,
                   "start_time", getTime(),
                   "nthreads", threads_per_worker);
        pthread_t thread;
        pthread_create(&thread, NULL, run_task, (void*)task);
        runningTasks.push_back(thread);
        availableWorkers--;
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
  logger.log("message", "ALL DONE");
}

void
AdvancedScheduler::runMaster(Task* root)
{
  /* TODO: Add logging to document when decisions are made,
   * ie, out of power, out of cores, could use more cores, etc
   */
  Logger logger{"scheduler.log"};
  logger.log("message", "cataloging configuration",
             "power_limit", power_limit_,
             "max_threads", numAvailableCores());
  std::cout << "Starting execution" << std::endl;

  std::list<pthread_t> runningTasks;
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
    for(auto thread = begin(runningTasks); thread != end(runningTasks);){
      void* retval_v;
      if(pthread_tryjoin_np(*thread, &retval_v) == 0){
        // read back the elapsed time from the runner and log to file
        ThreadStats* retval = (ThreadStats*) retval_v;
        double elapsed_seconds = retval->elapsed_seconds;
        int nthreads = retval->nthreads;
        Task* task = retval->task;
        logger.log("end_task", "",
                   "name", Names[task->typeID()],
                   "elapsed_seconds", elapsed_seconds,
                   "nthreads", nthreads,
                   "end_time", getTime(),
                   "released_listeners", task->getNumListeners());
        task->clearListeners(pendingTasks);
        for(const auto& cpu : taskCpuAssignments[task]){
          returnCpu(cpu);
        }
        available_power_ += Powers[task->typeID()][nthreads];
        taskCpuAssignments[task].clear();
        runningTasks.erase(thread++);
      } else {
        thread++;
      }
    }

    bool increment_tick = true;
    if(!pendingTasks.empty()){
      // TODO: walk down pareto frontiers, ensuring both constraints
      std::map<Task*, std::vector<std::tuple<int,double,double>>::iterator> task_assignments;
      int sum_cores = 0;
      double sum_power = 0.0;
      for(const auto& t : pendingTasks){
        task_assignments[t] = begin(Paretos[t->typeID()]);
        sum_cores += std::get<0>(*task_assignments[t]);
        sum_power += std::get<2>(*task_assignments[t]);
      }
      while(sum_cores > numAvailableCores() ||
            sum_power > available_power_){
        // 1. find assignment that reduces execution time the least
        Task* loser = nullptr;
        double delta_t = std::numeric_limits<double>::infinity();
        for(auto& assignment : task_assignments){
          if(std::get<0>(*assignment.second) == 0){
            // this thread already get's no threads, so we can move on
            continue;
          }
          auto task = assignment.first;
          auto cur_makespan = task->estimateTime() -
                              std::get<1>(Paretos[task->typeID()][0]) +
                              std::get<1>(*assignment.second);
          auto next_makespan = task->estimateTime() -
                               std::get<1>(Paretos[task->typeID()][0]) +
                               std::get<1>(*(assignment.second + 1));
          auto delta = next_makespan - cur_makespan;
          if(delta <= delta_t){
            loser = task;
            delta_t = delta;
          }
        }
        // 2. update the assignments
        if(loser){
          sum_cores -= std::get<0>(*task_assignments[loser]);
          sum_power -= std::get<2>(*task_assignments[loser]);
          task_assignments[loser]++;
          sum_cores += std::get<0>(*task_assignments[loser]);
          sum_power += std::get<2>(*task_assignments[loser]);
        } else {
          // unable to find anything that would lower execution time (ie, nothing can run)
          logger.log("message", "Could not find loser, breaking");
          break;
        }
      }

      logger.log("message","reducing available power by the sum of launched tasks",
                 "avail_before", available_power_,
                 "sum_power", sum_power,
                 "avail_after", available_power_ - sum_power);
      available_power_ -= sum_power;
      assert(available_power_ >= 0 && "Error: trying to schedule with more power than we have available");

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
        auto nthreads = std::get<0>(*task_assignments[task]);
        if(nthreads == 0){
          /*
          logger.log("message", "Attempting to schedule task with zero threads",
                     "task", Names[task->typeID()],
                     "tick", tick_number);
                     */
          continue;
        }
        pendingTasks.erase(removed);
        // assign task_thread_assignments[task] threads
        for(int i = 0; i < std::get<0>(*task_assignments[task]); ++i){
          // get an available cpu
          int cpu = claimCpu();
          // add it to the task
          task->addCpu(cpu);
          taskCpuAssignments[task].insert(cpu);
        }
        if(increment_tick){
          ++tick_number;
          increment_tick = false;
        }
        logger.log("start_task", "",
                   "name", Names[task->typeID()],
                   "tick", tick_number,
                   "start_time", getTime(),
                   "nthreads", std::get<0>(*task_assignments[task]));
        pthread_t thread;
        pthread_create(&thread, NULL, run_task, (void*)task);
        runningTasks.push_back(thread);
      }
    }
  } while (!runningTasks.empty() || !pendingTasks.empty());

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
  logger.log("message", "ALL DONE");
}

void 
ProfilingScheduler::runMaster(Task* root)
{
  do_profiling_ = true;
  SequentialScheduler::runMaster(root);
}

void
BaselineScheduler::runMaster(Task* root)
{
  /* TODO: Add logging to document when decisions are made,
   * ie, out of power, out of cores, could use more cores, etc
   */
  Logger logger{"scheduler.log"};
  logger.log("message", "cataloging configuration",
             "max_threads", numAvailableCores());
  std::cout << "Starting execution" << std::endl;

  std::list<pthread_t> runningTasks;
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
    for(auto thread = begin(runningTasks); thread != end(runningTasks);){
      void* retval_v;
      if(pthread_tryjoin_np(*thread, &retval_v) == 0){
        // read back the elapsed time from the runner and log to file
        ThreadStats* retval = (ThreadStats*) retval_v;
        double elapsed_seconds = retval->elapsed_seconds;
        int nthreads = retval->nthreads;
        Task* task = retval->task;
        logger.log("end_task", "",
                   "name", Names[task->typeID()],
                   "elapsed_seconds", elapsed_seconds,
                   "nthreads", nthreads,
                   "end_time", getTime(),
                   "released_listeners", task->getNumListeners());
        task->clearListeners(pendingTasks);
        for(const auto& cpu : taskCpuAssignments[task]){
          returnCpu(cpu);
        }
        taskCpuAssignments[task].clear();
        runningTasks.erase(thread++);
      } else {
        thread++;
      }
    }

    bool increment_tick = true;
    if(!pendingTasks.empty()){

      std::map<Task*,int> task_thread_assignments;
      int sum_s = 0;
      for(const auto& t : pendingTasks){
        int minthreads = std::get<0>(Paretos[t->typeID()][0]);
        task_thread_assignments[t] = minthreads;
        sum_s += minthreads;
      }
      if(sum_s > numAvailableCores()){
        logger.log("message", "Scaling desired threads to those available");
        int new_sum_s = 0;
        for(auto& s : task_thread_assignments){
          s.second = std::floor((double)numAvailableCores() / sum_s * s.second);
          new_sum_s += s.second;
        }
        sum_s = new_sum_s;
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
          est_makespan = t->estimateTime() - std::get<1>(Paretos[t->typeID()][0]) +
            Times[t->typeID()][task_thread_assignments[t]];
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
          est_makespans[min] = min->estimateTime() - std::get<1>(Paretos[min->typeID()][0]) +
            Times[min->typeID()][task_thread_assignments[min]];
        }
        if(task_thread_assignments[max] != 0){
          est_makespans[max] = max->estimateTime() - std::get<1>(Paretos[max->typeID()][0]) +
            Times[max->typeID()][task_thread_assignments[max]];
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
                     "task", Names[task->typeID()],
                     "tick", tick_number);
          continue;
        }
        pendingTasks.erase(removed);
        // assign task_thread_assignments[task] threads
        for(int i = 0; i < task_thread_assignments[task]; ++i){
          // get an available cpu
          int cpu = claimCpu();
          // add it to the task
          task->addCpu(cpu);
          taskCpuAssignments[task].insert(cpu);
        }
        if(increment_tick){
          ++tick_number;
          increment_tick = false;
        }
        logger.log("start_task", "",
                   "name", Names[task->typeID()],
                   "tick", tick_number,
                   "start_time", getTime(),
                   "nthreads", task_thread_assignments[task]);
        pthread_t thread;
        pthread_create(&thread, NULL, run_task, (void*)task);
        runningTasks.push_back(thread);
      }
    }
  } while (!runningTasks.empty() || !pendingTasks.empty());

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
  logger.log("message", "ALL DONE");
}

