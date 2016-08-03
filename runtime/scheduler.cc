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


#define __GNU_SOURCE
#include <sched.h>

#define error(...) fprintf(stderr, __VA_ARGS__)

#define DEFAULT_NUM_THREADS 228

Scheduler* Scheduler::global = 0;

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
  ThreadStats(double e, int n) : elapsed_seconds{e}, nthreads{n} {}
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
  return new ThreadStats(elapsed_seconds, ncores);
}

void
SequentialScheduler::runMaster(Task* root)
{
  Logger logger{"scheduler.log"};
  logger.log("message", "cataloging configuration",
             "max_threads", numAvailableCores());

  std::list<Task*> pendingTasks;

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
               "name", Names[task->typeID()]);
    pthread_t thread;
    pthread_create(&thread, NULL, run_task, (void*)task);
    void* retval_v;
    pthread_join(thread, &retval_v);
    ThreadStats* retval = (ThreadStats*) retval_v;
    // read back the elapsed time from the runner and log to file
    double elapsed_seconds = retval->elapsed_seconds;
    int nthreads = retval->nthreads;
    logger.log("end_task", "",
               "name", Names[task->typeID()],
               "elapsed_seconds", elapsed_seconds,
               "nthreads", nthreads,
               "released_listeners", task->getNumListeners());
    task->clearListeners(pendingTasks);
    delete retval;
  }
  logger.log("message", "ALL DONE");
}
