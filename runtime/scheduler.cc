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
  power_limit_(std::numeric_limits<double>::infinity()),
  current_power_(0.0),
  do_profiling_(false),
  logfile_{"scheduler.log"},
  tick_number_(0),
  estimated_max_power_(0.0)
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
  logfile_ << std::fixed;
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
  }
  available_cores_ = CPUList(num_threads);
  putenv("MKL_DYNAMIC=FALSE");
}

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
Scheduler::launchTask(Task* task, int nthreads)
{
  if(nthreads == 0){
    return;
  }
  auto task_it = std::find(std::begin(pendingTasks_),
                           std::end(pendingTasks_),
                           task);
  assert(task_it != std::end(pendingTasks_) && "Error: trying to launch task that doesn't exist");
  pendingTasks_.erase(task_it);

  for(int i = 0; i < nthreads; ++i){
    // get an available cpu
    int cpu = claimCpu();
    // add it to the task
    task->addCpu(cpu);
    taskCpuAssignments_[task].insert(cpu);
  }
  current_power_ += Powers[task->typeID()][nthreads];
  log("start_task", "",
      "name", Names[task->typeID()],
      "tick", tick_number_,
      "start_time", getTime(),
      "nthreads", nthreads);
  pthread_t thread;
  pthread_create(&thread, NULL, run_task, (void*)task);
  runningTasks_.push_back(thread);
}

void
Scheduler::run(Task* root)
{
  log("message","cataloging configuration",
      "max_threads", numAvailableCores(),
      "power_limit", power_limit_);
  // initialize power measurement on this rank
  cumulative_power_ = 0;
  num_power_samples_ = 0;
  max_power_ = 0;
  struct sigaction sa;
  sa.sa_sigaction = overflow;
  sa.sa_flags = SA_SIGINFO;
  if(sigaction(SIGALRM, &sa, nullptr) != 0){
    log("error", "Unable to set up signal handler");
    abort();
  }
  struct itimerval work_time;
  work_time.it_value.tv_sec = 0;
  work_time.it_value.tv_usec = 50000; // sleep for 50 ms
  work_time.it_interval.tv_sec = 0;
  work_time.it_interval.tv_usec = 50000; // sleep for 50 ms
  setitimer(ITIMER_REAL, &work_time, NULL);

  double start_time = getTime();

  pendingTasks_.push_back(root);
  bool task_completed = true;
  do{
    // finish completed tasks
    for(auto thread = begin(runningTasks_);
        thread != end(runningTasks_);){
      void* retval_v;
      if(pthread_tryjoin_np(*thread, &retval_v) == 0){
        task_completed = true;
        // read back the elapsed time from the runner and log to file
        ThreadStats* retval = (ThreadStats*) retval_v;
        double elapsed_seconds = retval->elapsed_seconds;
        int nthreads = retval->nthreads;
        Task* task = retval->task;
        auto predicted_time = Times[task->typeID()][nthreads];
        auto percent_error = fabs(predicted_time - elapsed_seconds) / elapsed_seconds * 100.0;
        auto avg_power_W = (double) cumulative_power_ / num_power_samples_ / 1000000.0;
        log("end_task", "",
            "name", Names[task->typeID()],
            "elapsed_seconds", elapsed_seconds,
            "predicted_seconds", predicted_time,
            "time_percent_error", percent_error,
            "nthreads", nthreads,
            "end_time", getTime(),
            "avg_power_W", avg_power_W,
            "max_power_W", max_power_ / 1000000.0,
            "avg_time_s", elapsed_seconds / task->getIters(),
            "released_listeners", task->getNumListeners());
        task->clearListeners(pendingTasks_);
        for(const auto& cpu : taskCpuAssignments_[task]){
          returnCpu(cpu);
        }
        current_power_ -= Powers[task->typeID()][nthreads];
        taskCpuAssignments_[task].clear();
        runningTasks_.erase(thread++);
      } else {
        thread++;
      }
    }
    // launch new tasks
    if(task_completed && !pendingTasks_.empty()){
      tick_number_++;
      task_completed = false;
      tick();
      if(current_power_ > estimated_max_power_){
        estimated_max_power_ = current_power_;
      }
      log("message", "tick stats",
          "tick", tick_number_,
          "estimated_current_power", current_power_);
    }
  } while (!runningTasks_.empty() || !pendingTasks_.empty());

  double end_time = getTime();
  double avg_power_W = (double) cumulative_power_ / num_power_samples_ / 1000000.0;

  // finalize power measurement on this rank
  struct sigaction sa2;
  sa2.sa_handler = SIG_IGN;
  sigaction(SIGALRM, &sa2, nullptr);
  work_time.it_value.tv_sec = 0;
  work_time.it_value.tv_usec = 0;
  setitimer(ITIMER_REAL, &work_time, NULL);

  log("est_max_power", estimated_max_power_);
  log("average_power_W", avg_power_W);
  log("max_power_W", max_power_ / 1000000.0);
  log("time_s", end_time - start_time);
  log("message", "ALL DONE");
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

void
SequentialScheduler::tick()
{
  int iter_num_threads = numAvailableCores();

  Task* task = pendingTasks_.front();
  max_power_ = 0;
  cumulative_power_ = 0;
  num_power_samples_ = 0;
  if(do_profiling_ && task->canProfile()){
    task->enableProfiling();
  }
  launchTask(task, iter_num_threads);
}
  
void
FairScheduler::tick()
{
  int threads_per_worker = numAvailableCores() / pendingTasks_.size();
  int remainder = numAvailableCores() % pendingTasks_.size();
  int numtasks = pendingTasks_.size();
  for(int i = 0; i < numtasks; i++){
    int nthreads = threads_per_worker;
    if(i < remainder){
      nthreads++;
    }
    if(nthreads != 0){
      auto task = pendingTasks_.front();
      launchTask(task, nthreads);
    }
  }
}

void
CoreConstrainedScheduler::tick()
{
  /* TODO: Add logging to document when decisions are made,
   * ie, out of power, out of cores, could use more cores, etc
   */

  std::map<Task*, std::vector<ParetoPoint>::iterator> task_assignments;
  int sum_cores = 0;
  for(const auto& t : pendingTasks_){
    task_assignments[t] = begin(Paretos[t->typeID()]);
    sum_cores += task_assignments[t]->nthreads;
  }
  while(sum_cores > numAvailableCores()){
    // 1. find assignment that minimizes the maximum execution time 
    Task* loser = nullptr;
    double max_t = std::numeric_limits<double>::infinity();
    for(auto& assignment : task_assignments){
      if(assignment.second->nthreads == 0){
        // this thread already get's no threads, so we can move on
        continue;
      }
      auto task = assignment.first;
      auto next_makespan = task->estimateMakespan((assignment.second + 1)->nthreads);
      if(next_makespan <= max_t){
        loser = task;
        max_t = next_makespan;
      }
    }
    // 2. update the assignments
    assert(loser && "Error: could not find losing task, but something should always be able to slow down");

    sum_cores -= task_assignments[loser]->nthreads;
    task_assignments[loser]++;
    sum_cores += task_assignments[loser]->nthreads;
  }

  for(auto& assignment : task_assignments){
    auto task = assignment.first;
    auto nthreads = assignment.second->nthreads;
    launchTask(task, nthreads);
  }
}

void
PowerAwareScheduler::tick()
{
  /* TODO: Add logging to document when decisions are made,
   * ie, out of power, out of cores, could use more cores, etc
   */

  // walk down pareto frontiers, ensuring both constraints
  std::map<Task*, std::vector<ParetoPoint>::iterator> task_assignments;
  int sum_cores = 0;
  double sum_power = 0.0;
  for(const auto& t : pendingTasks_){
    task_assignments[t] = begin(Paretos[t->typeID()]);
    sum_cores += task_assignments[t]->nthreads;
    sum_power += task_assignments[t]->power;
  }
  while(sum_cores > numAvailableCores() ||
        sum_power + current_power_ > power_limit_){
    // 1. find assignment that minimizes the maximum execution time 
    Task* loser = nullptr;
    double max_t = std::numeric_limits<double>::infinity();
    for(auto& assignment : task_assignments){
      if(assignment.second->nthreads == 0){
        // this thread already get's no threads, so we can move on
        continue;
      }
      auto task = assignment.first;
      auto next_makespan = task->estimateMakespan((assignment.second + 1)->nthreads);
      if(next_makespan <= max_t){
        loser = task;
        max_t = next_makespan;
      }
    }
    // 2. update the assignments
    assert(loser && "Error: could not find losing task, but something should always be able to slow down");

    sum_cores -= task_assignments[loser]->nthreads;
    sum_power -= task_assignments[loser]->power;
    task_assignments[loser]++;
    sum_cores += task_assignments[loser]->nthreads;
    sum_power += task_assignments[loser]->power;
  }

  // perform additional power savings if we remain within our power budget
  leverageSlack(task_assignments);

  for(auto& assignment : task_assignments){
    auto nthreads = assignment.second->nthreads;
    launchTask(assignment.first, nthreads);
  }
}

void SlackAwareScheduler::leverageSlack(std::map<Task*,
                                        std::vector<ParetoPoint>::iterator>& task_assignments)
{
  double max_t = 0.0;
  for(auto& assignment : task_assignments){
    if(assignment.second->nthreads == 0){
      continue;
    }
    auto time = assignment.first->estimateMakespan(assignment.second->nthreads);
    if(time > max_t){
      max_t = time;
    }
  }
  for(auto& assignment : task_assignments){
    if(assignment.second->nthreads == 0){
      continue;
    }
    while(1){
      auto new_time = assignment.first->estimateMakespan((assignment.second + 1)->nthreads);
      if(new_time <= max_t){
        log("leveraging_slack","",
            "task", Names[assignment.first->typeID()],
            "time_limit", max_t,
            "new_time", new_time);
        assignment.second++;
      } else {
        break;
      }
    }
  }
}

