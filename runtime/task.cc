#include "task.h"
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include <limits>
#include <mkl.h>
#include <offload.h>
#include <sys/timerfd.h>
#include <unistd.h>

std::map<int,const char*> Names;
std::map<int,std::vector<double> > Times;
std::map<int,double> MinTimes;
std::map<int,int> MinThreads;
std::map<int,std::vector<double> > Powers;

int get_next_least_powerful_num_threads(int id,
                                        int cur_num_threads) {
  double cur_power = Powers[id][cur_num_threads];
  int new_num_threads = cur_num_threads;
  for(; new_num_threads > 0; --new_num_threads){
    if(Powers[id][new_num_threads] < cur_power){
      return new_num_threads;
    }
  }
  return 0;
}

Task::Task(int typeID, bool isMut) :  
  numIters(0),
  cpu_state_(0),
  nthread_(0),
  typeID_(typeID),
  doProfiling(false),
  isMutating(isMut)
{
  CPU_ZERO(&cpumask_);
  CPU_ZERO(&full_mask_);
  for(int i = 0; i < NUM_THREADS; i++){
    CPU_SET(i, &full_mask_);
  }
}

void
Task::addCpu(int c)
{
  CPU_SET(c, &cpumask_);
}

void
Task::addCpus(int cOffset, int ncores)
{
  int maxCore = cOffset + ncores;
  for (int c=cOffset; c < maxCore; ++c){
    addCpu(c);
  }
}

void
Task::setup()
{
  nthread_ = CPU_COUNT(&cpumask_);
  omp_set_num_threads_target(TARGET_MIC, 0, nthread_);
}

void
Task::clearListeners(std::list<Task*>& ready)
{
  for (auto t : listeners_){
    int join_counter = t->clearDependency(this);
    if (join_counter == 0){
      ready.push_back(t);
    }
  }
  listeners_.clear();
}

int 
Task::getNumThreads() const
{
  return CPU_COUNT(&cpumask_);
}

double Task::estimateTime() const
{
  double my_time = MinTimes[typeID_];
  double children_time = 0.0;
  if(!listeners_.empty()){
    double min_time = std::numeric_limits<double>::infinity();
    for(const auto& listener : listeners_){
      double listener_time = listener->estimateTime();
      if(listener_time < min_time){
        min_time = listener_time;
      }
    }
    children_time = min_time;
  }
  return my_time + children_time;
}

void Task::runProfiling()
{
  setup();
  
  auto fd = timerfd_create(CLOCK_MONOTONIC, TFD_NONBLOCK);
  if(fd == -1){
    std::cerr << "Error: Unable to create timerfd" << std::endl;
  }
  itimerspec it;
  it.it_value.tv_sec = 0;
  it.it_value.tv_nsec = 500000000; // sleep for 500 ms
  it.it_interval.tv_sec = 0;
  it.it_interval.tv_nsec = 0;
  if(timerfd_settime(fd, 0, &it, NULL) == -1){
    std::cerr << "Error: Unable to set timerfd" << std::endl;
  }

  uint64_t buf;
  while(read(fd, &buf, sizeof(buf)) != sizeof(buf)){
    do_run();
  }
  close(fd);
}

