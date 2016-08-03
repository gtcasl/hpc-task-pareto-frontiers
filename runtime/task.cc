#include "task.h"
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include <limits>
#include <mkl.h>

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

Task::Task(int typeID) :  
  cpu_state_(0),
  nthread_(0),
  typeID_(typeID),
  doProfiling(false)
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
  mkl_set_num_threads(nthread_);
#pragma offload target(mic:0)
{
  //mkl_set_dynamic(0);
  sched_setaffinity(0, sizeof(cpu_set_t), &cpumask_);
}

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

bool impl::done;

void impl::wakeup(int signum, siginfo_t*, void*){
  if(signum == SIGALRM){
    done = true;
  }
}

