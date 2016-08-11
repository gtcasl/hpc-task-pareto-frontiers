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
std::map<int,std::vector<double> > Powers;
std::map<int,std::vector<ParetoPoint>> Paretos;

int get_next_least_powerful_num_threads(int id,
                                        int cur_num_threads) {
  auto current_power = Powers[id][cur_num_threads];
  if(current_power == 0){
    return 0; // already at least powerful config
  }
  auto idx = std::find_if(begin(Paretos[id]), end(Paretos[id]),
                          [&](const ParetoPoint& e){
                            return e.power < current_power;});
  assert(idx != end(Paretos[id]) && "Error: trying to find negative power?");
  return idx->nthreads;
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
  double my_time = Paretos[typeID_][0].time;
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

std::vector<ParetoPoint>
make_pareto(const std::vector<double>& time,
            const std::vector<double>& power){
  std::vector<ParetoPoint> res;

  std::vector<int> idx(time.size());
  std::iota(begin(idx), end(idx), 0);
  std::sort(begin(idx), end(idx),
            [&](int a, int b){ return time[a] < time[b]; });
  res.emplace_back(idx[0], time[idx[0]], power[idx[0]]);
  for(const auto i : idx){
    if(power[i] < res.back().power){
      res.emplace_back(i, time[i], power[i]);
    }
  }

  return res;
}

