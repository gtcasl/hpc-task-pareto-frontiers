#include "task.h"
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include <limits>

std::map<int,TaskRunner*> TaskRunner::runners_;
std::map<int,const char*> TaskRunner::names_;
std::map<int,std::vector<double> > TaskRunner::times_;
std::map<int,double> TaskRunner::min_times_;
std::map<int,int> TaskRunner::min_threads_;;
std::map<int,std::vector<double> > TaskRunner::powers_;

int TaskRunner::get_next_least_powerful_num_threads(int id,
                                                    int cur_num_threads) {
  double cur_power = powers_[id][cur_num_threads];
  int new_num_threads = cur_num_threads;
  for(; new_num_threads > 0; --new_num_threads){
    if(powers_[id][new_num_threads] < cur_power){
      return new_num_threads;
    }
  }
  return 0;
}

Task::Task(int mySize, int typeID) :  
  mySize_(mySize),
  cpu_state_(0),
  nthread_(0),
  typeID_(typeID),
  worker_(0),
  done_(false)
{
  CPU_ZERO(&cpumask_);
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
Task::waitDone()
{
  if (done_) return;
  MPI_Wait(&done_request_, MPI_STATUS_IGNORE);
  done_ = true;
}

bool
Task::checkDone()
{
  if (done_) return true;

  int flag;
  MPI_Test(&done_request_, &flag, MPI_STATUS_IGNORE);

  if (flag){
    done_ = true;
    return true;
  } else {
    return false;
  }
}

void
Task::setup()
{
  nthread_ = CPU_COUNT(&cpumask_);
  sched_setaffinity(0, sizeof(cpu_set_t), &cpumask_);
#ifdef _OPENMP
  omp_set_num_threads(nthread_);
#endif
}

void
Task::notifyDone()
{
  int rc = 0;
  int parent = 0;
  MPI_Send(&rc, 1, MPI_INT, parent, done_tag, MPI_COMM_WORLD);
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

double
Task::getTime() const 
{
#if no_timespec
  struct timeval t;
  gettimeofday(&t,0);
  return (t.tv_sec + 1e-6*t.tv_usec);
#else
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  return (t.tv_sec + 1e-9*t.tv_nsec);
#endif
  
}

void
Task::run(int worker, int start_tick)
{
  dep_debug("Running task %p\n", this);
  start_tick_ = start_tick;
  start_ = getTime();
  worker_ = worker;
  int rank = worker + 1;
  MPI_Send(this, mySize_, MPI_BYTE, rank, enqueue_tag, MPI_COMM_WORLD);
  MPI_Irecv(&rc_, 1, MPI_INT, rank, done_tag, MPI_COMM_WORLD, &done_request_);
}

int 
Task::getNumThreads() const
{
  return CPU_COUNT(&cpumask_);
}

double Task::estimateTime() const
{
  double my_time = TaskRunner::get_min_time(typeID_);
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

