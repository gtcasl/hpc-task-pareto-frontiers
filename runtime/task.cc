#include "task.h"
#include <omp.h>
#include <iostream>
#include <limits>

std::map<int,TaskRunner*> TaskRunner::runners_;
std::map<int,char*> TaskRunner::names_;
std::map<int,std::vector<double> > TaskRunner::times_;
std::map<int,double> TaskRunner::min_times_;
std::map<int,int> TaskRunner::min_threads_;;

Task::Task(int mySize, int typeID) :  
  nthread_(0),
  cpu_state_(0),
  done_(false),
  mySize_(mySize),
  typeID_(typeID)
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
  int num_threads = CPU_COUNT(&cpumask_);
  if(num_threads == 0){
    std::cerr << "Error: task never given thread assignment\n";
    return;
  }
  sched_setaffinity(0, sizeof(cpu_set_t), &cpumask_);
  omp_set_num_threads(num_threads);
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

void
Task::run(int worker, int start_tick)
{
  start_tick_ = start_tick;
  clock_gettime(CLOCK_REALTIME, &start_);
  worker_ = worker;
  int rank = worker + 1;
  MPI_Send(this, mySize_, MPI_BYTE, rank, enqueue_tag, MPI_COMM_WORLD);
  MPI_Irecv(&rc_, 1, MPI_INT, rank, done_tag, MPI_COMM_WORLD, &done_request_);
}

struct timespec
Task::getTime() const
{
  return start_;
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

