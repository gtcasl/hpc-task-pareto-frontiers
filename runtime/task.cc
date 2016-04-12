#include "task.h"
#include <sys/time.h>
#include <omp.h>
#include <iostream>

std::map<int,TaskRunner*> TaskRunner::runners_;
std::map<int,char*> TaskRunner::names_;

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
  nthread_ = CPU_COUNT(&cpumask_);
  sched_setaffinity(0, sizeof(cpu_set_t), &cpumask_);
  omp_set_num_threads(nthread_);
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
  start_tick_ = start_tick;
  start_ = getTime();
  worker_ = worker;
  int rank = worker + 1;
  MPI_Send(this, mySize_, MPI_BYTE, rank, enqueue_tag, MPI_COMM_WORLD);
  MPI_Irecv(&rc_, 1, MPI_INT, rank, done_tag, MPI_COMM_WORLD, &done_request_);
}

