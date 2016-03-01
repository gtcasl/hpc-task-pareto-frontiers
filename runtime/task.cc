#include "task.h"

std::map<int,TaskRunner*> TaskRunner::runners_;

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
  sched_setaffinity(0, sizeof(cpu_set_t), &cpumask_);
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
  //cpumask stuff would go here
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
Task::run(int worker)
{
  worker_ = worker;
  int rank = worker + 1;
  MPI_Send(this, mySize_, MPI_BYTE, rank, enqueue_tag, MPI_COMM_WORLD);
  MPI_Irecv(&rc_, 1, MPI_INT, rank, done_tag, MPI_COMM_WORLD, &done_request_);
}

