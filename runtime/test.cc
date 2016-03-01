#include "task.h"
#include <cstdio>
#include <cassert>
#include <list>

void test(int a, int b, int c)
{
  printf("a=%d, b=%d, c=%d\n", a, b, c);
}

template <class T> std::map<int,T> impl::FunctionMap<T>::functions_;

enum function_id {
  test_id,
};

char taskBuffer[1024];

static const int terminate_tag = 42;

void
runLoop()
{
  int parent = 0;
  MPI_Status stat;
  while (1){
    MPI_Recv(taskBuffer, sizeof(taskBuffer), MPI_BYTE, parent, 
      MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

    //int count; MPI_Get_count(&stat, MPI_BYTE, &count);
    //printf("got count %d for tag %d\n", count, stat.MPI_TAG);

    if (stat.MPI_TAG == terminate_tag){
      return;
    } else {
      auto t = reinterpret_cast<Task*>(taskBuffer);
      auto runner = TaskRunner::get(t->typeID());
      runner->run(t);
    }
  }
}

void
terminate()
{
  int nproc; MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  for (int dst=1; dst < nproc; ++dst){
    int dummy = 42;
    MPI_Send(&dummy, 1, MPI_INT, dst, terminate_tag, MPI_COMM_WORLD);
  }
}

Task* 
dummyTask(int a, int b, int c)
{
  Task* t = task(test, test_id, std::make_tuple(a,b,c));
  return t;
}

void
runMaster(int nworkers)
{
  //diamond dag
  Task* t1 = dummyTask(0,1,2); 
  Task* t2 = dummyTask(1,3,4);
  Task* t3 = dummyTask(2,4,5);
  Task* t4 = dummyTask(3,7,5);

  t2->dependsOn(t1);
  t3->dependsOn(t1);
  t4->dependsOn(t3);
  t4->dependsOn(t2);

  std::list<int> availableWorkers(nworkers);
  for (int i=1; i < nworkers; ++i){ //leave off 0
    availableWorkers.push_back(i);
  }

  std::list<Task*> runningTasks;
  std::list<Task*> pendingTasks;

  t1->run(0); //run the first task on worker 0
  runningTasks.push_back(t1);
  while (!runningTasks.empty()){
    std::list<Task*>::iterator tmp,
      it = runningTasks.begin(),
      end = runningTasks.end();
    while (it != end){
      tmp = it++;
      Task* t = *tmp;
      if (t->checkDone()){
        availableWorkers.push_front(t->worker());
        runningTasks.erase(tmp);
        t->clearListeners(pendingTasks);
      }
    }

    while (!pendingTasks.empty()){
      if (availableWorkers.empty()){
        break;
      }

      int worker = availableWorkers.front(); 
      availableWorkers.pop_front();

      Task* t = pendingTasks.front();
      pendingTasks.pop_front();

      t->run(worker);
      runningTasks.push_back(t);
    }
  }
}

int main(int argc, char** argv)
{
  RegisterTask(test,test_id,
    void, int, int, int);

  MPI_Init(&argc, &argv);
  int nproc, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  assert(nproc > 1);

  if (rank == 0){
    runMaster(nproc-1);
    terminate();
  } else {
    runLoop();
  }

  return 0;
}

