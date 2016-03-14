#include <task.h>
#include <scheduler.h>
#include <cassert>


char taskBuffer[1024];

void 
Scheduler::init(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc_);
  assert(nproc_ > 1);
  nworkers_ = nproc_ - 1;
}

void
Scheduler::terminateWorkers()
{
  int nproc; MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  for (int dst=1; dst < nproc; ++dst){
    int dummy = 42;
    MPI_Send(&dummy, 1, MPI_INT, dst, terminate_tag, MPI_COMM_WORLD);
  }
}

void
Scheduler::runWorker()
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
BasicScheduler::runMaster(Task* root)
{
  std::list<int> availableWorkers;
  for (int i=1; i < nworkers(); ++i){ //leave off 0
    availableWorkers.push_back(i);
  }

  std::list<Task*> runningTasks;
  std::list<Task*> pendingTasks;

  root->run(0); //run the first task on worker 0
  runningTasks.push_back(root);
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

