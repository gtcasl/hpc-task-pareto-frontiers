#include "test.h"

typedef Buffer<double> DoublePtr;

enum fxn_id {
  ddot_id = 0
} test_ids;

void my_ddot(int n, Buffer<double> a, Buffer<double> b, Buffer<double> c)
{
  printf("Running ddot %p = %p * %p\n", c.buffer, a.buffer, b.buffer);
  for (int i=0; i < n; ++i){
    c[i] = a[i]*b[i];
  }
}

Task* initDag(int nelems, DoublePtr a, DoublePtr b, DoublePtr c, DoublePtr d)
{
  Task* root = task(my_ddot, ddot_id, std::make_tuple(nelems,a,a,b));
  Task* next = task(my_ddot, ddot_id, std::make_tuple(nelems,a,b,c));
  Task* final = task(my_ddot, ddot_id, std::make_tuple(nelems,b,b,d));
  next->dependsOn(root);
  final->dependsOn(root);
  return root;
}

int run_ddot(int argc, char** argv)
{
  RegisterTask(my_ddot, ddot_id, 
    void, //return
    int, Buffer<double>, Buffer<double>, Buffer<double>); //params

  int nelems = 100;
  int ncopies = 2;

  Scheduler* sch = new BasicScheduler;
  sch->init(argc, argv);

  Buffer<double> a(nelems), b(nelems), c(nelems), d(nelems);
  sch->addNeededBuffers(a, b, c, d);

  sch->allocateHeap(ncopies);
  for (int i=0; i < ncopies; ++i, sch->nextIter()){
    Task* root = 0;
    if (sch->rank() == 0){
      printf("Assigning buffers for iter %d\n", i);
      sch->assignBuffers(a, b, c, d);
      root = initDag(nelems,a,b,c,d);
    }
    sch->run(root);
  }
  sch->deallocateHeap();
  return 0;
}

RegisterTest("ddot", run_ddot);

