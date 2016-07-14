#include <test.h>

typedef Buffer<double> DoublePtr;

enum fxn_id {
  my_ddot_id = 0
};

void my_ddot(int n, Buffer<double> a, Buffer<double> b, Buffer<double> c)
{
  printf("Running ddot %p = %p * %p\n", c.buffer, a.buffer, b.buffer);
  for (int i=0; i < n; ++i){
    c[i] = a[i]*b[i];
  }
}

Task* initDag(int nelems, DoublePtr a, DoublePtr b, DoublePtr c, DoublePtr d)
{
  Task* root = new_task(my_ddot, nelems,a,a,b);
  Task* next = new_task(my_ddot, nelems,a,b,c);
  Task* final = new_task(my_ddot, nelems,b,b,d);
  next->dependsOn(root);
  final->dependsOn(root);
  return root;
}

int run_ddot(Scheduler* sch, int argc, char** argv)
{

  RegisterTask(my_ddot,
    void, //return
    int, Buffer<double>, Buffer<double>, Buffer<double>); //params

  int nelems = 100;
  int ncopies = 2;

  Buffer<double> a(nelems), b(nelems), c(nelems), d(nelems);

  sch->allocateHeap(ncopies);
  for (int i=0; i < ncopies; ++i, sch->nextIter()){
    Task* root = 0;
    if (sch->rank() == 0){
      root = initDag(nelems,a,b,c,d);
    }
    sch->run(root);
  }
  sch->deallocateHeap();
  return 0;
}

RegisterTest("ddot", run_ddot);

