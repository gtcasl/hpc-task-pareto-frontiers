#include <test.h>

enum fxn_id {
  my_ddot_id = 0
};

void my_ddot(int n, double* a, double* b, double* c)
{
  printf("Running ddot %p = %p * %p\n", c, a, b);
  for (int i=0; i < n; ++i){
    c[i] = a[i]*b[i];
  }
}

Task* initDag(int nelems, double* a, double* b, double* c, double* d)
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

  int nelems = 100;
  double a[nelems], b[nelems], c[nelems], d[nelems];

  Task* root = 0;
  root = initDag(nelems,a,b,c,d);
  sch->run(root);

  return 0;
}

RegisterTest("ddot", run_ddot);

