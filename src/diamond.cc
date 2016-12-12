#include <test.h>

#include <cstdio>
#include <cassert>
#include <list>



void test(int a, int b, int c)
{
#pragma offload target(mic:0) \
  in(a) \
  in(b) \
  in(c)
{
  printf("a=%d, b=%d, c=%d\n", a, b, c);
}
}

enum function_id {
  test_id,
};

Task* 
dummyTask(int a, int b, int c)
{
  Task* t = new_task(test, a, b, c);
  return t;
}

Task*
initDag()
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

  return t1;
}

int diamond(Scheduler* sch, int argc, char** argv)
{
  Task* root = initDag();

  sch->run(root);

  return 0;
}

RegisterTest("diamond", diamond);

