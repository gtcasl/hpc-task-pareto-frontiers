#include <test.h>

#include <cstdio>
#include <cassert>
#include <list>



void test(int a, int b, int c)
{
  printf("a=%d, b=%d, c=%d\n", a, b, c);
}

enum function_id {
  test_id,
};

Task* 
dummyTask(int a, int b, int c)
{
  Task* t = task(test, test_id, std::make_tuple(a,b,c));
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

int diamond(int argc, char** argv)
{
  RegisterTask(test,test_id,
    void, int, int, int);

  Scheduler* scheduler = new BasicScheduler;
  
  scheduler->init(argc,argv);
  Task* root = initDag();

  scheduler->run(root);

  return 0;
}

RegisterTest("diamond", diamond);

