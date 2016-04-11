#include <test.h>

#include <cstdio>
#include <cassert>
#include <list>



void wide(int a, int b, int c)
{
  printf("a=%d, b=%d, c=%d\n", a, b, c);
}

enum function_id {
  wide_id,
};

extern
Task* 
dummyTask(int a, int b, int c);

Task*
initDag(int width)
{
  Task* root = dummyTask(0,1,2);
  for(int i = 0; i < width; ++i){
    Task* t = dummyTask(i+1,i+2,i+3);
    t->dependsOn(root);
  }
  return root;
}

int wideDAG(int argc, char** argv)
{
  RegisterTask(wide,
    void, int, int, int);

  Scheduler* scheduler = new BasicScheduler;
  
  if(argc != 3){
    std::cerr << "Usage: " << argv[1] << " <# concurrent tasks>" << std::endl;
    return -1;
  }

  int width = atoi(argv[2]);

  scheduler->init(argc,argv);
  Task* root = initDag(width);

  scheduler->run(root);

  return 0;
}

RegisterTest("wide", wideDAG);

