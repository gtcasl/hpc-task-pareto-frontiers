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

int wideDAG(Scheduler* sch, int argc, char** argv)
{
  if(argc != 2){
    std::cerr << "Usage: " << argv[0] << " <# concurrent tasks>" << std::endl;
    return -1;
  }

  int width = atoi(argv[1]);

  Task* root = initDag(width);

  sch->run(root);

  return 0;
}

RegisterTest("wide", wideDAG);

