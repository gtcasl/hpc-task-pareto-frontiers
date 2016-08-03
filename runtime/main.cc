#include <iostream>
#include <string>
#include <memory>
#include <unistd.h>
#include <test.h>

int main(int argc, char** argv)
{
  std::unique_ptr<Scheduler> scheduler(new SequentialScheduler());
  int opt;
  while((opt = getopt(argc, argv, "s:h")) != -1){
    switch(opt) {
      case 's':{
      } break;
      case 'h':
      default:{
        std::cout << "Usage: " << argv[0] << " [-s scheduler-type] <test name> [test params]" << std::endl;
        return 0;
      }
    }
  }

  if (optind == argc){
    std::cerr << "need a test name as first command-line parameter" << std::endl;
    return 1;
  }
  
  Test::fxn f = Test::get(argv[optind]);

  scheduler->init(argc, argv);
  auto res = (*f)(scheduler.get(), argc - optind, argv + optind);
  scheduler->finalize();
  return res;
}

