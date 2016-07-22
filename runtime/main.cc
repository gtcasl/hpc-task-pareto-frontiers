#include <iostream>
#include <string>
#include <memory>
#include <unistd.h>
#include <test.h>

int main(int argc, char** argv)
{
  std::unique_ptr<Scheduler> scheduler = nullptr;
  int opt;
  while((opt = getopt(argc, argv, "s:h")) != -1){
    switch(opt) {
      case 's':{
        std::string argstring{optarg};
        if(argstring == "baseline"){
          scheduler.reset(new BaselineScheduler);
        } else if(argstring == "sequential"){
          scheduler.reset(new SequentialScheduler);
        } else if(argstring == "advanced"){
          scheduler.reset(new AdvancedScheduler);
        } else if(argstring == "profile"){
          scheduler.reset(new ProfilingScheduler);
        } else {
          std::cerr << "Invalid scheduler type. Try 'basic', 'advanced', or 'profile'." << std::endl;
          return -1;
        }
      } break;
      case 'h':
      default:{
        std::cout << "Usage: " << argv[0] << " [-s scheduler-type] <test name> [test params]" << std::endl;
        return 0;
      }
    }
  }

  if(!scheduler){
    scheduler.reset(new AdvancedScheduler);
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

