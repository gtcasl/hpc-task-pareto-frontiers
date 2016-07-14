#include <map>
#include <scheduler.h>
#include <task.h>


class Test { 
 public:
  typedef int (*fxn)(Scheduler*,int,char**);

  static void 
  registerNew(const char* name, fxn f);

  static fxn get(const char* name);

 private:
  typedef std::map<std::string, fxn> fxn_map;
  static fxn_map* fxns_;
};

class TestRegistration {
 public:
  TestRegistration(const char* name, Test::fxn fxn){
    Test::registerNew(name, fxn);
  }
};

#define RegisterTest(name, fxn) TestRegistration register_##fxn(name, fxn)

