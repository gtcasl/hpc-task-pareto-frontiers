#include "test.h"
#include <iostream>

Test::fxn_map* Test::fxns_ = 0;

void
Test::registerNew(const char* name, fxn f)
{
  if (!fxns_) fxns_ = new fxn_map;
  (*fxns_)[name]  = f;
}

Test::fxn
Test::get(const char* name)
{
  std::map<std::string, fxn>::iterator it = fxns_->find(name);
  if (it == fxns_->end()){
    std::cerr << "invalid test name: " << name << std::endl;
    abort();
  }
  return it->second;
}

