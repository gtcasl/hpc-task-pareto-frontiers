#include <exception>
#include <iostream>
#include <cxxabi.h>

const char* demangle(const char* name)
{
  int     status;
  char   *realname;

  realname = abi::__cxa_demangle(name, 0, 0, &status);
  return realname;
  //std::cout << name << "\t=> " << realname << "\t: " << status << '\n';
   
  return 0;
}

