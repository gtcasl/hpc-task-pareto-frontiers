#include <exception>
#include <iostream>
#include <cxxabi.h>


int main(int argc, char** argv)
{
  int     status;
  char   *realname;

  for (int i=1; i < argc; ++i){
    const char* name = argv[i];
    realname = abi::__cxa_demangle(name, 0, 0, &status);
    std::cout << name << "\t=> " << realname << "\t: " << status << '\n';
    free(realname);
  }

  return 0;
}
