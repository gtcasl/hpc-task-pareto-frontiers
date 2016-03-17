#include <iostream>
#include <test.h>

int main(int argc, char** argv)
{
  if (argc < 2){
    std::cerr << "need a test name as first command-line parameter" << std::endl;
    return 1;
  }
  
  Test::fxn f = Test::get(argv[1]);

  return (*f)(argc, argv);
}

