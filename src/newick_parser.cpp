#include <iostream>
#include "driver.hpp"

int main(int argc, char *argv[]) {
  int res = 0;
  Driver drv;
  for (int i = 1; i < argc; ++i)
    if (argv[i] == std::string("-p"))
      drv.trace_parsing_ = true;
    else if (argv[i] == std::string("-s"))
      drv.trace_scanning_ = true;
    else
      drv.ParseFile(argv[i]);
  return res;
}
