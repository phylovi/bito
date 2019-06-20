#include <iostream>
#include "driver.hpp"

int main(int argc, char *argv[]) {
  Driver drv;
  bool print_trees = true;
  for (int i = 1; i < argc; ++i)
    if (argv[i] == std::string("-p"))
      drv.trace_parsing_ = true;
    else if (argv[i] == std::string("-s"))
      drv.trace_scanning_ = true;
    else if (argv[i] == std::string("-q"))
      print_trees = false;
    else {
      auto trees = drv.ParseFile(argv[i]);
      if (print_trees) {
        for (auto tree : *trees) {
          std::cout << tree->ToNewick() << std::endl;
        }
      } else {
        std::cout << "Parsed " << trees->size() << " trees.\n";
      }
    }
  return 0;
}
