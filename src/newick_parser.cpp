#include <chrono>
#include <iostream>
#include "driver.hpp"

int main(int argc, char *argv[]) {
  Driver drv;
  bool print_trees = true;
  bool print_timing = false;
  for (int i = 1; i < argc; ++i)
    if (argv[i] == std::string("-p")) {
      drv.trace_parsing_ = true;
    } else if (argv[i] == std::string("-s")) {
      drv.trace_scanning_ = true;
    } else if (argv[i] == std::string("-t")) {
      print_timing = true;
    } else if (argv[i] == std::string("-q")) {
      print_trees = false;
    } else {
      auto now = std::chrono::high_resolution_clock::now;
      auto t_start = now();
      auto trees = drv.ParseFile(argv[i]);
      auto t_parse = now();

      for (auto tree : *trees) {
        std::vector<unsigned int> leaves(tree->MaxLeafID());
        tree->PreOrder([&leaves](Node* node) {
          if (node->IsLeaf()) leaves.push_back(node->MaxLeafID());
        });
      }
      auto t_traverse = now();
      auto get_duration = [](auto t) {
        return std::chrono::duration<double>(t).count();
      };
      if (print_timing) {
        std::cout << "Parse time: " << get_duration(t_parse - t_start)
                  << "\nTraverse time: " << get_duration(t_traverse - t_parse)
                  << std::endl;
      }
      if (print_trees) {
        for (auto tree : *trees) {
          std::cout << tree->Newick() << std::endl;
        }
      } else {
        std::cout << "Parsed " << trees->size() << " trees.\n";
      }
    }
  return 0;
}
