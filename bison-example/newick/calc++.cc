#include <iostream>
#include "driver.hh"

int main() {

driver drv;


int result = drv.parse_string("(('hi there',x2),x3);");

std::cout << drv.result << '\n';
for (auto& x: drv.taxa) {
    std::cout << x.first << " => " << x.second << '\n';
}

drv.clear();

result = drv.parse_string("(('ddddddddddddddddd',x2),x3);");

std::cout << drv.result << '\n';
for (auto& x: drv.taxa) {
    std::cout << x.first << " => " << x.second << '\n';
}

return result;
//  int res = 0;
//  driver drv;
//  for (int i = 1; i < argc; ++i)
//    if (argv[i] == std::string("-p"))
//      drv.trace_parsing = true;
//    else if (argv[i] == std::string("-s"))
//      drv.trace_scanning = true;
//    else if (!drv.parse(argv[i])) {
//      std::cout << drv.result << '\n';
//      for (auto& x: drv.taxa) {
//          std::cout << x.first << " => " << x.second << '\n';
//      }
//    }
//    else
//      res = 1;
//  return res;
}
