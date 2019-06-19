#include "driver.hh"
#include "parser.hh"

driver::driver() : trace_parsing(false), trace_scanning(false) {
  id_counter = 0;
}

int driver::parse(const std::string &f) {
  file = f;
  location.initialize(&file);
  scan_begin();
  yy::parser parserObject(*this);
  parserObject.set_debug_level(trace_parsing);
  int res = parserObject();
  scan_end();
  return res;
}
