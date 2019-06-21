#include "libsbn.hpp"

#include <string>

#include "sbn.h"


extern "C" {
SBNInstance* sbn_NewInstance() { return new SBNInstance; }
void sbn_DeleteInstance(SBNInstance* inst) { delete inst; }

void sbn_ParseFile(SBNInstance* inst, const char* fname) {
  std::string cpp_fname(fname);
  inst->ParseFile(cpp_fname);
}

void sbn_PrintStatus(SBNInstance* inst) { inst->PrintStatus(); }
}
