#include <stdio.h>
#include "munit.h"
#include "sbn.h"

int main() {
  struct SBNInstance* inst = sbn_NewInstance();
  sbn_PrintStatus(inst);
  sbn_ParseFile(inst, "_ignore/ex.nwk");
  sbn_PrintStatus(inst);
  // munit_assert_int(sbn_MaxLeafID(n), ==, 0);
  sbn_DeleteInstance(inst);
  printf("Successfully deleted instance.\n");
}
