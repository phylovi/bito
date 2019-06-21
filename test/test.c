#include <stdio.h>
#include "munit.h"
#include "sbn.h"

int main(int argc, char** argv) {
  struct Node* n = sbn_newNode();
  munit_assert_int(sbn_MaxLeafID(n), ==, 0);
  sbn_deleteNode(n);
  printf("Tests pass\n");
}
