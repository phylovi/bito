#include <stdio.h>
#include "munit.h"
#include "sbn.h"

int main(int argc, char** argv) {
  struct MyClass* c = newMyClass();
  MyClass_int_set(c, 3);
  munit_assert_int(MyClass_int_get(c), ==, 3);
  deleteMyClass(c);
  printf("Tests pass\n");
}
