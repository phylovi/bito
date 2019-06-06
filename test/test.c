#include "munit.h"
#include "sbn.h"

int main(int argc, char **argv) {
  munit_assert_int(Factorial(1), ==, 1);
  munit_assert_int(Factorial(2), ==, 2);
  munit_assert_int(Factorial(3), ==, 6);
  munit_assert_int(Factorial(10), ==, 3628800);
}
