#include "sbn.h"
#include "sbn.hpp"

extern "C" {
MyClass* newMyClass() { return new MyClass(); }

void MyClass_int_set(MyClass* v, int i) { v->int_set(i); }

int MyClass_int_get(MyClass* v) { return v->int_get(); }

void deleteMyClass(MyClass* v) { delete v; }
}
