#ifndef __MYWRAPPER_H
#define __MYWRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct MyClass MyClass;

MyClass* newMyClass();

void MyClass_int_set(MyClass* v, int i);

int MyClass_int_get(MyClass* v);

void deleteMyClass(MyClass* v);

#ifdef __cplusplus
}
#endif
#endif

