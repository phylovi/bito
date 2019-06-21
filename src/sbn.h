#ifndef __MYWRAPPER_H
#define __MYWRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SBNInstance SBNInstance;

SBNInstance* sbn_NewInstance();
void sbn_ParseFile(SBNInstance* inst, const char* fname);
void sbn_PrintStatus(SBNInstance* inst);
void sbn_DeleteInstance(SBNInstance* inst);

#ifdef __cplusplus
}
#endif
#endif

