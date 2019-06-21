#ifndef SRC_SBN_H_
#define SRC_SBN_H_

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
#endif  // SRC_SBN_H_
