#ifndef __MYWRAPPER_H
#define __MYWRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Node Node;

Node* sbn_newNode();
unsigned int sbn_MaxLeafID(Node* n);
void sbn_deleteNode(Node* n);

#ifdef __cplusplus
}
#endif
#endif

