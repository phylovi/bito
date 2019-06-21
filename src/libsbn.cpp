#include "sbn.h"
#include "sbn.hpp"

extern "C" {
Node* sbn_newNode() { return new Node(0); }

unsigned int sbn_MaxLeafID(Node* n) { return n->MaxLeafID(); }

void sbn_deleteNode(Node* n) { delete n; }
}
