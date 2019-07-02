# libsbn
C-interface library implementing subsplit Bayes networks for phylogenetic posterior density estimation


## Dependencies

### System

* compiler that can handle C++14
* flex (it should install the `libfl` library)


### Suggest installing with conda

* scons
* pybind11
* pytest
* beagle (`conda install -c bioconda beagle-lib`)


## Building

* If you just want to build the library and binaries, use `scons`.
* If you have modified the lexer and parser, and want to recompile, you can use `make`. This assumes that you have installed Bison > 2.6.
* If you want to run the test preparation scripts, you need ete3. (`conda install -c etetoolkit ete3`)


## To discuss:

* look through abort-- how to handle? also look for cassert
* unique_ptr? Are we copying things in parsing?
* "NodeId" type rather than unsigned int?
* Bitset is limited to `size_t`, but intpacking is in terms of int32...
* Conflict between `python_convention` and `CPPConvention`.


## Terminology

PCSS stands for parent-child subsplit.
They are represented as bitsets in three "chunks", which are sub-bit-sets.
For example, `100011001` is composed of the chunks `100`, `011` and `001`.
If the taxa are x0, x1, and x2 then this means the parent subsplit is (A, BC), and the child subsplit is (B,C).

* The first chunk is called the "uncut parent" because it is not further split apart by the child subsplit.
* The second chunk is called the "cut parent" because it is further split apart by the child subsplit.
* The third chunk is called the "child," and it's well defined relative to the cut parent: the other part of the subsplit is the cut parent setminus the child.
