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
