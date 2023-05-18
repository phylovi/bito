# bito

[![Build with cmake](https://github.com/phylovi/bito/actions/workflows/build-with-cmake.yml/badge.svg)](https://github.com/phylovi/bito/actions/workflows/build-with-cmake.yml)

`bito`, or "Bayesian Inference of Trees via Optimization", is a Python-interface C++ library for phylogenetic variational inference so that you can express interesting parts of your phylogenetic model in Python/TensorFlow/PyTorch/etc and let bito handle the tree structure and likelihood computations for you.
"Bito" is also the name of a [tree](https://www.merriam-webster.com/dictionary/bito) native to Africa that produces medicinal oil.
We pronounce "bito" with a long /e/ sound ("bito" rhymes with "burrito").

This library is in an experimental state.
This library was formerly known as "libsbn".


## Dependencies

* If you are on linux, install gcc >= 7.5, which is standard in Debian Buster and Ubuntu 18.04
* If you are on OS X, use a recent version of Xcode and install command line tools

We suggest using [anaconda](https://docs.conda.io/en/latest/miniconda.html) and the associated conda environment file, which will nicely install relevant dependencies:

    conda env create -f environment.yml
    conda activate bito

(Very optional) The notebooks require R, IRKernel, rpy2 >=3.1.0, and some R packages such as ggplot and cowplot.


## Building

For your first build, do

* `git submodule update --init --recursive`
* `make`

This will install the `bito` Python module.

You can build and run tests using `make test` and `make fasttest` (the latter excludes some slow tests).

Note that `make` accepts `-j` flags for multi-core builds: e.g. `-j20` will build with 20 jobs.

* (Optional) If you modify the lexer and parser, call `make bison`. This assumes that you have installed Bison >= 3.4 (`conda install -c conda-forge bison`).
* (Optional) If you modify the test preparation scripts, call `make prep`. This assumes that you have installed ete3 (`conda install -c etetoolkit ete3`).


## Understanding

The following two papers will explain what this repository is about:

* Zhang & Matsen IV, NeurIPS 2018. [_Generalizing Tree Probability Estimation via Bayesian Networks_](http://papers.nips.cc/paper/7418-generalizing-tree-probability-estimation-via-bayesian-networks.pdf); 👉🏽 [blog post](https://matsen.fredhutch.org/general/2018/12/05/sbn.html).
* Zhang & Matsen IV, ICLR 2019. [_Variational Bayesian Phylogenetic Inference_](https://openreview.net/pdf?id=SJVmjjR9FX_); 👉🏽 [blog post](https://matsen.fredhutch.org/general/2019/08/24/vbpi.html).

Our documentation consists of:

* [Online documentation](https://phylovi.github.io/bito/)
* Derivations in `doc/tex`, which explain what's going on in the code.


## Contributing

We welcome your contributions!
Please see our detailed [contribution guidelines](CONTRIBUTING.md).


## Contributors

* Erick Matsen ([@matsen](https://github.com/matsen)): implementation, design, janitorial duties
* Dave H. Rich ([@DaveRich](https://github.com/davidrich27)): core developer
* Ognian Milanov ([@ognian-](https://github.com/ognian-)): core developer
* Mathieu Fourment ([@4ment](https://github.com/4ment)): implementation of substitution models and likelihoods/gradients, design
* Seong-Hwan Jun ([@junseonghwan](https://github.com/junseonghwan)): generalized pruning design and implementation, implementation of SBN gradients, design
* Hassan Nasif ([@hrnasif](https://github.com/hrnasif)): hot start for generalized pruning; gradient descent for generalized pruning
* Anna Kooperberg ([@annakooperberg](https://github.com/annakooperberg)): refactoring the subsplit DAG
* Sho Kiami ([@shokiami](https://github.com/shokiami)): refactoring the subsplit DAG
* Tanvi Ganapathy ([@tanviganapathy](https://github.com/tanviganapathy)): refactoring the subsplit DAG
* Lucy Yang ([@lucyyang01](https://github.com/lucyyang01)): subsplit DAG visualization
* Cheng Zhang ([@zcrabbit](https://github.com/zcrabbit)): concept, design, algorithms
* Christiaan Swanepoel ([@christiaanjs](https://github.com/christiaanjs)): design
* Xiang Ji ([@xji3](https://github.com/xji3)): gradient expertise and node height code
* Marc Suchard ([@msuchard](https://github.com/msuchard)): gradient expertise and node height code
* Michael Karcher ([@mdkarcher](https://github.com/mdkarcher)): SBN expertise
* Eric J. Isaac ([@EricJIsaac](https://github.com/EricJIsaac)): C++ wisdom

## Citations

If you are citing this library, please cite the NeurIPS and ICLR papers listed above.
We require BEAGLE, so please also cite these papers:

* [BEAGLE 3](http://dx.doi.org/10.1093/sysbio/syz020)
* [Linear-time, high-dimensional gradient in BEAGLE](https://doi.org/10.1093/molbev/msaa130)


## Acknowledgements

* Jaime Huerta-Cepas: several tree traversal functions are copied from [ete3](https://github.com/etetoolkit/ete)
* Thomas Junier: parts of the parser are copied from [newick\_utils](https://github.com/tjunier/newick_utils)
* The parser driver is derived from the [Bison C++ example](https://www.gnu.org/software/bison/manual/html_node/Calc_002b_002b-Parsing-Driver.html#Calc_002b_002b-Parsing-Driver)

In addition to the packages mentioned above we also employ:

* [cxx-prettyprint](https://github.com/louisdx/cxx-prettyprint) STL container pretty printing
* [Eigen](https://gitlab.com/libeigen/eigen)
* [fast-cpp-csv-parser](https://github.com/ben-strasser/fast-cpp-csv-parser)
* [Progress-CPP](https://github.com/prakhar1989/progress-cpp) progress bar
