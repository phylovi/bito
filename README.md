# libsbn

[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/phylovi/libsbn.svg)](https://hub.docker.com/r/phylovi/libsbn) &nbsp;
[![Travis CI status](https://travis-ci.org/phylovi/libsbn.svg?branch=master)](https://travis-ci.org/phylovi/libsbn)

We are building a Python-interface C++ library for phylogenetic variational inference so that you can express interesting parts of your phylogenetic model in Python/TensorFlow/PyTorch/etc and let libsbn handle the tree structure and likelihood computations for you.


## Dependencies

libsbn is written in C++17.
The associated Python module, `vip`, is targeting Python 3.7.

* If you are on linux, install gcc >= 8, which is standard in Debian Buster and Ubuntu 18.04
* If you are on OS X, use a recent version of Xcode and install command line tools

Then, install **the `hmc-clock` branch** of [BEAGLE](https://github.com/beagle-dev/beagle-lib).
This will require a from-source installation, as in their docs, but you have to do a full `git clone` (no `--depth=1`).
You can see a full installation procedure by taking a look at the [conda-beagle Dockerfile](https://github.com/matsengrp/conda-beagle/blob/master/Dockerfile).

To install additional dependencies, use the associated conda environment file:

    conda env create -f environment.yml
    conda activate libsbn

If you want to specify your compiler manually, set the `CC` and `CXX` shell variables to your desired compiler command.

The notebooks require R, IRKernel, rpy2 >=3.1.0, and some R packages such as ggplot and cowplot.
Do not install R via conda.
Doing so will install the conda compiler toolchain, this will mess up our compilation.


## Building

For your first build, do

* `git submodule update --init --recursive`
* `scons`
* Respond to interactive prompts about where `hmc-clock` BEAGLE is installed
* `conda activate libsbn`
* `make`

After these steps `make` will build, run tests, and install the Python packages, and this should be the only command you need to run after modifying the code.

The build process will modify the conda environment to point `[DY]LD_LIBRARY_PATH` to where BEAGLE is installed.
If you get an error about missing BEAGLE, just `conda activate libsbn` again and you should be good.
If you want to modify your desired BEAGLE installation location, do `unset BEAGLE_PREFIX` and start the steps above again starting at `scons`.

* (Optional) If you modify the lexer and parser, call `make bison`. This assumes that you have installed Bison > 3.4 (`conda install -c conda-forge bison`).
* (Optional) If you modify the test preparation scripts, call `make prep`. This assumes that you have installed ete3 (`conda install -c etetoolkit ete3`).


## Understanding

The following two papers will explain what this repository is about:

* Zhang & Matsen IV, NeurIPS 2018. [_Generalizing Tree Probability Estimation via Bayesian Networks_](http://papers.nips.cc/paper/7418-generalizing-tree-probability-estimation-via-bayesian-networks.pdf); 👉🏽 [blog post](https://matsen.fredhutch.org/general/2018/12/05/sbn.html).
* Zhang & Matsen IV, ICLR 2019. [_Variational Bayesian Phylogenetic Inference_](https://openreview.net/pdf?id=SJVmjjR9FX_); 👉🏽 [blog post](https://matsen.fredhutch.org/general/2019/08/24/vbpi.html).

Our documentation consists of:

* [Online documentation](https://phylovi.github.io/libsbn/)
* Derivations in `doc/tex`, which explain what's going on in the code.


## Contributing

We welcome your contributions!
Please see our detailed [contributing guidelines](CONTRIBUTING.md).


## Contributors

* Erick Matsen ([@matsen](https://github.com/matsen)): implementation, design, janitorial duties
* Mathieu Fourment ([@4ment](https://github.com/4ment)): implementation of substitution models and likelihoods/gradients, design
* Seong-Hwan Jun ([@junseonghwan](https://github.com/junseonghwan)): implementation of SBN gradients, design
* Cheng Zhang ([@zcrabbit](https://github.com/zcrabbit)): concept, design, algorithms
* Christiaan Swanepoel ([@christiaanjs](https://github.com/christiaanjs)): design
* Xiang Ji ([@xji3](https://github.com/xji3)): gradient expertise and node height code
* Marc Suchard ([@msuchard](https://github.com/msuchard)): gradient expertise and node height code
* Michael Karcher ([@mdkarcher](https://github.com/mdkarcher)): SBN expertise
* Ognian Milanov ([@ognian-](https://github.com/ognian-)): C++ wisdom and optimization


## Citations

If you are citing this library, please cite the NeurIPS and ICLR papers listed above.
We require BEAGLE, so please also cite these papers:

* [BEAGLE 3](http://dx.doi.org/10.1093/sysbio/syz020)
* [Linear-time, high-dimensional gradient in BEAGLE](http://arxiv.org/abs/1905.12146)


## Acknowledgements

* Jaime Huerta-Cepas: several tree traversal functions are copied from [ete3](https://github.com/etetoolkit/ete)
* Thomas Junier: parts of the parser are copied from [newick\_utils](https://github.com/tjunier/newick_utils)
* The parser driver is derived from the [Bison C++ example](https://www.gnu.org/software/bison/manual/html_node/Calc_002b_002b-Parsing-Driver.html#Calc_002b_002b-Parsing-Driver)

In addition to the packages mentioned above we also employ:

* [cxx-prettyprint](https://github.com/louisdx/cxx-prettyprint) STL container pretty printing
* [Eigen](https://gitlab.com/libeigen/eigen)
* [Progress-CPP](https://github.com/prakhar1989/progress-cpp) progress bar
