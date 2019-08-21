# libsbn

[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/phylovi/libsbn.svg)](https://cloud.docker.com/u/phylovi/repository/docker/phylovi/libsbn/general) &nbsp;
[![Travis CI status](https://travis-ci.org/phylovi/libsbn.svg?branch=master)](https://travis-ci.org/phylovi/libsbn)

We are building a Python-interface C++ library so that you can express interesting parts of your phylogenetic model in Python/TensorFlow/PyTorch/etc and let libsbn handle the tree structure and likelihood computations for you.


## Dependencies

To install dependencies, use the associated conda environment file:

```
conda env create -f environment.yml
conda activate libsbn
```

**However, you also need to install platform-specific compiler packages as follows.**

* if you are on linux, use `conda install -y gxx_linux-64`
* if you are on OS X, use `conda install -y clangxx_osx-64`


## Building

`make` will build and run tests.

On OS X the build process will also modify the conda environment to point `DYLD_LIBRARY_PATH` to where BEAGLE is installed.
If you get an error about missing BEAGLE, just `conda activate libsbn` again and you should be good.

* (Optional) If you modify the lexer and parser, call `make bison`. This assumes that you have installed Bison > 3.4 (`conda install -c conda-forge bison`).
* (Optional) If you modify the test preparation scripts, call `make prep`. This assumes that you have installed ete3 (`conda install -c etetoolkit ete3`).


## Understanding

The following two papers will explain what this repository is about:

* Zhang & Matsen IV, NeurIPS 2018. [_Generalizing Tree Probability Estimation via Bayesian Networks_](http://papers.nips.cc/paper/7418-generalizing-tree-probability-estimation-via-bayesian-networks.pdf). Also see [blog post](https://matsen.fredhutch.org/general/2018/12/05/sbn.html).
* Zhang & Matsen IV, ICLR 2018. [_Variational Bayesian Phylogenetic Inference_](https://openreview.net/pdf?id=SJVmjjR9FX_).

In the off chance that you are citing this library, don't forget to cite the [BEAGLE paper](http://dx.doi.org/10.1093/sysbio/syz020) too, as we use BEAGLE!


## Contributing

libsbn is written in C++14.

### Style

We want the code to be:

1. correct, so we write tests
1. efficient in an algorithmic sense, so we consider algorithms carefully
1. clear to read and understand, so we write code with readers in mind and use code standards
1. fast, so we do profiling to find and eliminate bottlenecks
1. robust, so we use immutable data structures and safe C++ practices
1. simple and beautiful, so we keep the code as minimal and DRY as we can without letting it get convoluted or over-technical

Also let's:

* Prefer a functional style: returning variables versus modifying them in place. Because of return value optimization, this doesn't have a performance penalty.
* [RAII](https://en.cppreference.com/w/cpp/language/raii). No `new`.
* Avoid classic/raw pointers except as const parameters to functions.
* Prefer [variable names and simple coding practices](https://blog.codinghorror.com/coding-without-comments/) to code comments.
  If that means having long identifier names, that's fine!
  If you can't make the code use and operation inherently obvious, please write documentation.
* Prefer GitHub issues to TODO comments in code.
* Always use curly braces for the body of conditionals and loops, even if they are one line.
* Prefer types of known size, such as `uint32_t`, to types that can vary across architectures.

The [C++ core guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines) are the authority for how to write C++, and we will follow them.
For issues not covered by these guidelines (especially naming conventions), we will use the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html) to the letter.
We use [cpplint](https://github.com/cpplint/cpplint) to check some aspects of this.

There are certainly violations of these guidelines in the code, so fix them when you see them!


### Formatting

Code gets formatted using [clang-format](https://clang.llvm.org/docs/ClangFormat.html).
See the Makefile for the invocation.


### Tests

Add a test for every new feature.

* For C++, we use [doctest](https://github.com/onqtam/doctest).
* For Python, we use [pytest](https://docs.pytest.org/en/latest/).


### [Git flow](https://guides.github.com/introduction/flow/)

* Code changes start by raising an issue proposing the changes, which often leads to a discussion
* Make a branch associated with the issue named with the issue number and a description, such as `4-efficiency-improvements` for a branch associated with issue #4 about efficiency improvements
* If you have another branch to push for the same issue (perhaps a fresh, alternate start), you can just name them consecutively `4-1-blah`, `4-2-etc`, and so on
* Push code to that branch
* Once the code is ready to merge, open a [pull request](https://help.github.com/articles/using-pull-requests/)
* Code review on GitHub
* [Squash and merge](https://help.github.com/en/articles/merging-a-pull-request), [closing the issue via the squash and merge commit message](https://help.github.com/articles/closing-issues-via-commit-messages/)
* Delete branch


## Terminology

* PCSS stands for parent-child subsplit.
  It's a general concept rather than a specific implementation of the concept.
  For example, see the documentation of PCSSFun (in `node.hpp`) and PCSS Bitsets (in `bitset.hpp`) for two different ways of using this concept.


## Contributors

* Erick Matsen ([@matsen](https://github.com/matsen)): coding, design, maintenance
* Cheng Zhang ([@zcrabbit](https://github.com/zcrabbit)): concept, design, algorithms
* Mathieu Fourment ([@4ment](https://github.com/4ment)): design, site pattern compression, BEAGLE gradients and scaling
* Christiaan Swanepoel ([@christiaanjs](https://github.com/christiaanjs)): design
