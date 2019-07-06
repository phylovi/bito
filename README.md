# libsbn

_This package is not yet useful. We're working hard to make it useful._

[![Travis CI status](https://travis-ci.org/matsengrp/libsbn.svg?branch=master)](https://travis-ci.org/matsengrp/libsbn)

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

Just use `make` to build and run tests.

* If you have modified the lexer and parser, use `make bison`. This assumes that you have installed Bison > 2.6 (`conda install -c conda-forge bison`).
* If you want to run the test preparation scripts, use `make prep`. This assumes that you have installed ete3 (`conda install -c etetoolkit ete3`).


## Contributing

libsbn is written in C++17.

### Style

* Prefer a functional style: returning variables versus modifying them in place. Because of return value optimization, this doesn't have a performance penalty.
* RAII. No "new."

We will use the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html) to the letter.
I use [cpplint](https://github.com/cpplint/cpplint) to check some aspects of this.

Notes:

* Always use curly braces for the body of conditionals and loops, even if they are one line.


### Formatting

Code gets formatted using [clang-format](https://clang.llvm.org/docs/ClangFormat.html).
See the Makefile for the invocation.


### Tests

Add a test for every new feature.

* For C++, use [doctest](https://github.com/onqtam/doctest)
* For Python, use [pytest](https://docs.pytest.org/en/latest/)


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

PCSS stands for parent-child subsplit.
They are represented as bitsets in three equal-sized "chunks", which are sub-bit-sets.
For example, `100011001` is composed of the chunks `100`, `011` and `001`.
If the taxa are x0, x1, and x2 then this means the parent subsplit is (A, BC), and the child subsplit is (B,C).

* The first chunk is called the "uncut parent" because it is not further split apart by the child subsplit.
* The second chunk is called the "cut parent" because it is further split apart by the child subsplit.
* The third chunk is called the "child," and it's well defined relative to the cut parent: the other part of the subsplit is the cut parent setminus the child.
