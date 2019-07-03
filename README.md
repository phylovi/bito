# libsbn

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

* If you have modified the lexer and parser, use `make bison`. This assumes that you have installed Bison > 2.6.
* If you want to run the test preparation scripts, use `make prep`. This assumes that you have installed ete3 (`conda install -c etetoolkit ete3`).


## Terminology

PCSS stands for parent-child subsplit.
They are represented as bitsets in three equal-sized "chunks", which are sub-bit-sets.
For example, `100011001` is composed of the chunks `100`, `011` and `001`.
If the taxa are x0, x1, and x2 then this means the parent subsplit is (A, BC), and the child subsplit is (B,C).

* The first chunk is called the "uncut parent" because it is not further split apart by the child subsplit.
* The second chunk is called the "cut parent" because it is further split apart by the child subsplit.
* The third chunk is called the "child," and it's well defined relative to the cut parent: the other part of the subsplit is the cut parent setminus the child.
