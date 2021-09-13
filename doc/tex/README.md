# Implementation notes

This tex document contains a full derivation of the algorithms implemented in the code.
Thus is only of interest for people wishing to understand the inner workings of the algorithms, such as to read and/or work on the code.

For API documentation see [phylovi.github.io/bito/](https://phylovi.github.io/bito/).


## Compilation

You need a tex distribution.
Type `scons --sconstruct SConstruct-noinkscape` in the command line if you don't have Inkscape all set up to work on the command line.
If you do, `scons` should work.
Use `scons continuous` to have continuous building after installing inotify.

