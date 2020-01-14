libsbn concepts
===============

Definitions
-----------

:index:`PCSS`
  PCSS stands for parent-child subsplit pair, that is, a pair consisting of a parent subsplit and a child subsplit.
  It's a general concept rather than a specific implementation of the concept.
  For example, see the documentation of PCSSFun (in ``node.hpp``) and PCSS Bitsets (in ``bitset.hpp``) for two different ways of using this concept.

:index:`PSP`
  PSP stands for primary subsplit pair, which is a way of parameterizing branch lengths on trees.
  See the 2019 ICLR paper for details.

:index:`SBN`
  PSP stands for subsplit Bayes network, which is how we parameterize variational distributions on tree topologies.
  See the 2018 NeurIPS paper for details.


Terminology and conventions
---------------------------

We follow different naming conventions for C++ and Python.
In C++ we follow the Google C++ conventions: ``PascalCase`` for types and functions, ``snake_case`` for variables, and ``a_trailing_underscore_`` for member variables.
In Python we just use ``snake_case`` for everything, with no trailing underscore.
Here we'll just use the Python versions, and the C++ hackers will have to extrapolate.

In Python bitsets will be expressed as strings, while in C++ they are Bitset objects.
Hashtables are implemented using ``unordered_map`` in C++ and dictionaries in Python, but we will use the Python terminology even if we are talking about a C++ map.


Implementation
--------------
All of the SBN parameters are held in a single ``sbn_parameters`` block of doubles that is available for reading and writing through the libsbn instance interface.
We index into that block using an "indexer" dictionary which maps from bitsets representing PCSSs to the index for that PCSS.
See the definition of PCSS bitset in ``bitset.hpp`` and the tests in ``libsbn.hpp`` to see this in action.

The ``sbn_parameters`` block is set up such that PCSSs that share a parent are laid out contiguously.
This is nice because if we want to sample from the child subsplit conditional on the parent, then we have all of the probabilities there in one place.
The ``parent_to_range`` dictionary maps parent bitsets to the range of indices corresponding to children descending from that parent.

We can ask for an "indexer representation" of a tree, which is described in ``sbn_maps.hpp`` and shown off in the use of ``StringIndexerRepresentationOf`` in ``libsbn.hpp``.
Basically, its the expression of the rootsplits and PCSSs for a tree, in index form.

PSPs are indexed using a different scheme, which uses a vector of three vectors.
The first vector describes the rootsplit, the second describes the child subsplit going "up" the tree, and the third the child subsplit going "down" the tree.
These directionality terms only make sense on rooted trees, but we represent unrooted trees as rooted trees in the typical way.
