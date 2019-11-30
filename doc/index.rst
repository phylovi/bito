Welcome to libsbn's documentation!
==================================

.. toctree::
   :caption: Contents:

See the README on the `libsbn GitHub repository <https://github.com/phylovi/libsbn>`_ for installation and general description.
This is the documentation for interacting with libsbn from Python.

From the perspective of a Python user, code is broken into two components: ``libsbn``, which is the direct interface to the C++ code, and ``vip``, which implements the Python components.


Modules
=======

.. autosummary::
   :toctree: modules

   libsbn
   vip
   vip.benchmark
   vip.branch_model
   vip.burrito
   vip.cli
   vip.optimizers
   vip.priors
   vip.sbn_model
   vip.scalar_model
   vip.sgd_server

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
