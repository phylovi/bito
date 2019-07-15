// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "libsbn.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

namespace py = pybind11;

PYBIND11_MODULE(sbn, m) {
  m.doc() = "libsbn bindings";
  py::class_<SBNInstance>(m, "instance")
      .def(py::init<const std::string &>())
      .def("tree_count", &SBNInstance::TreeCount)
      .def("read_newick_file", &SBNInstance::ReadNewickFile)
      .def("read_fasta_file", &SBNInstance::ReadFastaFile)
      .def("print_status", &SBNInstance::PrintStatus)
      .def("split_supports", &SBNInstance::SplitSupports)
      .def("add_beagle_instances", &SBNInstance::AddBeagleInstances)
      .def("tree_log_likelihoods", &SBNInstance::TreeLogLikelihoods);
  m.def("f", &SBNInstance::f, "test");
}
