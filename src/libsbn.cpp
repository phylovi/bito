// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "libsbn.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

namespace py = pybind11;

// PYBIND11_MAKE_OPAQUE(Node::NodePtr);
// PYBIND11_MAKE_OPAQUE(Node::NodePtrVec);
// PYBIND11_MAKE_OPAQUE(Node::NodePtrVecPtr);

PYBIND11_MODULE(sbn, m) {
  m.doc() = "libsbn bindings";
  py::class_<SBNInstance>(m, "instance")
      .def(py::init<const std::string &>())
      // .def_readonly("trees", &SBNInstance::trees_)
      .def("tree_count", &SBNInstance::TreeCount)
      .def("parse_file", &SBNInstance::ParseFile)
      .def("read_fasta", &SBNInstance::ReadFasta)
      .def("print_status", &SBNInstance::PrintStatus)
      .def("rootsplit_support", &SBNInstance::RootsplitSupport)
      .def("subsplit_support", &SBNInstance::SubsplitSupport);
  m.def("f", &SBNInstance::f, "test");
}
