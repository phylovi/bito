#include "libsbn.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

// PYBIND11_MAKE_OPAQUE(Node::NodePtr);
// PYBIND11_MAKE_OPAQUE(Node::NodePtrVec);
// PYBIND11_MAKE_OPAQUE(Node::NodePtrVecPtr);

PYBIND11_MODULE(sbn, m) {
  m.doc() = "libsbn bindings";
  py::class_<SBNInstance>(m, "instance")
      .def(py::init<const std::string &>())
      //.def_readonly("trees", &SBNInstance::trees_)
      .def("tree_count", &SBNInstance::TreeCount)
      .def("parse_file", &SBNInstance::ParseFile)
      .def("print_status", &SBNInstance::PrintStatus)
      .def("g", &SBNInstance::g)
      .def("rootsplits", &SBNInstance::Rootsplits);
  m.def("f", &SBNInstance::f, "test");
}
