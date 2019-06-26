#include "libsbn.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(sbn, m) {
  m.doc() = "libsbn bindings";
  py::class_<SBNInstance>(m, "instance")
      .def(py::init<const std::string &>())
      .def_readwrite("indexer", &SBNInstance::indexer_)
      .def("tree_count", &SBNInstance::TreeCount)
      .def("parse_file", &SBNInstance::ParseFile)
      .def("init_indexer", &SBNInstance::InitIndexer)
      .def("print_status", &SBNInstance::PrintStatus);
  m.def("f", &SBNInstance::f, "test");
}
