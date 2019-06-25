#include <pybind11/pybind11.h>
#include "libsbn.hpp"

namespace py = pybind11;

PYBIND11_MODULE(sbn, m) {
  py::class_<SBNInstance>(m, "instance")
      .def(py::init<const std::string &>())
      .def("parse_file", &SBNInstance::ParseFile)
      .def("print_status", &SBNInstance::PrintStatus);
}
