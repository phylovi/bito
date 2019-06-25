#include <pybind11/pybind11.h>
#include "libsbn.hpp"

namespace py = pybind11;

PYBIND11_MODULE(sbn, m) {
  py::class_<SBNInstance>(m, "instance")
      .def(py::init<const std::string &>())
      .def("parse_file", &SBNInstance::ParseFile)
      .def("print_status", &SBNInstance::PrintStatus);
}


// int add(int i, int j) {
//   struct SBNInstance* inst = sbn_NewInstance();
//   sbn_PrintStatus(inst);
//   return i + j;
// }
//
// PYBIND11_MODULE(example, m) {
//     m.doc() = "pybind11 example plugin"; // optional module docstring
//
//     m.def("add", &add, "A function which adds two numbers");
// }
