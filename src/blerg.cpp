// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include "Eigen/Dense"

namespace py = pybind11;

struct Blerg {
  Eigen::Matrix4d matrix;

  Eigen::Ref<Eigen::Matrix4d> GetMatrix() { return matrix; }
};

PYBIND11_MODULE(blerg, m) {
  m.doc() = "blerg bindings";
  py::class_<Blerg>(m, "Blerg")
      .def(py::init<>())
      .def("get_matrix", &Blerg::GetMatrix)
      .def_readwrite("matrix", &Blerg::matrix);
}
