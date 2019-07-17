// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "libsbn.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

namespace py = pybind11;

class Matrix {
 public:
  Matrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols) {
    m_data = new double[rows * cols];
  }
  double *data() { return m_data; }
  size_t rows() const { return m_rows; }
  size_t cols() const { return m_cols; }

 private:
  size_t m_rows, m_cols;
  double *m_data;
};

Matrix make_matrix() {
  Matrix x(4, 3);
  return x;
}

double get00(Matrix m) { return m.data()[0]; }

PYBIND11_MODULE(sbn, m) {
  m.doc() = "libsbn bindings";
  py::class_<SBNInstance>(m, "instance")
      .def(py::init<const std::string &>())
      .def("tree_count", &SBNInstance::TreeCount)
      .def("read_newick_file", &SBNInstance::ReadNewickFile)
      .def("read_nexus_file", &SBNInstance::ReadNexusFile)
      .def("read_fasta_file", &SBNInstance::ReadFastaFile)
      .def("print_status", &SBNInstance::PrintStatus)
      .def("split_supports", &SBNInstance::SplitSupports)
      .def("make_beagle_instances", &SBNInstance::MakeBeagleInstances)
      .def("tree_log_likelihoods", &SBNInstance::TreeLogLikelihoods);
  m.def("f", &SBNInstance::f, "test");
  py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
      .def_buffer([](Matrix &m) -> py::buffer_info {
        return py::buffer_info(
            m.data(),                                /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python struct-style
                                                       format descriptor */
            2,                                       /* Number of dimensions */
            {m.rows(), m.cols()},                    /* Buffer dimensions */
            {sizeof(double) * m.rows(), /* Strides (in bytes) for each index */
             sizeof(double)});
      });
  m.def("make_matrix", &make_matrix, "test");
  m.def("get00", &get00, "test");
}

