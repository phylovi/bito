// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "libsbn.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <random>
#include <string>

namespace py = pybind11;

std::vector<double> make_vector() {
  std::vector<double> x(5);
  for (size_t i = 0; i < x.size(); i++) {
    x[i] = 0;
  }
  return x;
}

double total_vector(std::vector<double> v) {
  double total = 0;
  for (const auto &x : v) {
    total += x;
  }
  return total;
}

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
      .def("tree_log_likelihoods", &SBNInstance::TreeLogLikelihoods)
      .def("build_indexer", &SBNInstance::BuildIndexer)
      .def("sbn_total_prob", &SBNInstance::SBNTotalProb)
      .def_readwrite("sbn_probs", &SBNInstance::sbn_probs_);
  py::class_<std::vector<double>>(m, "vector_double", py::buffer_protocol())
      .def_buffer([](std::vector<double> &v) -> py::buffer_info {
        return py::buffer_info(
            v.data(),                                /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python
                                                       struct-style format
                                                       descriptor */
            1,                                       /* Number of dimensions */
            {v.size()},                              /* Buffer dimensions */
            {sizeof(double)});
      });
  m.def("make_vector", &make_vector, "test");
  m.def("total_vector", &total_vector, "test");
}
