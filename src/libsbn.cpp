// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "libsbn.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

namespace py = pybind11;

// In order to make vector<double>s available to numpy, we take two steps.
// First, we make them opaque to pybind11, so that it doesn't do its default
// conversion of STL types.
PYBIND11_MAKE_OPAQUE(std::vector<double>);

PYBIND11_MODULE(sbn, m) {
  m.doc() = "libsbn bindings";
  // Second, we expose vector<double> as a buffer object so that we can use it
  // as an in-place numpy array with np.array(v, copy=False). See
  // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
  py::class_<std::vector<double>>(m, "vector_double", py::buffer_protocol())
      .def_buffer([](std::vector<double> &v) -> py::buffer_info {
        return py::buffer_info(
            v.data(),                                 // Pointer to buffer
            sizeof(double),                           // Size of one scalar
            py::format_descriptor<double>::format(),  // See docs
            1,                                        // Number of dimensions
            {v.size()},                               // Buffer dimensions
            {sizeof(double)});                        // Stride
      });
  // Expose trees.
  py::class_<Tree>(m, "Tree", py::buffer_protocol())
      .def_static("of_index_vector", &Tree::OfIndexVector)
      .def_buffer([](Tree &tree) -> py::buffer_info {
        return py::buffer_info(
            tree.branch_lengths_.data(),              // Pointer to buffer
            sizeof(double),                           // Size of one scalar
            py::format_descriptor<double>::format(),  // See docs
            1,                                        // Number of dimensions
            {tree.branch_lengths_.size()},            // Buffer dimensions
            {sizeof(double)});                        // Stride
      });

  // Now we set things up our SBNInstance class.
  py::class_<SBNInstance>(m, "instance")
      // Constructors
      .def(py::init<const std::string &>())
      // Methods
      .def("tree_count", &SBNInstance::TreeCount)
      .def("read_newick_file", &SBNInstance::ReadNewickFile)
      .def("read_nexus_file", &SBNInstance::ReadNexusFile)
      .def("read_fasta_file", &SBNInstance::ReadFastaFile)
      .def("print_status", &SBNInstance::PrintStatus)
      .def("split_counters", &SBNInstance::SplitCounters)
      .def("make_beagle_instances", &SBNInstance::MakeBeagleInstances)
      .def("log_likelihoods", &SBNInstance::LogLikelihoods)
      .def("branch_gradients", &SBNInstance::BranchGradients)
      .def("process_loaded_trees", &SBNInstance::ProcessLoadedTrees)
      .def("get_indexers", &SBNInstance::GetIndexers)
      .def("sbn_total_prob", &SBNInstance::SBNTotalProb)
      // Member Variables
      .def_readwrite("sbn_probs", &SBNInstance::sbn_probs_);
}
