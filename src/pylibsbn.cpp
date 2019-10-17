// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include "libsbn.hpp"

namespace py = pybind11;

// In order to make vector<double>s available to numpy, we take two steps.
// First, we make them opaque to pybind11, so that it doesn't do its default
// conversion of STL types.
PYBIND11_MAKE_OPAQUE(std::vector<double>);

PYBIND11_MODULE(libsbn, m) {
  m.doc() = "libsbn bindings";
  // Second, we expose them as buffer objects so that we can use them
  // as in-place numpy arrays with np.array(v, copy=False). See
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
  // Tree
  py::class_<Tree>(m, "Tree", py::buffer_protocol())
      .def("parent_id_vector", &Tree::ParentIdVector)
      .def_static("of_parent_id_vector", &Tree::OfParentIdVector)
      .def_readwrite("branch_lengths", &Tree::branch_lengths_);
  // TreeCollection
  py::class_<TreeCollection>(m, "TreeCollection")
      .def(py::init<Tree::TreeVector>())
      .def(py::init<Tree::TreeVector, TagStringMap>())
      .def(py::init<Tree::TreeVector, const std::vector<std::string> &>())
      .def("newick", &TreeCollection::Newick)
      .def_readwrite("trees", &TreeCollection::trees_);
  // PSPIndexer
  py::class_<PSPIndexer>(m, "PSPIndexer").def("details", &PSPIndexer::Details);
  // SBNInstance
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
      .def("parameter_gradients", &SBNInstance::ParameterGradients)
      .def("process_loaded_trees", &SBNInstance::ProcessLoadedTrees)
      .def("get_indexers", &SBNInstance::GetIndexers)
      .def("sample_trees", &SBNInstance::SampleTrees)
      .def("get_indexer_representations",
           &SBNInstance::GetIndexerRepresentations)
      .def("get_psp_indexer_representations",
           &SBNInstance::GetPSPIndexerRepresentations)
      .def("split_lengths", &SBNInstance::SplitLengths)
      .def("set_beagle_subst_models", &SBNInstance::SetBeagleSubstModels)
      // Member Variables
      .def_readonly("psp_indexer", &SBNInstance::psp_indexer_)
      .def_readwrite("sbn_parameters", &SBNInstance::sbn_parameters_)
      .def_readwrite("tree_collection", &SBNInstance::tree_collection_)
      // TODO
      .def_readwrite("evec", &SBNInstance::evec)
      .def_readwrite("ivec", &SBNInstance::ivec)
      .def_readwrite("eval", &SBNInstance::eval)
      .def_readwrite("freqs", &SBNInstance::freqs)
      .def_readwrite("q_differential", &SBNInstance::q_differential)
    ;
  // If you want to be sure to get all of the stdout and cerr messages, put your
  // Python code in a context like so:
  // `with libsbn.ostream_redirect(stdout=True, stderr=True):`
  // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/utilities.html#capturing-standard-output-from-ostream
  py::add_ostream_redirect(m, "ostream_redirect");
}
