// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include <pybind11/eigen.h>
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
  m.doc() = "Python interface to libsbn.";
  // Second, we expose them as buffer objects so that we can use them
  // as in-place numpy arrays with np.array(v, copy=False). See
  // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
  py::class_<std::vector<double>>(m, "vector_double",
                                  "A wrapper for vector<double>.",
                                  py::buffer_protocol())
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
  py::class_<Tree>(m, "Tree", "A tree with branch lengths.",
                   py::buffer_protocol())
      .def("parent_id_vector", &Tree::ParentIdVector)
      .def_static("of_parent_id_vector", &Tree::OfParentIdVector)
      .def_readwrite("branch_lengths", &Tree::branch_lengths_);
  // TreeCollection
  py::class_<TreeCollection>(m, "TreeCollection", R"raw(
  A collection of trees.

  In addition to the methods, TreeCollection also offers direct access to
  the trees through the ``trees`` member variable.
  )raw")
      .def(py::init<Tree::TreeVector>())
      .def(py::init<Tree::TreeVector, TagStringMap>())
      .def(py::init<Tree::TreeVector, const std::vector<std::string> &>())
      .def("erase", &TreeCollection::Erase,
           "Erase the specified range from the current tree collection.")
      .def("newick", &TreeCollection::Newick,
           "Get the current set of trees as a big Newick string.")
      .def_readwrite("trees", &TreeCollection::trees_);
  // PSPIndexer
  py::class_<PSPIndexer>(m, "PSPIndexer", "The primary split pair indexer.")
      .def("details", &PSPIndexer::Details);
  // PhyloModelSpecification
  py::class_<PhyloModelSpecification>(m, "PhyloModelSpecification",
                                      R"raw(
    Phylogenetic model specification.

    This is how we specify phylogenetic models, with strings for the substitution
    model, the site model, and the clock model.
      )raw")
      .def(py::init<const std::string &, const std::string &,
                    const std::string &>(),
           py::arg("substitution"), py::arg("site"), py::arg("clock"));
  // SBNInstance
  py::class_<SBNInstance>(m, "instance",
                          "A wrapper for the all of the C++-side state.")
      // ** Initialization and status
      .def(py::init<const std::string &>())
      .def("tree_count", &SBNInstance::TreeCount,
           "Return the number of trees that are currently stored in the "
           "instance.")
      .def("print_status", &SBNInstance::PrintStatus,
           "Print information about the instance.")
      // ** SBN-related items
      .def("process_loaded_trees", &SBNInstance::ProcessLoadedTrees, R"raw(
          Process the trees currently stored in the instance.

          Specifically, parse them and build the indexers and the ``sbn_parameters`` vector.
      )raw")
      .def("get_indexers", &SBNInstance::GetIndexers,
           "Return the indexer and parent_to_range as string-keyed maps.")
      .def("sample_trees", &SBNInstance::SampleTrees,
           "Sample trees from the SBN and store them internally.",
           py::arg("count"))
      .def("get_indexer_representations",
           &SBNInstance::GetIndexerRepresentations)
      .def("get_psp_indexer_representations",
           &SBNInstance::GetPSPIndexerRepresentations)
      .def("split_counters", &SBNInstance::SplitCounters,
           "A testing method to count splits.")
      .def("split_lengths", &SBNInstance::SplitLengths)
      // ** Phylogenetic likelihood
      .def("prepare_for_phylo_likelihood",
           &SBNInstance::PrepareForPhyloLikelihood,
           "Prepare instance for phylogenetic likelihood computation.",
           py::arg("specification"), py::arg("thread_count"),
           py::arg("tree_count_option") = std::nullopt)
      .def("get_phylo_model_params", &SBNInstance::GetPhyloModelParams)
      .def("get_phylo_model_param_block_map",
           &SBNInstance::GetPhyloModelParamBlockMap)
      .def("resize_phylo_model_params", &SBNInstance::ResizePhyloModelParams,
           "Resize phylo_model_params.",
           py::arg("tree_count_option") = std::nullopt)
      .def("log_likelihoods", &SBNInstance::LogLikelihoods)
      .def("branch_gradients", &SBNInstance::BranchGradients)
      // ** I/O
      .def("read_newick_file", &SBNInstance::ReadNewickFile,
           "Read trees from a Newick file.")
      .def("read_nexus_file", &SBNInstance::ReadNexusFile,
           "Read trees from a Nexus file.")
      .def("read_fasta_file", &SBNInstance::ReadFastaFile,
           "Read a sequence alignment from a FASTA file.")
      // Member Variables
      .def_readonly("psp_indexer", &SBNInstance::psp_indexer_)
      .def_readwrite("sbn_parameters", &SBNInstance::sbn_parameters_)
      .def_readwrite("tree_collection", &SBNInstance::tree_collection_);
  // If you want to be sure to get all of the stdout and cerr messages, put your
  // Python code in a context like so:
  // `with libsbn.ostream_redirect(stdout=True, stderr=True):`
  // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/utilities.html#capturing-standard-output-from-ostream
  py::add_ostream_redirect(m, "ostream_redirect");
}
