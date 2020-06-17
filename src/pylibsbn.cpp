// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include "rooted_gradient_transforms.hpp"
#include "rooted_sbn_instance.hpp"
#include "unrooted_sbn_instance.hpp"

namespace py = pybind11;

// This is how we can have Eigen objects be directly mutable from Python. See
// https://github.com/eacousineau/repro/blob/f4ba595d077af7363f501f6c85d3d2449219f04a/python/pybind11/custom_tests/test_tmp.cc#L16-L38
// Thanks to @eacousineau!
template <typename PyClass, typename C, typename D>
void def_read_write_mutable(PyClass &cls, const char *name, D C::*pm) {
  cls.def_property(name, [pm](C & self) -> auto & { return self.*pm; },
                   [pm](C &self, const D &value) { self.*pm = value; });
}

// In order to make vector<double>s available to numpy, we take two steps.
// First, we make them opaque to pybind11, so that it doesn't do its default
// conversion of STL types.
PYBIND11_MAKE_OPAQUE(std::vector<double>);

// MODULE
PYBIND11_MODULE(libsbn, m) {
  m.doc() = R"raw(Python interface to libsbn.)raw";
  // Second, we expose them as buffer objects so that we can use them
  // as in-place numpy arrays with np.array(v, copy=False). See
  // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
  py::class_<std::vector<double>>(m, "vector_double", "A wrapper for vector<double>.",
                                  py::buffer_protocol())
      .def_buffer([](std::vector<double> &v) -> py::buffer_info {
        return py::buffer_info(v.data(),        // Pointer to buffer
                               sizeof(double),  // Size of one scalar
                               py::format_descriptor<double>::format(),  // See docs
                               1,                  // Number of dimensions
                               {v.size()},         // Buffer dimensions
                               {sizeof(double)});  // Stride
      });

  // CLASS
  // RootedTree
  py::class_<RootedTree>(m, "RootedTree", "A rooted tree with branch lengths.",
                         py::buffer_protocol())
      .def("parent_id_vector", &RootedTree::ParentIdVector)
      .def("set_node_heights_via_height_ratios",
           &RootedTree::SetNodeHeightsViaHeightRatios)
      .def_static("example", &RootedTree::Example)
      .def_static("of_parent_id_vector", &RootedTree::OfParentIdVector)
      .def_readwrite("branch_lengths", &RootedTree::branch_lengths_)
      .def_readwrite("height_ratios", &RootedTree::height_ratios_)
      .def_readwrite("node_heights", &RootedTree::node_heights_)
      .def_readwrite("node_bounds", &RootedTree::node_bounds_)
      .def_readwrite("rates", &RootedTree::rates_);

  // CLASS
  // RootedTreeCollection
  py::class_<RootedTreeCollection>(m, "RootedTreeCollection", R"raw(
  A collection of rooted trees.

  In addition to the methods, RootedTreeCollection also offers direct access to
  the trees through the ``trees`` member variable.
  )raw")
      .def(py::init<RootedTree::RootedTreeVector>(), "The empty constructor.")
      .def(py::init<RootedTree::RootedTreeVector, TagStringMap>(),
           "Constructor from a vector of trees and a tags->taxon names map.")
      .def(py::init<RootedTree::RootedTreeVector, const std::vector<std::string> &>(),
           "Constructor from a vector of trees and a vector of taxon names.")
      .def("erase", &RootedTreeCollection::Erase,
           "Erase the specified range from the current tree collection.")
      .def("drop_first", &RootedTreeCollection::DropFirst,
           "Drop the first ``fraction`` trees from the tree collection.",
           py::arg("fraction"))
      .def("newick", &RootedTreeCollection::Newick,
           "Get the current set of trees as a big Newick string.")
      .def_readwrite("trees", &RootedTreeCollection::trees_);

  // CLASS
  // RootedTreeGradient
  py::class_<RootedTreeGradient>(m, "RootedTreeGradient",
                                 R"raw(A rooted tree phylogenetic gradient.)raw")
      .def_readonly("log_likelihood", &RootedTreeGradient::log_likelihood_)
      .def_readonly("site_model", &RootedTreeGradient::site_model_)
      .def_readonly("substitution_model", &RootedTreeGradient::substitution_model_)
      .def_readonly("branch_lengths", &RootedTreeGradient::branch_lengths_)
      .def_readonly("clock_model", &RootedTreeGradient::clock_model_)
      .def_readonly("ratios_root_height", &RootedTreeGradient::ratios_root_height_);

  // CLASS
  // UnrootedTree
  py::class_<UnrootedTree>(m, "UnrootedTree", "An unrooted tree with branch lengths.",
                           py::buffer_protocol())
      .def("parent_id_vector", &UnrootedTree::ParentIdVector)
      .def_static("of_parent_id_vector", &UnrootedTree::OfParentIdVector)
      .def_readwrite("branch_lengths", &UnrootedTree::branch_lengths_);

  // CLASS
  // UnrootedTreeCollection
  py::class_<UnrootedTreeCollection>(m, "UnrootedTreeCollection", R"raw(
  A collection of unrooted trees.

  In addition to the methods, UnrootedTreeCollection also offers direct access to
  the trees through the ``trees`` member variable.
  )raw")
      .def(py::init<UnrootedTree::UnrootedTreeVector>(), "The empty constructor.")
      .def(py::init<UnrootedTree::UnrootedTreeVector, TagStringMap>(),
           "Constructor from a vector of trees and a tags->taxon names map.")
      .def(py::init<UnrootedTree::UnrootedTreeVector,
                    const std::vector<std::string> &>(),
           "Constructor from a vector of trees and a vector of taxon names.")
      .def("erase", &UnrootedTreeCollection::Erase,
           "Erase the specified range from the current tree collection.")
      .def("drop_first", &UnrootedTreeCollection::DropFirst,
           "Drop the first ``fraction`` trees from the tree collection.",
           py::arg("fraction"))
      .def("newick", &UnrootedTreeCollection::Newick,
           "Get the current set of trees as a big Newick string.")
      .def_readwrite("trees", &UnrootedTreeCollection::trees_);

  // UnrootedTreeGradient
  py::class_<UnrootedTreeGradient>(m, "UnrootedTreeGradient",
                                   R"raw(An unrooted tree phylogenetic gradient.)raw")
      .def_readonly("log_likelihood", &UnrootedTreeGradient::log_likelihood_)
      .def_readonly("site_model", &UnrootedTreeGradient::site_model_)
      .def_readonly("substitution_model", &UnrootedTreeGradient::substitution_model_)
      .def_readonly("branch_lengths", &UnrootedTreeGradient::branch_lengths_);

  // CLASS
  // PSPIndexer
  py::class_<PSPIndexer>(m, "PSPIndexer", "The primary split pair indexer.")
      .def("details", &PSPIndexer::Details);

  // CLASS
  // PhyloModelSpecification
  py::class_<PhyloModelSpecification>(m, "PhyloModelSpecification",
                                      R"raw(
    Phylogenetic model specification.

    This is how we specify phylogenetic models, with strings for the substitution
    model, the site model, and the clock model.
      )raw")
      .def(py::init<const std::string &, const std::string &, const std::string &>(),
           py::arg("substitution"), py::arg("site"), py::arg("clock"));

  // CLASS
  // SBNInstance
  // So, we'd like to have functionality be shared here between RootedSBNInstance and
  // UnrootedSBNInstance, and there doesn't appear to be a way to have that happen
  // without actually giving it a name. We just give it an underscore-name.
  py::class_<SBNInstance> sbn_instance_class(
      m, "_instance",
      "A general SBNInstance class (that just holds shared functionality between the "
      "two types of SBN instances). Don't use directly, but read these docs to learn "
      "about the relevant functionality in the various subclasses.");
  sbn_instance_class.def(py::init<const std::string &>())
      .def("get_phylo_model_params", &SBNInstance::GetPhyloModelParams)
      .def("get_phylo_model_param_block_map",
           &UnrootedSBNInstance::GetPhyloModelParamBlockMap)
      .def("prepare_for_phylo_likelihood",
           &UnrootedSBNInstance::PrepareForPhyloLikelihood,
           R"raw(
            Prepare instance for phylogenetic likelihood computation.

            See the ``libsbn.beagle_flags`` online documentation to learn about the allowable flags.

            ``use_tip_states`` tells BEAGLE if it should use tip states (versus tip partials).
            Note that libsbn currently treats degenerate nucleotides as gaps irrespective of this setting.

            ``tree_count_option`` tells libsbn for how many trees you will be asking for the likelihood
            or gradient at a time. If not specified, this is set to the number of trees currently loaded
            into the instance. This allocates the correct number of slots in the phylogenetic model
            parameter matrices, and it's up to the user to set those model parameters after calling
            this function.
            Note that this tree count need not be the same as the number of threads (and is typically bigger).
           )raw",
           py::arg("model_specification"), py::arg("thread_count"),
           py::arg("beagle_flags") = std::vector<BeagleFlags>(),
           py::arg("use_tip_states") = true,
           py::arg("tree_count_option") = std::nullopt)
      .def("resize_phylo_model_params", &UnrootedSBNInstance::ResizePhyloModelParams,
           "Resize phylo_model_params.", py::arg("tree_count_option") = std::nullopt)
      .def("read_fasta_file", &SBNInstance::ReadFastaFile,
           "Read a sequence alignment from a FASTA file.")
      // Member Variables
      .def_readonly("psp_indexer", &SBNInstance::psp_indexer_)
      .def_readonly("taxon_names", &SBNInstance::taxon_names_);
  def_read_write_mutable(sbn_instance_class, "sbn_parameters",
                         &SBNInstance::sbn_parameters_);

  // CLASS
  // RootedSBNInstance
  py::class_<RootedSBNInstance, SBNInstance> rooted_sbn_instance_class(
      m, "rooted_instance", R"raw(A rooted SBN instance.)raw");
  rooted_sbn_instance_class
      .def(py::init<const std::string &>())

      // ** Initialization and status
      .def("print_status", &RootedSBNInstance::PrintStatus,
           "Print information about the instance.")
      .def("tree_count", &RootedSBNInstance::TreeCount,
           "Return the number of trees that are currently stored in the "
           "instance.")

      // ** Phylogenetic likelihood
      .def("log_likelihoods", &RootedSBNInstance::LogLikelihoods,
           "Calculate log likelihoods for the current set of trees.")
      .def("set_rescaling", &RootedSBNInstance::SetRescaling,
           "Set whether BEAGLE's likelihood rescaling is used.")
      .def("phylo_gradients", &RootedSBNInstance::PhyloGradients,
           "Calculate gradients of parameters for the current set of trees.")

      // ** I/O
      .def("read_newick_file", &RootedSBNInstance::ReadNewickFile,
           "Read trees from a Newick file.")
      .def("read_nexus_file", &RootedSBNInstance::ReadNexusFile,
           "Read trees from a Nexus file.")

      // ** Member variables
      .def_readwrite("tree_collection", &RootedSBNInstance::tree_collection_);

  // CLASS
  // UnrootedSBNInstance
  py::class_<UnrootedSBNInstance, SBNInstance> unrooted_sbn_instance_class(
      m, "unrooted_instance", R"raw(A unrooted SBN instance.)raw");
  unrooted_sbn_instance_class
      .def(py::init<const std::string &>())

      // ** Initialization and status
      .def("print_status", &UnrootedSBNInstance::PrintStatus,
           "Print information about the instance.")
      .def("tree_count", &UnrootedSBNInstance::TreeCount,
           "Return the number of trees that are currently stored in the "
           "instance.")

      // ** SBN-related items
      .def("process_loaded_trees", &UnrootedSBNInstance::ProcessLoadedTrees, R"raw(
          Process the trees currently stored in the instance.

          Specifically, parse them and build the indexers and the ``sbn_parameters`` vector.
      )raw")
      .def("get_indexers", &UnrootedSBNInstance::GetIndexers,
           "Return the indexer and parent_to_range as string-keyed maps.")
      .def("train_simple_average", &UnrootedSBNInstance::TrainSimpleAverage,
           R"raw(
           Train the SBN using the "simple average" estimator.

           This is described in the "Maximum Lower Bound Estimates" section of the 2018
           NeurIPS paper, and is later referred to as the "SBN-SA" estimator.
           )raw")
      .def("train_expectation_maximization",
           &UnrootedSBNInstance::TrainExpectationMaximization,
           R"raw(
           Train the SBN using the expectation-maximization estimator.

           This is described in the "Expectation Maximization" section of the 2018
           NeurIPS paper, and is later referred to as the "SBN-EM" estimator.

           Here we can supply alpha, the absolute maxiumum number of iterations, and
           a score-based termination criterion for EM. EM will stop if the scaled
           score increase is less than the provided ``score_epsilon``.
           )raw",
           py::arg("alpha"), py::arg("max_iter"), py::arg("score_epsilon") = 0.)
      .def("calculate_sbn_probabilities",
           &UnrootedSBNInstance::CalculateSBNProbabilities,
           R"raw(Get the SBN probabilities of the currently loaded trees.)raw")
      .def("sample_trees", &UnrootedSBNInstance::SampleTrees,
           "Sample trees from the SBN and store them internally.", py::arg("count"))
      .def("make_indexer_representations",
           &UnrootedSBNInstance::MakeIndexerRepresentations,
           R"raw(
            Make the indexer representation of each currently stored tree.

            See the comment for ``IndexerRepresentationOf`` in ``sbn_maps.hpp`` to learn about what that means.

            Note: any rootsplit or a PCSS that is not contained in the subsplit support is given an index equal
            to the length of ``sbn_parameters``. No warning is given.
           )raw")
      .def("make_psp_indexer_representations",
           &UnrootedSBNInstance::MakePSPIndexerRepresentations, R"raw(
            Make the PSP indexer representation of each currently stored tree.

            See the comments in ``psp_indexer.hpp`` to understand the layout.
           )raw")
      .def("split_lengths", &UnrootedSBNInstance::SplitLengths,
           "Get the lengths of the current set of trees, indexed by splits.")
      .def("split_counters", &UnrootedSBNInstance::SplitCounters,
           "A testing method to count splits.")

      // ** Phylogenetic likelihood
      .def("log_likelihoods", &UnrootedSBNInstance::LogLikelihoods,
           "Calculate log likelihoods for the current set of trees.")
      .def("set_rescaling", &UnrootedSBNInstance::SetRescaling,
           "Set whether BEAGLE's likelihood rescaling is used.")
      .def("gradients", &UnrootedSBNInstance::PhyloGradients,
           "Calculate gradients of parameters for the current set of trees.")
      .def("topology_gradients", &UnrootedSBNInstance::TopologyGradients,
           R"raw(Calculate gradients of SBN parameters for the current set of trees.
           Should be called after sampling trees and setting branch lengths.)raw")

      // ** I/O
      .def("read_newick_file", &UnrootedSBNInstance::ReadNewickFile,
           "Read trees from a Newick file.")
      .def("read_nexus_file", &UnrootedSBNInstance::ReadNexusFile,
           "Read trees from a Nexus file.")

      // ** Member variables
      .def_readwrite("tree_collection", &UnrootedSBNInstance::tree_collection_);

  // FUNCTIONS
  m.def("ratio_gradient_of_height_gradient", &RatioGradientOfHeightGradientEigen,
        "Obtain a ratio gradient from a height gradient.");

  // If you want to be sure to get all of the stdout and cerr messages, put your
  // Python code in a context like so:
  // `with libsbn.ostream_redirect(stdout=True, stderr=True):`
  // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/utilities.html#capturing-standard-output-from-ostream
  py::add_ostream_redirect(m, "ostream_redirect");

  // MODULE
  py::module beagle_flags = m.def_submodule("beagle_flags",
                                            R"raw(
      Flags that can be passed to BEAGLE.

      They are used in Python like ``beagle_flags.PROCESSOR_GPU``.

      Note that we expose only a subset of the BEAGLE flags on purpose.
      )raw");
  py::enum_<BeagleFlags>(beagle_flags, "beagle_flag")
      .value("PRECISION_SINGLE", BEAGLE_FLAG_PRECISION_SINGLE,
             "Single precision computation")
      .value("PRECISION_DOUBLE", BEAGLE_FLAG_PRECISION_DOUBLE,
             "Double precision computation")
      .value("COMPUTATION_SYNCH", BEAGLE_FLAG_COMPUTATION_SYNCH,
             "Synchronous computation (blocking)")
      .value("COMPUTATION_ASYNCH", BEAGLE_FLAG_COMPUTATION_ASYNCH,
             "Asynchronous computation (non-blocking)")
      .value("VECTOR_SSE", BEAGLE_FLAG_VECTOR_SSE, "SSE computation")
      .value("VECTOR_NONE", BEAGLE_FLAG_VECTOR_NONE, "No vector computation")
      .value("THREADING_CPP", BEAGLE_FLAG_THREADING_CPP, "C++11 threading")
      .value("THREADING_OPENMP", BEAGLE_FLAG_THREADING_OPENMP, "OpenMP threading")
      .value("THREADING_NONE", BEAGLE_FLAG_THREADING_NONE, "No threading (default)")
      .value("PROCESSOR_CPU", BEAGLE_FLAG_PROCESSOR_CPU, "Use CPU as main processor")
      .value("PROCESSOR_GPU", BEAGLE_FLAG_PROCESSOR_GPU, "Use GPU as main processor")
      .value("FRAMEWORK_CUDA", BEAGLE_FLAG_FRAMEWORK_CUDA,
             "Use CUDA implementation with GPU resources")
      .value("FRAMEWORK_OPENCL", BEAGLE_FLAG_FRAMEWORK_OPENCL,
             "Use OpenCL implementation with GPU resources")
      .value("FRAMEWORK_CPU", BEAGLE_FLAG_FRAMEWORK_CPU, "Use CPU implementation")
      .value("PARALLELOPS_STREAMS", BEAGLE_FLAG_PARALLELOPS_STREAMS,
             "Operations in updatePartials may be assigned to separate device streams")
      .value("PARALLELOPS_GRID", BEAGLE_FLAG_PARALLELOPS_GRID,
             "Operations in updatePartials may be folded into single kernel launch "
             "(necessary for partitions; typically performs better for problems with "
             "fewer pattern sites)")
      .export_values();
}
