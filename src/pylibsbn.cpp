// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

#include "gp_instance.hpp"
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
      .def("initialize_time_tree_using_height_ratios",
           &RootedTree::InitializeTimeTreeUsingHeightRatios)
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
  // RootedPhyloGradient
  py::class_<RootedPhyloGradient>(m, "RootedPhyloGradient",
                                  R"raw(A rooted tree phylogenetic gradient.)raw")
      .def_readonly("log_likelihood", &RootedPhyloGradient::log_likelihood_)
      .def_readonly("site_model", &RootedPhyloGradient::site_model_)
      .def_readonly("substitution_model", &RootedPhyloGradient::substitution_model_)
      .def_readonly("branch_lengths", &RootedPhyloGradient::branch_lengths_)
      .def_readonly("clock_model", &RootedPhyloGradient::clock_model_)
      .def_readonly("ratios_root_height", &RootedPhyloGradient::ratios_root_height_);

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

  // UnrootedPhyloGradient
  py::class_<UnrootedPhyloGradient>(m, "UnrootedPhyloGradient",
                                    R"raw(An unrooted tree phylogenetic gradient.)raw")
      .def_readonly("log_likelihood", &UnrootedPhyloGradient::log_likelihood_)
      .def_readonly("site_model", &UnrootedPhyloGradient::site_model_)
      .def_readonly("substitution_model", &UnrootedPhyloGradient::substitution_model_)
      .def_readonly("branch_lengths", &UnrootedPhyloGradient::branch_lengths_);

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

  // ** SBNInstance variants

  const char prepare_for_phylo_likelihood_docstring[] =
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
           )raw";

  const char process_loaded_trees_docstring[] = R"raw(
          Process the trees currently stored in the instance.

          Specifically, parse them and build the subsplit support and allocate (but do not train) the corresponding SBN
          parameters.
      )raw";

  const char read_sbn_parameters_from_csv_docstring[] = R"raw(
        Read SBN parameters from a CSV mapping a string representation of the GPCSP to its probability in linear (not
        log) space.

        Any GPCSPs that are not found in the supplied CSV will be assigned a probability of 0.
      )raw";

  // CLASS
  // RootedSBNInstance
  py::class_<PreRootedSBNInstance>(m, "PreRootedSBNInstance");
  py::class_<RootedSBNInstance, PreRootedSBNInstance> rooted_sbn_instance_class(
      m, "rooted_instance", R"raw(
      A rooted SBN instance.

      The intent of this class is primarily to support rooted time trees, however some functionality works without
      dates.

      If you are using time trees, you will need to assign tip dates using one of the available methods and then
      initialize the time trees, either using branch lengths or per-tree height ratios.

      If you don't do this, then trying to calculate a phylogenetic gradient will raise an exception.
      )raw");
  rooted_sbn_instance_class
      .def(py::init<const std::string &>())

      // ** BEGIN DUPLICATED CODE BLOCK between this and UnrootedSBNInstance
      .def("get_phylo_model_params", &RootedSBNInstance::GetPhyloModelParams)
      .def("get_phylo_model_param_block_map",
           &RootedSBNInstance::GetPhyloModelParamBlockMap)
      .def(
          "prepare_for_phylo_likelihood", &RootedSBNInstance::PrepareForPhyloLikelihood,
          prepare_for_phylo_likelihood_docstring, py::arg("model_specification"),
          py::arg("thread_count"), py::arg("beagle_flags") = std::vector<BeagleFlags>(),
          py::arg("use_tip_states") = true, py::arg("tree_count_option") = std::nullopt)
      .def("resize_phylo_model_params", &RootedSBNInstance::ResizePhyloModelParams,
           "Resize phylo_model_params.", py::arg("tree_count_option") = std::nullopt)
      .def("read_fasta_file", &RootedSBNInstance::ReadFastaFile,
           "Read a sequence alignment from a FASTA file.")
      .def("taxon_names", &RootedSBNInstance::TaxonNames,
           "Return a list of taxon names.")

      // ** Initialization and status
      .def("print_status", &RootedSBNInstance::PrintStatus,
           "Print information about the instance.")
      .def("tree_count", &RootedSBNInstance::TreeCount,
           "Return the number of trees that are currently stored in the "
           "instance.")

      // ** SBN-related items
      .def("process_loaded_trees", &RootedSBNInstance::ProcessLoadedTrees,
           process_loaded_trees_docstring)
      .def("train_simple_average", &RootedSBNInstance::TrainSimpleAverage,
           R"raw(
           Train the SBN using the "simple average" estimator.

           For rooted trees this training is simpler than in the unrooted case:
           we simply take the normalized frequency of PCSPs.
           )raw")
      .def("sbn_parameters_to_csv", &RootedSBNInstance::SBNParametersToCSV,
           R"raw(Write "pretty" formatted SBN parameters to a CSV.)raw")
      .def("read_sbn_parameters_from_csv", &RootedSBNInstance::ReadSBNParametersFromCSV,
           read_sbn_parameters_from_csv_docstring)
      .def("calculate_sbn_probabilities", &RootedSBNInstance::CalculateSBNProbabilities,
           R"raw(Calculate the SBN probabilities of the currently loaded trees.)raw")
      // ** END DUPLICATED CODE BLOCK between this and UnrootedSBNInstance
      .def("unconditional_subsplit_probabilities_to_csv",
           &RootedSBNInstance::UnconditionalSubsplitProbabilitiesToCSV,
           "Write out the overall probability of seeing each subsplit when we sample a "
           "tree from the SBN.")

      // ** Tip dates
      .def("set_dates_to_be_constant", &RootedSBNInstance::SetDatesToBeConstant,
           "Set tip dates to be constant.",
           py::arg("initialize_time_trees_using_branch_lengths"))
      .def("parse_dates_from_taxon_names", &RootedSBNInstance::ParseDatesFromTaxonNames,
           "Take dates to be the numbers after an underscore in the taxon names.",
           py::arg("initialize_time_trees_using_branch_lengths"))
      .def(
          "parse_dates_from_csv", &RootedSBNInstance::ParseDatesFromCSV,
          "Parse dates from a headerless 2-column CSV of quoted taxon names and dates.",
          py::arg("csv_path"), py::arg("initialize_time_trees_using_branch_lengths"))

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
  py::class_<PreUnrootedSBNInstance>(m, "PreUnrootedSBNInstance");
  py::class_<UnrootedSBNInstance, PreUnrootedSBNInstance> unrooted_sbn_instance_class(
      m, "unrooted_instance", R"raw(A unrooted SBN instance.)raw");
  unrooted_sbn_instance_class
      .def(py::init<const std::string &>())

      // ** BEGIN DUPLICATED CODE BLOCK between this and RootedSBNInstance
      .def("get_phylo_model_params", &UnrootedSBNInstance::GetPhyloModelParams)
      .def("get_phylo_model_param_block_map",
           &UnrootedSBNInstance::GetPhyloModelParamBlockMap)
      .def(
          "prepare_for_phylo_likelihood",
          &UnrootedSBNInstance::PrepareForPhyloLikelihood,
          prepare_for_phylo_likelihood_docstring, py::arg("model_specification"),
          py::arg("thread_count"), py::arg("beagle_flags") = std::vector<BeagleFlags>(),
          py::arg("use_tip_states") = true, py::arg("tree_count_option") = std::nullopt)
      .def("resize_phylo_model_params", &UnrootedSBNInstance::ResizePhyloModelParams,
           "Resize phylo_model_params.", py::arg("tree_count_option") = std::nullopt)
      .def("read_fasta_file", &UnrootedSBNInstance::ReadFastaFile,
           "Read a sequence alignment from a FASTA file.")
      .def("taxon_names", &UnrootedSBNInstance::TaxonNames,
           "Return a list of taxon names.")

      // ** Initialization and status
      .def("print_status", &UnrootedSBNInstance::PrintStatus,
           "Print information about the instance.")
      .def("tree_count", &UnrootedSBNInstance::TreeCount,
           "Return the number of trees that are currently stored in the "
           "instance.")

      // ** SBN-related items
      .def("process_loaded_trees", &UnrootedSBNInstance::ProcessLoadedTrees,
           process_loaded_trees_docstring)
      .def("train_simple_average", &UnrootedSBNInstance::TrainSimpleAverage,
           R"raw(
           Train the SBN using the "simple average" estimator.

           This is described in the "Maximum Lower Bound Estimates" section of the 2018
           NeurIPS paper, and is later referred to as the "SBN-SA" estimator.
           )raw")
      .def("sbn_parameters_to_csv", &UnrootedSBNInstance::SBNParametersToCSV,
           R"raw(Write "pretty" formatted SBN parameters to a CSV.)raw")
      .def("read_sbn_parameters_from_csv",
           &UnrootedSBNInstance::ReadSBNParametersFromCSV,
           read_sbn_parameters_from_csv_docstring)
      .def("calculate_sbn_probabilities",
           &UnrootedSBNInstance::CalculateSBNProbabilities,
           R"raw(Calculate the SBN probabilities of the currently loaded trees.)raw")
      // ** END DUPLICATED CODE BLOCK between this and RootedSBNInstance

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
      .def("sample_trees", &UnrootedSBNInstance::SampleTrees,
           "Sample trees from the SBN and store them internally.", py::arg("count"))
      .def("make_indexer_representations",
           &UnrootedSBNInstance::MakeIndexerRepresentations,
           R"raw(
            Make the indexer representation of each currently stored tree.

            See the comment for ``IndexerRepresentationOf`` in ``sbn_maps.hpp`` to learn about what that means.

            Note: any rootsplit or a PCSP that is not contained in the subsplit support is given an index equal
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
      .def("phylo_gradients", &UnrootedSBNInstance::PhyloGradients,
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
      .def_readonly("psp_indexer", &UnrootedSBNInstance::psp_indexer_)
      .def_readwrite("tree_collection", &UnrootedSBNInstance::tree_collection_);

  def_read_write_mutable(unrooted_sbn_instance_class, "sbn_parameters",
                         &UnrootedSBNInstance::sbn_parameters_);

  // FUNCTIONS
  m.def("ratio_gradient_of_height_gradient", &RatioGradientOfHeightGradientEigen,
        "Obtain a ratio gradient from a height gradient.");

  // CLASS
  // GPInstance
  py::class_<GPInstance> gp_instance_class(m, "gp_instance",
                                           R"raw(A generalized pruning instance.)raw");
  gp_instance_class.def(py::init<const std::string &>())
      .def("print_status", &GPInstance::PrintStatus,
           "Print information about the instance.")
      .def("print_dag", &GPInstance::PrintDAG, "Print the generalized pruning DAG.")

      // ** I/O
      .def("read_newick_file", &GPInstance::ReadNewickFile,
           "Read trees from a Newick file.")
      .def("read_nexus_file", &GPInstance::ReadNexusFile,
           "Read trees from a Nexus file.")
      .def("read_fasta_file", &GPInstance::ReadFastaFile,
           "Read a sequence alignment from a FASTA file.")
      .def("sbn_parameters_to_csv", &GPInstance::SBNParametersToCSV,
           R"raw(Write "pretty" formatted SBN parameters to a CSV.)raw")
      .def("sbn_prior_to_csv", &GPInstance::SBNPriorToCSV,
           R"raw(Write "pretty" formatted SBN parameters for the prior to a CSV.)raw")
      .def("branch_lengths_to_csv", &GPInstance::BranchLengthsToCSV,
           R"raw(Write "pretty" formatted branch lengths to a CSV.)raw")
      .def(
          "export_trees_with_a_pcsp", &GPInstance::ExportTreesWithAPCSP,
          R"raw(Write out trees with a given PCSP string to a Newick file (using current
          GP branch lengths).)raw",
          py::arg("pcsp_string"), py::arg("newick_path"))

      // ** Estimation
      .def("make_engine", &GPInstance::MakeEngine, "Prepare for optimization.",
           py::arg("rescaling_threshold") = GPEngine::default_rescaling_threshold_)
      .def("hot_start_branch_lengths", &GPInstance::HotStartBranchLengths,
           "Use given trees to initialize branch lengths.")
      .def("estimate_sbn_parameters", &GPInstance::EstimateSBNParameters,
           "EstimateSBNParameters.")
      .def("estimate_branch_lengths", &GPInstance::EstimateBranchLengths,
           "Estimate branch lengths for the GPInstance.", py::arg("tol"),
           py::arg("max_iter"), py::arg("quiet") = false);

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
