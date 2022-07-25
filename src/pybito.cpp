// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#pragma GCC diagnostic pop

#include <string>

#include "gp_instance.hpp"
#include "phylo_flags.hpp"
#include "rooted_gradient_transforms.hpp"
#include "rooted_sbn_instance.hpp"
#include "unrooted_sbn_instance.hpp"

namespace py = pybind11;

// This is how we can have Eigen objects be directly mutable from Python. See
// https://github.com/eacousineau/repro/blob/f4ba595d077af7363f501f6c85d3d2449219f04a/python/pybind11/custom_tests/test_tmp.cc#L16-L38
// Thanks to @eacousineau!
template <typename PyClass, typename C, typename D>
void def_read_write_mutable(PyClass &cls, const char *name, D C::*pm) {
  cls.def_property(
      name, [pm](C & self) -> auto & { return self.*pm; },
      [pm](C &self, const D &value) { self.*pm = value; });
}

// Helper for adding definitions to class.
template <typename PyClass, typename... PyArgTypes, typename CppClass, typename RetType,
          typename... ArgTypes>
void def_template(PyClass pyclass, const char *name, const char *description,
                  RetType (CppClass::*func)(ArgTypes...),
                  std::tuple<PyArgTypes...> pyargs) {
  std::apply(
      [&pyclass, &name, &description, &func](auto &&...pyargs) {
        pyclass.def(
            name,
            [func](CppClass &self, ArgTypes... args) { return (self.*func)(args...); },
            description, pyargs...);
      },
      pyargs);
}

// Define pyclass function for all function overloads (non-const methods).
template <typename PyClass, typename... PyArgTypes, typename CppFunc,
          typename... OtherCppFuncs>
void def_overload(PyClass pyclass, const char *name, const char *description,
                  std::tuple<CppFunc, std::tuple<PyArgTypes...>> overload_def,
                  OtherCppFuncs... other_overloads) {
  // Add definition to class.
  auto &[func, pyargs] = overload_def;
  def_template(pyclass, name, description, func, pyargs);
  // Get next function from template list.
  if constexpr (sizeof...(OtherCppFuncs) > 0) {
    def_overload(pyclass, name, description, other_overloads...);
  }
}

// Use same function name, description, pyargs for multiple functions from multiple
// classes (non-const methods).
template <typename PyClass, typename... PyArgTypes, typename CppFunc,
          typename... OtherClassDefs>
void def_multiclass(const char *name, const char *description,
                    std::tuple<PyArgTypes...> pyargs,
                    std::tuple<PyClass, CppFunc> class_def,
                    OtherClassDefs... other_defs) {
  // Add definition to class.
  auto &[pyclass, func] = class_def;
  def_template(pyclass, name, description, func, pyargs);
  // Get next function from template list.
  if constexpr (sizeof...(OtherClassDefs) > 0) {
    def_multiclass(name, description, pyargs, other_defs...);
  }
}

// In order to make vector<double>s available to numpy, we take two steps.
// First, we make them opaque to pybind11, so that it doesn't do its default
// conversion of STL types.
PYBIND11_MAKE_OPAQUE(std::vector<double>);

// MODULE
PYBIND11_MODULE(bito, m) {
  m.doc() = R"raw(Python interface to bito.)raw";

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

  // PhyloGradient
  py::class_<PhyloGradient>(m, "PhyloGradient", R"raw(A phylogenetic gradient.)raw")
      .def_readonly("log_likelihood", &PhyloGradient::log_likelihood_)
      .def_readonly("gradient", &PhyloGradient::gradient_);

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

            See the ``bito.beagle_flags`` online documentation to learn about the allowable flags.

            ``use_tip_states`` tells BEAGLE if it should use tip states (versus tip partials).
            Note that bito currently treats degenerate nucleotides as gaps irrespective of this setting.

            ``tree_count_option`` tells bito for how many trees you will be asking for the likelihood
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
      .def("load_duplicates_of_first_tree",
           &RootedSBNInstance::LoadDuplicatesOfFirstTree,
           "Replace all of the loaded trees with duplicates of the first tree.",
           py::arg("number_of_times"))
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
      .def("log_det_jacobian_of_height_transform",
           &RootedSBNInstance::LogDetJacobianHeightTransform,
           "Calculate the log det jacobian of the node height transform.")
      .def("set_rescaling", &RootedSBNInstance::SetRescaling,
           "Set whether BEAGLE's likelihood rescaling is used.")

      // ** Phylogenetic gradients
      .def("gradient_log_det_jacobian_of_height_transform",
           &RootedSBNInstance::GradientLogDeterminantJacobian,
           "Obtain the log determinant of the gradient")

      // ** I/O
      .def("read_newick_file", &RootedSBNInstance::ReadNewickFile,
           "Read trees from a Newick file.")
      .def("read_nexus_file", &RootedSBNInstance::ReadNexusFile,
           "Read trees from a Nexus file.")

      // ** Member variables
      .def_readwrite("tree_collection", &RootedSBNInstance::tree_collection_);

  def_overload(
      rooted_sbn_instance_class, "phylo_gradients",
      "Calculate gradients of parameters for the current set of trees.",
      std::tuple(static_cast<std::vector<PhyloGradient> (RootedSBNInstance::*)(
                     std::optional<PhyloFlags>)>(&RootedSBNInstance::PhyloGradients),
                 std::tuple(py::arg("phylo_flags") = std::nullopt)),
      std::tuple(&RootedSBNInstance::PhyloGradients<StringVector>,
                 std::tuple(py::arg("flag_names"), py::arg("use_defaults") = true)),
      std::tuple(
          &RootedSBNInstance::PhyloGradients<StringBoolVector>,
          std::tuple(py::arg("flag_names_and_set"), py::arg("use_defaults") = true)),
      std::tuple(
          &RootedSBNInstance::PhyloGradients<StringDoubleVector>,
          std::tuple(py::arg("flag_names_and_values"), py::arg("use_defaults") = true)),
      std::tuple(&RootedSBNInstance::PhyloGradients<StringBoolDoubleVector>,
                 std::tuple(py::arg("flag_names_and_set_and_values"),
                            py::arg("use_defaults") = true)));
  def_overload(
      rooted_sbn_instance_class, "log_likelihoods",
      "Calculate log likelihoods for the current set of trees.",
      std::tuple(static_cast<std::vector<double> (RootedSBNInstance::*)(
                     std::optional<PhyloFlags>)>(&RootedSBNInstance::LogLikelihoods),
                 std::tuple(py::arg("phylo_flags") = std::nullopt)),
      std::tuple(&RootedSBNInstance::LogLikelihoods<StringVector>,
                 std::tuple(py::arg("flag_names"), py::arg("use_defaults") = true)),
      std::tuple(
          &RootedSBNInstance::LogLikelihoods<StringBoolVector>,
          std::tuple(py::arg("flag_names_and_set"), py::arg("use_defaults") = true)),
      std::tuple(
          &RootedSBNInstance::LogLikelihoods<StringDoubleVector>,
          std::tuple(py::arg("flag_names_and_values"), py::arg("use_defaults") = true)),
      std::tuple(&RootedSBNInstance::LogLikelihoods<StringBoolDoubleVector>,
                 std::tuple(py::arg("flag_names_and_set_and_values"),
                            py::arg("use_defaults") = true)));

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
      .def("load_duplicates_of_first_tree",
           &UnrootedSBNInstance::LoadDuplicatesOfFirstTree,
           "Replace all of the loaded trees with duplicates of the first tree.",
           py::arg("number_of_times"))
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
      .def("set_rescaling", &UnrootedSBNInstance::SetRescaling,
           "Set whether BEAGLE's likelihood rescaling is used.")

      // ** Phylogenetic gradients
      .def("phylo_gradients",
           static_cast<std::vector<PhyloGradient> (UnrootedSBNInstance::*)(
               std::optional<PhyloFlags>)>(&UnrootedSBNInstance::PhyloGradients),
           "Calculate gradients of parameters for the current set of trees.",
           py::arg("phylo_flags") = std::nullopt)
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

  def_overload(
      unrooted_sbn_instance_class, "phylo_gradients",
      "Calculate gradients of parameters for the current set of trees.",
      std::tuple(static_cast<std::vector<PhyloGradient> (UnrootedSBNInstance::*)(
                     std::optional<PhyloFlags>)>(&UnrootedSBNInstance::PhyloGradients),
                 std::tuple(py::arg("phylo_flags") = std::nullopt)),
      std::tuple(&UnrootedSBNInstance::PhyloGradients<StringVector>,
                 std::tuple(py::arg("flag_names"), py::arg("use_defaults") = true)),
      std::tuple(
          &UnrootedSBNInstance::PhyloGradients<StringBoolVector>,
          std::tuple(py::arg("flag_names_and_set"), py::arg("use_defaults") = true)),
      std::tuple(
          &UnrootedSBNInstance::PhyloGradients<StringDoubleVector>,
          std::tuple(py::arg("flag_names_and_values"), py::arg("use_defaults") = true)),
      std::tuple(&UnrootedSBNInstance::PhyloGradients<StringBoolDoubleVector>,
                 std::tuple(py::arg("flag_names_and_set_and_values"),
                            py::arg("use_defaults") = true)));
  def_overload(
      unrooted_sbn_instance_class, "log_likelihoods",
      "Calculate log likelihoods for the current set of trees.",
      std::tuple(static_cast<std::vector<double> (UnrootedSBNInstance::*)(
                     std::optional<PhyloFlags>)>(&UnrootedSBNInstance::LogLikelihoods),
                 std::tuple(py::arg("phylo_flags") = std::nullopt)),
      std::tuple(&UnrootedSBNInstance::LogLikelihoods<StringVector>,
                 std::tuple(py::arg("flag_names"), py::arg("use_defaults") = true)),
      std::tuple(
          &UnrootedSBNInstance::LogLikelihoods<StringBoolVector>,
          std::tuple(py::arg("flag_names_and_set"), py::arg("use_defaults") = true)),
      std::tuple(
          &UnrootedSBNInstance::LogLikelihoods<StringDoubleVector>,
          std::tuple(py::arg("flag_names_and_values"), py::arg("use_defaults") = true)),
      std::tuple(&UnrootedSBNInstance::LogLikelihoods<StringBoolDoubleVector>,
                 std::tuple(py::arg("flag_names_and_set_and_values"),
                            py::arg("use_defaults") = true)));

  // ** PhyloFlags -- for RootedSBNInstance and UnrootedSBNInstance
  def_multiclass(
      "init_phylo_flags", "Create a PhyloFlags object for instance.", std::tuple<>(),
      std::tuple(unrooted_sbn_instance_class, &PreRootedSBNInstance::MakePhyloFlags),
      std::tuple(rooted_sbn_instance_class, &PreRootedSBNInstance::MakePhyloFlags));
  def_multiclass("set_phylo_defaults", "Set whether to use flag defaults.",
                 std::tuple(py::arg("use_defaults") = true),
                 std::tuple(unrooted_sbn_instance_class,
                            &PreRootedSBNInstance::SetPhyloFlagDefaults),
                 std::tuple(rooted_sbn_instance_class,
                            &PreRootedSBNInstance::SetPhyloFlagDefaults));
  def_multiclass(
      "clear_phylo_flags", "Unset all flag settings.", std::tuple<>(),
      std::tuple(unrooted_sbn_instance_class, &PreRootedSBNInstance::ClearPhyloFlags),
      std::tuple(rooted_sbn_instance_class, &PreRootedSBNInstance::ClearPhyloFlags));
  unrooted_sbn_instance_class.def(
      "set_phylo_flag", &PreUnrootedSBNInstance::SetPhyloFlag,
      "Set function flag for given option.", py::arg("flag_name"),
      py::arg("set_to") = true, py::arg("set_value") = 1.0);
  rooted_sbn_instance_class.def("set_phylo_flag", &PreRootedSBNInstance::SetPhyloFlag,
                                "Set function flag for given option.",
                                py::arg("flag_name"), py::arg("set_to") = true,
                                py::arg("set_value") = 1.0);

  // FUNCTIONS
  m.def("ratio_gradient_of_height_gradient",
        &RootedGradientTransforms::RatioGradientOfHeightGradientEigen,
        "Obtain a ratio gradient from a height gradient.");
  m.def("log_det_jacobian_of_height_transform",
        &RootedGradientTransforms::LogDetJacobianHeightTransform,
        "Obtain the log determinant jacobian of a height transform.");
  m.def("gradient_log_det_jacobian_of_height_transform",
        &RootedGradientTransforms::GradientLogDeterminantJacobian,
        "Obtain the log determinant jacobian of the gradient");

  // CLASS
  // GPInstance
  py::class_<GPInstance> gp_instance_class(m, "gp_instance",
                                           R"raw(A generalized pruning instance.)raw");
  gp_instance_class.def(py::init<const std::string &>())
      .def("print_status", &GPInstance::PrintStatus,
           "Print information about the instance.")
      .def("dag_summary_statistics", &GPInstance::DAGSummaryStatistics,
           "Return summary statistics about the DAG.")
      .def("make_dag", &GPInstance::MakeDAG, "Build subsplit DAG.")
      .def("print_dag", &GPInstance::PrintDAG, "Print the subsplit DAG.")

      // ** I/O
      .def("read_newick_file", &GPInstance::ReadNewickFile,
           "Read trees from a Newick file.")
      .def("read_newick_file_gz", &GPInstance::ReadNewickFileGZ,
           "Read trees from a gzip-ed Newick file.")
      .def("read_nexus_file", &GPInstance::ReadNexusFile,
           "Read trees from a Nexus file.")
      .def("read_nexus_file_gz", &GPInstance::ReadNexusFileGZ,
           "Read trees from a gzip-ed Nexus file.")
      .def("read_fasta_file", &GPInstance::ReadFastaFile,
           "Read a sequence alignment from a FASTA file.")
      .def("sbn_parameters_to_csv", &GPInstance::SBNParametersToCSV,
           R"raw(Write "pretty" formatted SBN parameters to a CSV.)raw")
      .def("sbn_prior_to_csv", &GPInstance::SBNPriorToCSV,
           R"raw(Write "pretty" formatted SBN parameters for the prior to a CSV.)raw")
      .def("branch_lengths_to_csv", &GPInstance::BranchLengthsToCSV,
           R"raw(Write "pretty" formatted branch lengths to a CSV.)raw")
      .def("per_gpcsp_llhs_to_csv", &GPInstance::PerGPCSPLogLikelihoodsToCSV,
           R"raw(Write "pretty" formatted per pcsp likelihoods to CSV.)raw")
      .def(
          "intermediate_bls_to_csv", &GPInstance::IntermediateBranchLengthsToCSV,
          R"raw(Write "pretty" formatted per pcsp branch lengths throughout optimization to CSV.)raw")
      .def(
          "intermediate_per_gpcsp_llhs_to_csv",
          &GPInstance::IntermediatePerGPCSPLogLikelihoodsToCSV,
          R"raw(Write "pretty" formatted per pcsp log likelihoods throughout optimization to CSV.)raw")
      .def("per_gpcsp_llh_surfaces_to_csv",
           &GPInstance::PerGPCSPLogLikelihoodSurfacesToCSV,
           R"raw(Write "pretty" formatted per pcsp log likelihood surfaces to CSV.)raw")
      .def(
          "tracked_optim_values_to_csv", &GPInstance::TrackedOptimizationValuesToCSV,
          R"raw(Write "pretty" formatted per pcsp branch lengths and llh values tracked from optimization to CSV.)raw")
      .def("export_trees", &GPInstance::ExportTrees,
           R"raw(Write out currently loaded trees to a Newick file
          (using current GP branch lengths).)raw",
           py::arg("out_path"))
      .def(
          "export_all_generated_topologies", &GPInstance::ExportAllGeneratedTopologies,
          R"raw(Write out all topologies spanned by the current SBN DAG to a Newick file.
            Doesn't require an Engine.)raw",
          py::arg("out_path"))
      .def("export_all_generated_trees", &GPInstance::ExportAllGeneratedTrees,
           R"raw(Write out all trees spanned by the current SBN DAG to a Newick file
          (using current GP branch lengths). Requires an Engine.)raw",
           py::arg("out_path"))
      .def(
          "export_trees_with_a_pcsp", &GPInstance::ExportTreesWithAPCSP,
          R"raw(Write out trees with a given PCSP string to a Newick file (using current
          GP branch lengths).)raw",
          py::arg("pcsp_string"), py::arg("newick_path"))
      .def("subsplit_dag_to_dot", &GPInstance::SubsplitDAGToDot,
           R"raw(Write the current subsplit DAG to a DOT format file.)raw")
      .def("get_branch_lengths", &GPInstance::GetBranchLengths,
           "Return branch lengths from the GPInstance.")
      .def(
          "build_edge_idx_to_pcsp_map",
          [](GPInstance &self) { return self.GetDAG().BuildInverseEdgeIndexer(); },
          "Build a map from SubsplitDAG edge index to its corresponding PCSP bitset.")

      // ** Estimation
      .def("use_gradient_optimization", &GPInstance::UseGradientOptimization,
           "Use gradients for branch length optimization?",
           py::arg("use_gradients") = false)
      .def("make_engine", &GPInstance::MakeEngine, "Prepare for optimization.",
           py::arg("rescaling_threshold") = GPEngine::default_rescaling_threshold_)
      .def("hot_start_branch_lengths", &GPInstance::HotStartBranchLengths,
           "Use given trees to initialize branch lengths.")
      .def("gather_branch_lengths", &GPInstance::GatherBranchLengths,
           "Gather branch lengths into a map keyed by PCSP index for a given tree "
           "sample.")
      .def("calculate_hybrid_marginals", &GPInstance::CalculateHybridMarginals,
           "Calculate hybrid marginals.")
      .def("estimate_sbn_parameters", &GPInstance::EstimateSBNParameters,
           "Estimate the SBN parameters based on current branch lengths.")
      .def("estimate_branch_lengths", &GPInstance::EstimateBranchLengths,
           "Estimate branch lengths for the GPInstance.", py::arg("tol"),
           py::arg("max_iter"), py::arg("quiet") = false,
           py::arg("track_intermediate_iterations") = false)
      .def("get_perpcsp_llh_surface", &GPInstance::GetPerGPCSPLogLikelihoodSurfaces,
           "Scan the likelihood surface for the pcsps in the GPInstance.",
           py::arg("steps"), py::arg("scale_min"), py::arg("scale_max"))
      .def("perturb_and_track_optimization_values",
           &GPInstance::PerturbAndTrackValuesFromOptimization,
           "Reinitiate optimization, perturbing one branch length at a time,  and "
           "track branch length and per pcsp likelihoods.")

      // ** NNI Engine
      .def("make_nni_engine", &GPInstance::MakeNNIEngine,
           R"raw(Initialize NNI Engine.)raw")
      .def(
          "sync_adjacent_nnis_with_dag",
          [](GPInstance &self) { self.GetNNIEngine().SyncAdjacentNNIsWithDAG(); },
          R"raw(Find all adjacent NNIs of DAG.)raw")
      .def(
          "adjacent_nni_count",
          [](GPInstance &self) { self.GetNNIEngine().GetAdjacentNNICount(); },
          R"raw(Get number of adjacent NNIs to current DAG.)raw");

  // If you want to be sure to get all of the stdout and cerr messages, put your
  // Python code in a context like so:
  // `with bito.ostream_redirect(stdout=True, stderr=True):`
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

  // ** Export Keys
  auto ExportFlagsToModuleAttributes = [](py::module &module,
                                          const PhyloFlagOptionSet &flag_set) {
    for (const auto &[name, flag] : flag_set.GetAllNames()) {
      module.attr(name.c_str()) = py::cast(std::string(flag));
    }
  };
  auto ExportMapkeysToModuleAttributes = [](py::module &module,
                                            const PhyloMapkeySet &mapkey_set) {
    for (const auto &[name, mapkey] : mapkey_set.GetAllNames()) {
      std::ignore = name;
      module.attr(mapkey.GetName().c_str()) = py::cast(std::string(mapkey.GetKey()));
    }
  };

  // * Export PhyloFlagOptions
  py::module phylo_flags = m.def_submodule("phylo_flags",
                                           R"raw(
        Option flags for functions such as ``SBNInstance::phylo_gradient`` and ``SBNInstanct::log_likelihood``.
      )raw");
  ExportFlagsToModuleAttributes(phylo_flags, PhyloGradientFlagOptions::set_);
  ExportFlagsToModuleAttributes(phylo_flags, LogLikelihoodFlagOptions::set_);

  // * Export PhyloMapkeys
  py::module phylo_model_mapkeys = m.def_submodule("phylo_model_mapkeys",
                                                   R"raw(
        Dict keys for accessing the PhyloModel, returned by ``SBNInstance::get_phylo_model_param_block_map``.
      )raw");
  ExportMapkeysToModuleAttributes(phylo_model_mapkeys, PhyloModelMapkeys::set_);
  py::module phylo_gradient_mapkeys = m.def_submodule("phylo_gradient_mapkeys",
                                                      R"raw(
        Dict keys for accessing the GradientMap, returned by ``SBNInstance::phylo_gradient``.
      )raw");
  ExportMapkeysToModuleAttributes(phylo_gradient_mapkeys, PhyloGradientMapkeys::set_);
}
