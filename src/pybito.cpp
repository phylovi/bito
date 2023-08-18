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
      .def("__eq__", [](const RootedTree &self,
                        const RootedTree &other) { return self == other; })
      .def("compare_by_topology",
           [](const RootedTree &self, const RootedTree &other) {
             return self.Topology() == other.Topology();
           })
      .def(
          "to_newick", [](const RootedTree &self) { return self.Newick(); },
          "Output to Newick string with branch lengths.")
      .def(
          "to_newick_topology",
          [](const RootedTree &self) { return self.NewickTopology(std::nullopt); },
          "Output to Newick string without branch lengths.")
      .def("parent_id_vector", &RootedTree::ParentIdVector)
      .def("initialize_time_tree_using_height_ratios",
           &RootedTree::InitializeTimeTreeUsingHeightRatios)
      .def_static("example", &RootedTree::Example)
      .def_static("of_parent_id_vector", &RootedTree::OfParentIdVector)
      .def_readwrite("branch_lengths", &RootedTree::branch_lengths_)
      .def_readwrite("height_ratios", &RootedTree::height_ratios_)
      .def_readwrite("node_heights", &RootedTree::node_heights_)
      .def_readwrite("node_bounds", &RootedTree::node_bounds_)
      .def_readwrite("rates", &RootedTree::rates_)
      .def(
          "topology", [](const RootedTree &self) { return self.Topology(); },
          py::return_value_policy::reference)
      .def(
          "id", [](const RootedTree &self) { return self.Topology()->Id(); },
          "Unique node id within topology.")
      .def(
          "to_leaves", [](const RootedTree &self) { return self.Topology()->Leaves(); },
          "Output node to leave bitset.")
      .def(
          "build_subsplit",
          [](const RootedTree &self) { return self.Topology()->BuildSubsplit(); },
          "Build subsplit node bitset of node.")
      .def(
          "build_pcsp",
          [](const RootedTree &self, const size_t child_id) {
            Assert(child_id < 2, "child_count must be 0 (left) or 1 (right).");
            auto child_clade =
                (child_id == 0) ? SubsplitClade::Left : SubsplitClade::Right;
            return self.Topology()->BuildPCSP(child_clade);
          },
          "Build PCSP edge bitset of edge below node.")
      .def(
          "build_set_of_subsplits",
          [](const RootedTree &self) { return self.Topology()->BuildSetOfSubsplits(); },
          "Build set of all subsplit bitsets for all nodes in topology.")
      .def(
          "build_set_of_pcsps",
          [](const RootedTree &self) { return self.Topology()->BuildSetOfPCSPs(); },
          "Build set of all PCSP edge bitsets for all edges in topology.");

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
      .def("__eq__", [](const UnrootedTree &self,
                        const UnrootedTree &other) { return self == other; })
      .def(
          "to_newick", [](const UnrootedTree &self) { return self.Newick(); },
          "Output to Newick string with branch lengths.")
      .def(
          "to_newick_topology",
          [](const UnrootedTree &self) { return self.NewickTopology(std::nullopt); },
          "Output to Newick string without branch lengths.")
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
      .def("read_newick_file", &RootedSBNInstance::ReadNewickFile, py::arg("path"),
           py::arg("sort_taxa") = true, "Read trees from a Newick file.")
      .def("read_nexus_file", &RootedSBNInstance::ReadNexusFile, py::arg("path"),
           py::arg("sort_taxa") = true, "Read trees from a Nexus file.")

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
      m, "unrooted_instance", R"raw(An unrooted SBN instance.)raw");
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
      .def("read_newick_file", &UnrootedSBNInstance::ReadNewickFile, py::arg("path"),
           py::arg("sort_taxa") = true, "Read trees from a Newick file.")
      .def("read_nexus_file", &UnrootedSBNInstance::ReadNexusFile, py::arg("path"),
           py::arg("sort_taxa") = true, "Read trees from a Nexus file.")

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
      .def("read_newick_file", &GPInstance::ReadNewickFile, py::arg("path"),
           py::arg("sort_taxa") = true, "Read trees from a Newick file.")
      .def("read_newick_file_gz", &GPInstance::ReadNewickFileGZ, py::arg("path"),
           py::arg("sort_taxa") = true, "Read trees from a gzip-ed Newick file.")
      .def("read_nexus_file", &GPInstance::ReadNexusFile, py::arg("path"),
           py::arg("sort_taxa") = true, "Read trees from a Nexus file.")
      .def("read_nexus_file_gz", &GPInstance::ReadNexusFileGZ, py::arg("path"),
           py::arg("sort_taxa") = true, "Read trees from a gzip-ed Nexus file.")
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
      .def("currently_loaded_trees_with_gp_branch_lengths",
           &GPInstance::CurrentlyLoadedTreesWithGPBranchLengths,
           "Collection of all rooted trees loaded into DAG.")
      .def("generate_complete_rooted_tree_collection",
           &GPInstance::GenerateCompleteRootedTreeCollection,
           "Generate collection of all rooted trees expressed in DAG.")
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
          "Build a map from DAG edge index to its corresponding PCSP bitset.")

      // ** Estimation
      .def("use_gradient_optimization", &GPInstance::UseGradientOptimization,
           "Use gradients for branch length optimization?",
           py::arg("use_gradients") = false)
      .def("hot_start_branch_lengths", &GPInstance::HotStartBranchLengths,
           "Use given trees to initialize branch lengths.")
      .def("gather_branch_lengths", &GPInstance::GatherBranchLengths,
           "Gather branch lengths into a map keyed by PCSP index for a given tree "
           "sample.")
      .def("calculate_hybrid_marginals", &GPInstance::CalculateHybridMarginals,
           "Calculate hybrid marginals.")
      .def("estimate_sbn_parameters", &GPInstance::EstimateSBNParameters,
           "Estimate the SBN parameters based on current branch lengths.")
      .def("hot_start_branch_length", &GPInstance::HotStartBranchLengths)
      .def("take_first_branch_length", &GPInstance::TakeFirstBranchLength)
      .def("estimate_branch_lengths", &GPInstance::EstimateBranchLengths,
           "Estimate branch lengths for the GPInstance.", py::arg("tol"),
           py::arg("max_iter"), py::arg("quiet") = false,
           py::arg("track_intermediate_iterations") = false,
           py::arg("optimization_method") = std::nullopt)
      .def("get_perpcsp_llh_surface", &GPInstance::GetPerGPCSPLogLikelihoodSurfaces,
           "Scan the likelihood surface for the pcsps in the GPInstance.",
           py::arg("steps"), py::arg("scale_min"), py::arg("scale_max"))
      .def("perturb_and_track_optimization_values",
           &GPInstance::PerturbAndTrackValuesFromOptimization,
           "Reinitiate optimization, perturbing one branch length at a time,  and "
           "track branch length and per pcsp likelihoods.")

      // ** GP Likelihoods
      .def("populate_plvs", &GPInstance::PopulatePLVs, "Populate PLVs.")
      .def("compute_likelihoods", &GPInstance::ComputeLikelihoods,
           "Compute Likelihoods.")
      .def("get_per_pcsp_log_likelihoods", &GPInstance::GetPerPCSPLogLikelihoods,
           "Get Per-PCSP Log Likelihoods.")

      // ** DAG
      .def("make_dag", &GPInstance::MakeDAG, "Initialize Subsplit DAG.")
      .def(
          "get_dag", [](GPInstance &self) -> GPDAG * { return &self.GetDAG(); },
          py::return_value_policy::reference, "Get Subsplit DAG.")

      // ** DAG Engines
      .def("make_gp_engine", &GPInstance::MakeGPEngine, "Initialize GP Engine.",
           py::arg("rescaling_threshold") = GPEngine::default_rescaling_threshold_,
           py::arg("use_gradients") = false)
      // .def("get_gp_engine", &GPInstance::GetGPEngine,
      //      py::return_value_policy::reference, "Get GP Engine.")
      .def(
          "get_gp_engine",
          [](GPInstance &self) -> GPEngine * { return &self.GetGPEngine(); },
          py::return_value_policy::reference, "Get GP Engine.")
      .def("make_nni_engine", &GPInstance::MakeNNIEngine, "Initialize NNI Engine.")
      .def(
          "get_nni_engine",
          [](GPInstance &self) -> NNIEngine * { return &self.GetNNIEngine(); },
          py::return_value_policy::reference, "Get Subsplit DAG.")
      .def("make_tp_engine", &GPInstance::MakeTPEngine, "Initialize TP Engine.")
      .def(
          "get_tp_engine",
          [](GPInstance &self) -> TPEngine * { return &self.GetTPEngine(); },
          py::return_value_policy::reference, "Get TP Engine.")
      .def("tp_engine_set_branch_lengths_by_taking_first",
           &GPInstance::TPEngineSetBranchLengthsByTakingFirst)
      .def("tp_engine_set_choice_map_by_taking_first",
           &GPInstance::TPEngineSetChoiceMapByTakingFirst,
           py::arg("use_subsplit_method") = true)

      // ** Tree Engines
      .def("get_likelihood_tree_engine", &GPInstance::GetLikelihoodTreeEngine,
           py::return_value_policy::reference)
      .def("get_parsimony_tree_engine", &GPInstance::GetParsimonyTreeEngine,
           py::return_value_policy::reference)
      .def("compute_tree_likelihood",
           [](const GPInstance &self, const RootedTree &tree) {
             auto beagle_pref_flags = BEAGLE_FLAG_VECTOR_SSE;
             PhyloModelSpecification model_spec{"JC69", "constant", "strict"};
             SitePattern site_pattern = self.MakeSitePattern();
             bool use_tip_states = true;
             auto tree_engine =
                 FatBeagle(model_spec, site_pattern, beagle_pref_flags, use_tip_states);
             return tree_engine.UnrootedLogLikelihood(tree);
           })
      .def("compute_tree_parsimony",
           [](const GPInstance &self, const RootedTree &tree) {
             auto site_pattern = self.MakeSitePattern();
             auto mmap_file_path = self.GetMMapFilePath() + ".sankoff";
             auto tree_engine = SankoffHandler(site_pattern, mmap_file_path);
             tree_engine.RunSankoff(tree.Topology());
             return tree_engine.ParsimonyScore();
           });

  // ** DAGs

  py::class_<GPDAG> dag_class(m, "dag", "Subsplit DAG for performing GPOperations.");
  dag_class.def("__eq__", [](const GPDAG &self, const GPDAG &other) { self == other; })
      .def("node_count", &GPDAG::NodeCount, "Get number of nodes contained in DAG.")
      .def("edge_count", &GPDAG::EdgeCountWithLeafSubsplits,
           "Get number of edges contained in DAG.")
      .def("taxon_count", &GPDAG::TaxonCount, "Get number of taxa in DAG.")
      .def("topology_count", &GPDAG::TopologyCount,
           "Get number of unique topologies contained in DAG.")
      .def("get_nni", &GPDAG::GetNNI, "Get NNI for the given DAG edge.")
      .def("get_node_id",
           [](const GPDAG &self, const Bitset &bitset) {
             return self.GetDAGNodeId(bitset);
           })
      .def("get_edge_id", [](const GPDAG &self,
                             const Bitset &bitset) { return self.GetEdgeIdx(bitset); })
      .def("get_edge_id", [](const GPDAG &self,
                             const NNIOperation &nni) { return self.GetEdgeIdx(nni); })
      .def("get_taxon_map", &GPDAG::GetTaxonMap,
           "Get map of taxon names contained in DAG.")
      .def("build_set_of_node_bitsets", &GPDAG::BuildSetOfNodeBitsets,
           "Build a set of node Subsplit bitsets contained in DAG.")
      .def("build_set_of_edge_bitsets", &GPDAG::BuildSetOfEdgeBitsets,
           "Build a set of edge PCSP bitsets contained in DAG.")
      .def("contains_node",
           [](const GPDAG &self, const Bitset &bitset) {
             return self.ContainsNode(bitset);
           })
      .def("contains_edge",
           [](const GPDAG &self, const Bitset &bitset) {
             return self.ContainsEdge(bitset);
           })
      .def("contains_nni", &GPDAG::ContainsNNI)
      .def("contains_tree", &GPDAG::ContainsTree, "Check whether DAG contains tree.",
           py::arg("tree"), py::arg("is_quiet") = true)
      .def("contains_topology", &GPDAG::ContainsTopology,
           "Check whether DAG contains topology.")
      .def("is_valid_add_node_pair", &GPDAG::IsValidAddNodePair,
           "Checks whether a given parent/child subsplit pair is valid to be added to "
           "the DAG.")
      .def(
          "add_node_pair",
          [](GPDAG &self, const Bitset &parent, const Bitset &child) {
            self.AddNodePair(parent, child);
          },
          "Add parent/child subsplit pair to DAG.")
      .def("add_nodes", &GPDAG::AddNodes)
      .def("add_edges", &GPDAG::AddEdges)
      .def("fully_connect", &GPDAG::FullyConnect,
           "Adds all valid edges with present nodes to the DAG.")
      // ** I/O
      .def("tree_to_newick_topology", &GPDAG::TreeToNewickTopology)
      .def("tree_to_newick_tree", &GPDAG::TreeToNewickTree)
      .def("topology_to_newick_topology", &GPDAG::TopologyToNewickTopology)
      .def("generate_all_topologies", &GPDAG::GenerateAllTopologies)
      .def("to_newick_of_all_topologies", &GPDAG::ToNewickOfAllTopologies)
      .def("generate_covering_topologies", &GPDAG::GenerateCoveringTopologies)
      .def("to_newick_of_covering_topologies", &GPDAG::ToNewickOfCoveringTopologies);

  py::class_<GraftDAG> graft_dag_class(m, "graft_dag",
                                       "Subsplit DAG for grafting nodes and edges.");
  graft_dag_class
      .def("compare_to_dag",
           [](const GraftDAG &self, const GPDAG &other) { self.CompareToDAG(other); })
      .def("graft_node_count", &GraftDAG::GraftNodeCount,
           "Get number of graft nodes appended to DAG.")
      .def("graft_edge_count", &GraftDAG::GraftEdgeCount,
           "Get number of graft edges appended to DAG.")
      .def("host_node_count", &GraftDAG::HostNodeCount,
           "Get number of host nodes contained in DAG.")
      .def("host_edge_count", &GraftDAG::HostEdgeCount,
           "Get number of host edges contained in DAG.")
      .def("get_host_dag", &GraftDAG::GetHostDAG, py::return_value_policy::reference,
           "Get underlying Host DAG.")
      .def("is_valid_add_node_pair", &GraftDAG::IsValidAddNodePair,
           "Checks whether a given parent/child subsplit pair is valid to be added to "
           "the DAG.")
      .def(
          "add_node_pair",
          [](GraftDAG &self, const Bitset &parent, const Bitset &child) {
            self.AddNodePair(parent, child);
          },
          "Add parent/child subsplit pair to DAG.");

  // ** Engines for DAGs

  py::class_<GPEngine> gp_engine_class(m, "gp_engine",
                                       "An engine for computing Generalized Pruning.");
  gp_engine_class.def("node_count", &GPEngine::GetNodeCount, "Get number of nodes.")
      .def("plv_count", &GPEngine::GetPLVCount, "Get number of PLVs.")
      .def("edge_count", &GPEngine::GetGPCSPCount, "Get number of edges.");

  py::class_<TPEngine> tp_engine_class(m, "tp_engine",
                                       "An engine for computing Top Pruning.");
  tp_engine_class.def("node_count", &TPEngine::GetNodeCount, "Get number of nodes.")
      .def("edge_count", &TPEngine::GetEdgeCount, "Get number of edges.")
      .def("get_top_tree_with_edge", &TPEngine::GetTopTreeWithEdge,
           "Output the top tree of tree containing given edge.")
      .def("get_top_tree_likelihood_with_edge", &TPEngine::GetTopTreeLikelihood,
           "Output the top tree likelihood containing given edge.")
      .def("get_top_tree_parsimony_with_edge", &TPEngine::GetTopTreeParsimony,
           "Output the top tree parsimony containing given edge.")
      .def("get_top_tree_topology_with_edge", &TPEngine::GetTopTopologyWithEdge,
           "Output the top tree of tree containing given edge.")
      // ** Branch Length Optimization
      .def("get_branch_lengths", [](TPEngine &self) { self.GetBranchLengths(); })
      .def("optimize_branch_lengths", &TPEngine::OptimizeBranchLengths,
           py::arg("check_branch_convergence") = std::nullopt)
      // ** I/O
      .def("build_map_of_tree_id_to_top_topologies",
           &TPEngine::BuildMapOfTreeIdToTopTopologies)
      .def("to_newick_of_top_topologies", &TPEngine::ToNewickOfTopTopologies)
      .def("to_newick_of_top_trees", &TPEngine::ToNewickOfTopTrees);

  py::class_<TPChoiceMap> tp_choice_map_class(
      m, "tp_choice_map", "An choice map for finding the top tree in TPEngine.");
  tp_choice_map_class.def("__str__", &TPChoiceMap::ToString)
      .def(
          "edge_choice_to_string",
          [](const TPChoiceMap &self, const EdgeId edge_id) {
            return self.EdgeChoiceToString(edge_id);
          },
          "Output the edge choice map for a given edge in DAG.");

  py::class_<NNIEngine> nni_engine_class(
      m, "nni_engine", "An engine for computing NNI Systematic Search.");
  nni_engine_class
      // Getters
      .def(
          "get_graft_dag", [](NNIEngine &self) { return self.GetGraftDAG(); },
          py::return_value_policy::reference, "Get the Graft DAG.")
      .def(
          "adjacent_nnis", [](NNIEngine &self) { return self.GetAdjacentNNIs(); },
          "Get NNIs adjacent to DAG.")
      .def(
          "accepted_nnis", [](NNIEngine &self) { return self.GetAcceptedNNIs(); },
          "Get NNIs accepted into DAG.")
      .def(
          "rejected_nnis", [](NNIEngine &self) { return self.GetRejectedNNIs(); },
          "Get NNIs rejected from DAG.")
      .def(
          "scored_nnis", [](const NNIEngine &self) { return self.GetScoredNNIs(); },
          "Get Scored NNIs of current iteration.")
      // Counts
      .def("adjacent_nni_count", &NNIEngine::GetAdjacentNNICount,
           "Get number of NNIs adjacent to DAG.")
      .def("accepted_nni_count", &NNIEngine::GetAcceptedNNICount,
           "Get number of adjacent NNIs were accepted by the filter on current "
           "iteration.")
      .def("rejected_nni_count", &NNIEngine::GetRejectedNNICount,
           "Get number of adjacent NNIs were rejected by the filter on current "
           "iteration.")
      .def("past_accepted_nni_count", &NNIEngine::GetPastAcceptedNNICount,
           "Get number of adjacent NNIs were accepted by the filter on all previous "
           "iterations.")
      .def("past_rejected_nni_count", &NNIEngine::GetPastRejectedNNICount,
           "Get number of adjacent NNIs were rejected by the filter on all previous "
           "iterations.")
      .def("past_scored_nnis", &NNIEngine::GetPastScoredNNIs,
           "Get scores from NNIs from previous iterations.")
      .def("iter_count", &NNIEngine::GetIterationCount,
           "Get number of iterations of NNI search run.")
      // Search primary routines
      .def("run", &NNIEngine::Run, "Primary runner for NNI systematic search.",
           py::arg("is_quiet") = false)
      .def("run_init", &NNIEngine::RunInit, "Run initialization step of NNI search.",
           py::arg("is_quiet") = false)
      .def("run_main_loop", &NNIEngine::RunMainLoop, "Run main loop of NNI search.",
           py::arg("is_quiet") = false)
      .def("run_post_loop", &NNIEngine::RunPostLoop, "Run post loop of NNI search.",
           py::arg("is_quiet") = false)
      // Search subroutines
      // Init
      .def("reset_all_nnis", &NNIEngine::ResetAllNNIs)
      .def("sync_adjacent_nnis_with_dag", &NNIEngine::SyncAdjacentNNIsWithDAG)
      .def("prep_eval_engine", &NNIEngine::PrepEvalEngine)
      .def("filter_init", &NNIEngine::FilterInit)
      // Main Loop
      .def("graft_adjacent_nnis_to_dag", &NNIEngine::GraftAdjacentNNIsToDAG)
      .def("filter_pre_update", &NNIEngine::FilterPreUpdate)
      .def("filter_eval_adjacent_nnis", &NNIEngine::FilterEvaluateAdjacentNNIs)
      .def("filter_post_update", &NNIEngine::FilterPostUpdate)
      .def("filter_process_adjacent_nnis", &NNIEngine::FilterProcessAdjacentNNIs)
      .def("remove_all_graft_nnis_from_dag", &NNIEngine::RemoveAllGraftedNNIsFromDAG)
      .def("add_accepted_nnis_to_dag", &NNIEngine::AddAcceptedNNIsToDAG)
      // Post Loop
      .def("update_adjacent_nnis", &NNIEngine::UpdateAdjacentNNIs)
      .def("update_accepted_nnis", &NNIEngine::UpdateAcceptedNNIs)
      .def("update_rejected_nnis", &NNIEngine::UpdateRejectedNNIs)
      .def("update_scored_nnis", &NNIEngine::UpdateScoredNNIs)
      // Filtering schemes
      .def("set_no_filter", &NNIEngine::SetNoFilter,
           "Set filter to either accept (True) or deny (False) all NNIs.",
           py::arg("set_all_nni_to_accept"))
      .def("set_gp_likelihood_cutoff_filtering_scheme",
           &NNIEngine::SetGPLikelihoodCutoffFilteringScheme,
           "Set filtering scheme to use Generalized Pruning based on constant score "
           "cutoff.")
      .def("set_tp_likelihood_cutoff_filtering_scheme",
           &NNIEngine::SetTPLikelihoodCutoffFilteringScheme,
           "Set filtering scheme to use Top Pruning with Likelihoods based on constant "
           "score cutoff.")
      .def("set_tp_parsimony_cutoff_filtering_scheme",
           &NNIEngine::SetTPParsimonyCutoffFilteringScheme,
           "Set filtering scheme to use Top Pruning with Parsimony based on constant "
           "score cutoff.")
      .def("set_gp_likelihood_drop_filtering_scheme",
           &NNIEngine::SetGPLikelihoodDropFilteringScheme,
           "Set filtering scheme to use Generalized Pruning based on based on drop "
           "from best score.")
      .def("set_tp_likelihood_drop_filtering_scheme",
           &NNIEngine::SetTPLikelihoodDropFilteringScheme,
           "Set filtering scheme to use Top Pruning with Likelihoods based on drop "
           "from best score.")
      .def("set_tp_parsimony_drop_filtering_scheme",
           &NNIEngine::SetTPParsimonyDropFilteringScheme,
           "Set filtering scheme to use Top Pruning with Parsimony based on drop from "
           "best score.")
      .def("set_top_n_score_filtering_scheme", &NNIEngine::SetTopNScoreFilteringScheme,
           "Set filter scheme that accepts the top N best-scoring NNIs.",
           py::arg("top_n"), py::arg("max_is_best") = true)
      // Options
      .def("set_include_rootsplits", &NNIEngine::SetIncludeRootsplitNNIs,
           "Set whether to include rootsplits in adjacent NNIs")
      .def("set_reevaluate_rejected_nnis", &NNIEngine::SetReevaluateRejectedNNIs,
           "Set whether to re-evaluate NNIs which were rejected in a previous pass.")
      // Scoring
      .def("get_score_by_nni", &NNIEngine::GetScoreByNNI, "Get score by NNI.")
      .def("get_score_by_edge", &NNIEngine::GetScoreByEdge, "Get score by EdgeId.");

  py::class_<SankoffHandler> parsimony_engine_class(
      m, "parsimony_tree_engine",
      "An engine that computes parsimonies for tree topologies.");
  parsimony_engine_class.def("compute_parsimony",
                             [](SankoffHandler &self, const RootedTree &tree) {
                               self.RunSankoff(tree.Topology());
                               return self.ParsimonyScore();
                             });

  py::class_<FatBeagle> likelihood_engine_class(
      m, "likelihood_tree_engine",
      "An engine that computes likelihoods for tree topologies.");
  likelihood_engine_class.def("compute_likelihood",
                              [](FatBeagle &self, const RootedTree &tree) {
                                return self.UnrootedLogLikelihood(tree);
                              });

  // ** Node Topology

  py::class_<Node::Topology> topology_class(
      m, "node_topology", "A node in a node topology representing a tree.");
  topology_class.def("__str__", [](const Node::Topology &self) { self->Leaves(); })
      .def(
          "id", [](const Node::Topology &self) { return self->Id(); },
          "Unique node id within topology.")
      .def(
          "to_leaves", [](const Node::Topology &self) { return self->Leaves(); },
          "Output node to leave bitset.")
      .def(
          "build_subsplit",
          [](const Node::Topology &self) { return self->BuildSubsplit(); },
          "Build subsplit node bitset of node.")
      .def(
          "build_pcsp",
          [](const Node::Topology &self, const size_t child_id) {
            Assert(child_id < 2, "child_count must be 0 (left) or 1 (right).");
            auto child_clade =
                (child_id == 0) ? SubsplitClade::Left : SubsplitClade::Right;
            return self->BuildPCSP(child_clade);
          },
          "Build PCSP edge bitset of edge below node.")
      .def(
          "build_set_of_subsplits",
          [](const Node::Topology &self) { return self->BuildSetOfSubsplits(); },
          "Build a vector of all subsplit bitsets for all nodes in topology.")
      .def(
          "build_set_of_pcsps",
          [](const Node::Topology &self) { return self->BuildSetOfPCSPs(); },
          "Build vector of all PCSP edge bitsets for all edges in topology.")
      .def(
          "to_newick", [](const Node::Topology &self) { return self->Newick(); },
          "Output to Newick string.");

  // ** Bitsets, Subsplits, PCSPs, NNIs, etc

  py::class_<Bitset> bitset_class(
      m, "bitset", "A bitset representing the taxon membership of a Subsplit or PCSP.");
  bitset_class.def(py::init<const std::string &>())
      .def("__str__", &Bitset::ToString)
      .def("__eq__",
           [](const Bitset &self, const Bitset &other) { return self == other; })
      .def("__hash__", &Bitset::Hash)
      .def("to_string", &Bitset::ToString, "Output to bitset string.")
      .def("clade_get_count", &Bitset::Count)
      .def("subsplit_get_clade",
           [](const Bitset &self, const size_t i) {
             SubsplitClade clade =
                 (i == 0) ? SubsplitClade::Left : SubsplitClade::Right;
             return self.SubsplitGetClade(clade);
           })
      .def("subsplit_to_string", &Bitset::SubsplitToString,
           "Output as Subsplit-style string.")
      .def("pcsp_to_string", &Bitset::PCSPToString, "Output as PCSP-style string.")
      .def("pcsp_get_parent_subsplit", &Bitset::PCSPGetParentSubsplit,
           "Get parent subsplit from PCSP.")
      .def("pcsp_get_child_subsplit", &Bitset::PCSPGetChildSubsplit,
           "Get child subsplit from PCSP.");
  m.def(
      "subsplit",
      [](const std::string &left_clade, const std::string &right_clade) {
        return Bitset::Subsplit(left_clade, right_clade);
      },
      "A Subsplit Bitset constructed from two Bitset Clades.");
  m.def(
      "pcsp",
      [](const Bitset &parent, const Bitset &child) {
        return Bitset::PCSP(parent, child);
      },
      "A PCSP Bitset constructed from two Bitset Subsplits.");

  py::class_<NodeId> node_id_class(m, "node_id",
                                   "An ID representing a unique node within a DAG.");
  node_id_class.def(py::init<const size_t>())
      .def("__str__", [](const NodeId self) { return self.ToString(); })
      .def("__eq__", [](const NodeId lhs, const NodeId rhs) { return lhs == rhs; })
      .def("value", [](const NodeId &self) -> int { return self.value_; });

  py::class_<EdgeId> edge_id_class(m, "edge_id",
                                   "An ID representing a unique edge within a DAG.");
  edge_id_class.def(py::init<const size_t>())
      .def("__str__", [](const EdgeId self) { return self.ToString(); })
      .def("__eq__", [](const EdgeId lhs, const EdgeId rhs) { return lhs == rhs; })
      .def("value", [](const EdgeId &self) -> int { return self.value_; });

  py::class_<TaxonId> taxon_id_class(m, "taxon_id",
                                     "An ID representing a unique edge within a DAG.");
  taxon_id_class.def(py::init<const size_t>())
      .def("__str__", [](const TaxonId self) { return self.ToString(); })
      .def("__eq__", [](const TaxonId lhs, const TaxonId rhs) { return lhs == rhs; })
      .def("value", [](const TaxonId &self) -> int { return self.value_; });

  py::class_<TreeId> tree_id_class(m, "tree_id", "An ID representing a unique tree.");
  tree_id_class.def(py::init<const size_t>())
      .def("__str__", [](const TreeId self) { return self.ToString(); })
      .def("__eq__", [](const TreeId lhs, const TreeId rhs) { return lhs == rhs; })
      .def("value", [](const TreeId &self) -> int { return self.value_; });

  py::class_<NNIOperation> nni_op_class(
      m, "nni_op",
      "A proposed NNI Operation for the DAG. Repesents the PCSP to be added.");
  nni_op_class.def(py::init<const std::string &, const std::string &>())
      .def("__str__", &NNIOperation::ToString)
      .def("__eq__",
           [](const NNIOperation &lhs, const NNIOperation &rhs) { return lhs == rhs; })
      .def("__hash__", &NNIOperation::Hash)
      .def("to_string", &NNIOperation::ToString, "Output to string")
      .def("get_parent", &NNIOperation::GetParent, "Get parent Subsplit of PCSP.")
      .def("get_child", &NNIOperation::GetChild, "Get child Subsplit of PCSP.")
      .def("get_central_edge_pcsp", &NNIOperation::GetCentralEdgePCSP,
           "Get central edge PCSP.")
      .def("is_valid", &NNIOperation::IsValid,
           "Checks that NNI Operation is a valid PCSP.");

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
