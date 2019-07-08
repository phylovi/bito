// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_LIBSBN_HPP_
#define SRC_LIBSBN_HPP_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>
#include "alignment.hpp"
#include "beagle.hpp"
#include "build.hpp"
#include "driver.hpp"
#include "tree.hpp"

namespace py = pybind11;

typedef std::unordered_map<std::string, float> StringFloatMap;
typedef std::unordered_map<std::string, uint32_t> StringUInt32Map;

StringFloatMap StringFloatMapOf(BitsetUInt32Map m) {
  StringFloatMap m_str;
  for (auto iter = m.begin(); iter != m.end(); ++iter) {
    m_str[iter->first.ToString()] = static_cast<float>(iter->second);
  }
  return m_str;
}

StringUInt32Map StringUInt32MapOf(BitsetUInt32Map m) {
  StringUInt32Map m_str;
  for (auto iter = m.begin(); iter != m.end(); ++iter) {
    m_str[iter->first.ToString()] = iter->second;
  }
  return m_str;
}

struct SBNInstance {
  std::string name_;
  Driver driver_;
  TreeCollection::TreeCollectionPtr tree_collection_;
  Alignment alignment_;
  CharIntMap symbol_table_;
  int beagle_instance_;

  explicit SBNInstance(const std::string &name)
      : name_(name),
        symbol_table_(beagle::GetSymbolTable()),
        beagle_instance_(-1)
  {}

  ~SBNInstance() { assert(beagleFinalizeInstance(beagle_instance_) == 0); }

  size_t TreeCount() const { return tree_collection_->TreeCount(); }

  void ReadNewickFile(std::string fname) {
    driver_.Clear();
    tree_collection_ = driver_.ParseFile(fname);
  }

  void ReadFastaFile(std::string fname) { alignment_.ReadFasta(fname); }

  void PrintStatus() {
    std::cout << "Status for instance '" << name_ << "':\n";
    if (tree_collection_) {
      std::cout << TreeCount() << " unique tree topologies loaded.\n";
    } else {
      std::cout << "No trees loaded.\n";
    }
    std::cout << alignment_.Data().size() << " sequences loaded.\n";
  }

  StringUInt32Map RootsplitSupport() {
    return StringUInt32MapOf(RootsplitSupportOf(tree_collection_->Trees()));
  }
  StringUInt32Map SubsplitSupport() {
    return StringUInt32MapOf(SubsplitSupportOf(tree_collection_->Trees()));
  }

  void BeagleCreate() {
    if (beagle_instance_ != -1) {
      std::cerr << "BEAGLE Instance already exists and needs to be freed.\n";
      abort();
    }
    if (TreeCount() == 0) {
      std::cerr << "Load some trees into your SBNInstance on which you wish to "
                   "calculate phylogenetic likelihoods.\n";
      abort();
    }
    if (alignment_.SequenceCount() == 0) {
      std::cerr << "Load an alignment into your SBNInstance on which you wish "
                   "to calculate phylogenetic likelihoods.\n";
      abort();
    }
    int tip_count = tree_collection_->FirstTree()->LeafCount();
    if (tip_count != static_cast<int>(alignment_.SequenceCount())) {
      std::cerr << "The number of tree tips doesn't match the alignment "
                   "sequence count!\n";
      abort();
    }
    BeagleInstanceDetails *return_info = new BeagleInstanceDetails();
    beagle_instance_ = beagle::CreateInstance(
        static_cast<int>(tip_count), static_cast<int>(alignment_.Length()),
        return_info);
    // TODO(erick) do something with return_info?
    // TODO(erick) free return_info?
  }

  void PrepareBeagleInstance() {
    beagle::SetTipStates(beagle_instance_, tree_collection_->TagTaxonMap(),
                         alignment_, symbol_table_);
    std::vector<double> pattern_weights(alignment_.Length(), 1.);
    beagleSetPatternWeights(beagle_instance_, pattern_weights.data());
    // Use uniform rates and weights.
    const double weights[1] = {1.0};
    const double rates[1] = {1.0};
    beagleSetCategoryWeights(beagle_instance_, 0, weights);
    beagleSetCategoryRates(beagle_instance_, rates);
  }

  void SetJCModel() {
    std::vector<double> freqs(4, 0.25);
    beagleSetStateFrequencies(beagle_instance_, 0, freqs.data());
    // an eigen decomposition for the JC69 model
    std::vector<double> evec = {1.0, 2.0, 0.0, 0.5,  1.0, -2.0, 0.5,  0.0,
                                1.0, 2.0, 0.0, -0.5, 1.0, -2.0, -0.5, 0.0};

    std::vector<double> ivec = {0.25,  0.25,   0.25, 0.25, 0.125, -0.125,
                                0.125, -0.125, 0.0,  1.0,  0.0,   -1.0,
                                1.0,   0.0,    -1.0, 0.0};

    std::vector<double> eval = {0.0, -1.3333333333333333, -1.3333333333333333,
                                -1.3333333333333333};
    beagleSetEigenDecomposition(beagle_instance_, 0, evec.data(), ivec.data(),
                                eval.data());
  }

  double TreeLogLikelihood(Tree::TreePtr tree) {
    int next_internal_index = tree->LeafCount();
    std::vector<int> node_indices;
    std::vector<double> branch_lengths;
    std::vector<BeagleOperation> operations;
    tree->Root()->PostOrder(
        [&tree, &next_internal_index, &node_indices, &branch_lengths,
         &operations](const Node *node, const std::vector<int> &below_indices) {
          if (node->IsLeaf()) {
            auto leaf_id = node->MaxLeafID();
            node_indices.push_back(leaf_id);
            branch_lengths.push_back(tree->BranchLength(node));
            return static_cast<int>(leaf_id);
          }
          // else
          int this_index = next_internal_index;
          node_indices.push_back(this_index);
          branch_lengths.push_back(tree->BranchLength(node));
          assert(below_indices.size() == 2);
          BeagleOperation op = {
              this_index,                          // dest
              BEAGLE_OP_NONE,   BEAGLE_OP_NONE,    // src and dest scaling
              below_indices[0], below_indices[0],  // src1 and matrix1
              below_indices[1], below_indices[1]   // src2 and matrix2
          };
          operations.push_back(op);
          next_internal_index++;
          return this_index;
        });
    beagleUpdateTransitionMatrices(beagle_instance_,
                                   0,  // eigenIndex
                                   node_indices.data(),
                                   NULL,  // firstDerivativeIndices
                                   NULL,  // secondDervativeIndices
                                   branch_lengths.data(),
                                   static_cast<int>(branch_lengths.size()));
    beagleUpdatePartials(beagle_instance_,
                         operations.data(),  // eigenIndex
                         static_cast<int>(operations.size()),
                         BEAGLE_OP_NONE);  // cumulative scale index
    double log_like = 0;
    std::vector<int> root_index = {node_indices.back()};
    std::vector<int> category_weight_index = {0};
    std::vector<int> state_frequency_index = {0};
    std::vector<int> cumulative_scale_index = {BEAGLE_OP_NONE};
    int count = 1;
    beagleCalculateRootLogLikelihoods(
        beagle_instance_, root_index.data(), category_weight_index.data(),
        state_frequency_index.data(), cumulative_scale_index.data(), count,
        &log_like);
    return log_like;
  }

  static void f(py::array_t<double> array) {
    py::buffer_info buf = array.request();
    std::cout << "You passed a " << buf.ndim << " dim array" << std::endl;
    double *ptr = reinterpret_cast<double *>(buf.ptr);
    for (auto idx = 0; idx < buf.shape[0]; idx++) {
      std::cout << ptr[idx] << std::endl;
    }
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("libsbn") {
  SBNInstance inst("charlie");
  inst.ReadNewickFile("data/hello.nwk");
  inst.ReadFastaFile("data/hello.fasta");
  inst.BeagleCreate();
  inst.PrepareBeagleInstance();
  inst.SetJCModel();
  CHECK_LT(abs(inst.TreeLogLikelihood(inst.tree_collection_->FirstTree()) -
               -84.852358),
           0.000001);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_LIBSBN_HPP_
