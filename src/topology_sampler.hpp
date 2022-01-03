// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "node.hpp"
#include "mersenne_twister.hpp"
#include "eigen_sugar.hpp"

class TopologySampler {
 public:
  class Input {
   public:
    virtual EigenConstVectorXdRef SBNParameters() const = 0;
    virtual size_t RootsplitCount() const = 0;
    virtual const Bitset &RootsplitsAt(size_t rootsplit_idx) const = 0;
    virtual const SizePair &ParentToRangeAt(const Bitset &parent) const = 0;
    virtual const Bitset &IndexToChildAt(size_t child_idx) const = 0;
  };

  Node::NodePtr SampleTopology(const Input &input, bool is_rooted) const;
  inline void SetSeed(uint64_t seed) { mersenne_twister_.SetSeed(seed); }

 private:
  using Range = std::pair<size_t, size_t>;

  size_t SampleIndex(const Input &input, Range range) const;
  Node::NodePtr SampleTopology(const Input &input, const Bitset &parent_subsplit) const;

  MersenneTwister mersenne_twister_;
};

template <typename T>
class SBNInstanceSamplerInput : public TopologySampler::Input {
 public:
  SBNInstanceSamplerInput(const T &inst) : inst_{inst} {}

  size_t RootsplitCount() const override { return inst_.SBNSupport().RootsplitCount(); }

  EigenConstVectorXdRef SBNParameters() const override { return inst_.SBNParameters(); }

  const Bitset &RootsplitsAt(size_t rootsplit_idx) const override {
    return inst_.SBNSupport().RootsplitsAt(rootsplit_idx);
  }

  const SizePair &ParentToRangeAt(const Bitset &parent) const override {
    return inst_.SBNSupport().ParentToRangeAt(parent);
  }

  const Bitset &IndexToChildAt(size_t child_idx) const override {
    return inst_.SBNSupport().IndexToChildAt(child_idx);
  }

 private:
  const T &inst_;
};

class SubsplitDAG;
class SubsplitDAGSamplerInput : public TopologySampler::Input {
 public:
  SubsplitDAGSamplerInput(const SubsplitDAG &dag, EigenConstVectorXdRef sbn_parameters)
      : dag_{dag}, sbn_parameters_{sbn_parameters} {}

  size_t RootsplitCount() const override;
  EigenConstVectorXdRef SBNParameters() const override;
  const Bitset &RootsplitsAt(size_t rootsplit_idx) const override;
  const SizePair &ParentToRangeAt(const Bitset &parent) const override;
  const Bitset &IndexToChildAt(size_t child_idx) const override;

 private:
  const SubsplitDAG &dag_;
  EigenConstVectorXdRef sbn_parameters_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

#include "doctest_constants.hpp"
#include "unrooted_sbn_instance.hpp"
#include "rooted_sbn_instance.hpp"
#include "subsplit_dag.hpp"

template <typename TSBNInstance, typename TIRCount>
static void TestSBNInstance(const TopologySampler& sampler, const TSBNInstance& inst, TIRCount&& count) {
  // Count the frequencies of rooted trees in a file.
  size_t rooted_tree_count_from_file = 0;
  RootedIndexerRepresentationSizeDict counter_from_file(0);
  for (const auto &indexer_representation : inst.MakeIndexerRepresentations()) {
    RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(counter_from_file,
                                                                indexer_representation);
    rooted_tree_count_from_file += count(indexer_representation);
  }
  // Count the frequencies of trees when we sample after training with
  // SimpleAverage.
  size_t sampled_tree_count = 1'000'000;
  RootedIndexerRepresentationSizeDict counter_from_sampling(0);
  ProgressBar progress_bar(sampled_tree_count / 1000);
  SBNInstanceSamplerInput<TSBNInstance> input{inst};
  for (size_t sample_idx = 0; sample_idx < sampled_tree_count; ++sample_idx) {
    const auto rooted_topology = sampler.SampleTopology(input, true);
    RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(
        counter_from_sampling,
        RootedSBNMaps::IndexerRepresentationOf(inst.SBNSupport().Indexer(),
                                               rooted_topology, out_of_sample_index));
    if (sample_idx % 1000 == 0) {
      ++progress_bar;
      progress_bar.display();
    }
  }
  // These should be equal in the limit when we're training with SA.
  for (const auto &[key, _] : counter_from_file) {
    std::ignore = _;
    double observed =
        static_cast<double>(counter_from_sampling.at(key)) / sampled_tree_count;
    double expected =
        static_cast<double>(counter_from_file.at(key)) / rooted_tree_count_from_file;
    CHECK_LT(fabs(observed - expected), 5e-3);
  }
  progress_bar.done();
}

TEST_CASE("TopologySampler: UnrootedSBNInstance tree sampling") {
  UnrootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/five_taxon_unrooted.nwk");
  inst.ProcessLoadedTrees();
  inst.TrainSimpleAverage();

  TestSBNInstance(TopologySampler{}, inst, [](auto&& indexer_representation){
    return indexer_representation.size();
  });
}

TEST_CASE("TopologySampler: RootedSBNInstance tree sampling") {
  RootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/five_taxon_rooted.nwk");
  inst.ProcessLoadedTrees();
  inst.TrainSimpleAverage();
  
  TestSBNInstance(TopologySampler{}, inst, [](auto&&){
    return 1;
  });
}

static void TestSubsplitDAG(const TopologySampler& sampler, const SubsplitDAG& dag,
    std::vector<RootedIndexerRepresentation>&& indexer_representations,
    EigenConstVectorXdRef normalized_sbn_parameters, std::map<SizeVector, double>&& expected_freq) {

  RootedIndexerRepresentationSizeDict counter_from_file(0);
  for (const auto &indexer_representation : indexer_representations) {
    RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(counter_from_file,
                                                                indexer_representation);
  }
  // Count the frequencies of trees when we sample after training with
  // SimpleAverage.
  size_t sampled_tree_count = 1'000'000;
  RootedIndexerRepresentationSizeDict counter_from_sampling(0);
  ProgressBar progress_bar(sampled_tree_count / 1000);
  auto indexer = dag.BuildGPCSPIndexer();
  SubsplitDAGSamplerInput input{dag, normalized_sbn_parameters};
  for (size_t sample_idx = 0; sample_idx < sampled_tree_count; ++sample_idx) {
    const auto rooted_topology = sampler.SampleTopology(input, true);
    RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(
        counter_from_sampling,
        RootedSBNMaps::IndexerRepresentationOf(indexer,
                                               rooted_topology, out_of_sample_index));
    if (sample_idx % 1000 == 0) {
      ++progress_bar;
      progress_bar.display();
    }
  }
  // These should be equal in the limit when we're training with SA.
  for (const auto &[key, _] : counter_from_file) {
    std::ignore = _;
    double observed =
        static_cast<double>(counter_from_sampling.at(key)) / sampled_tree_count;
    double expected = expected_freq.at(key);
    CHECK_LT(fabs(observed - expected), 5e-3);
  }
  progress_bar.done();
}

static std::vector<RootedIndexerRepresentation> MakeIndexerRepresentations(const SubsplitDAG& dag, const RootedTreeCollection& tree_collection) {
  auto indexer = dag.BuildGPCSPIndexer();
  std::vector<RootedIndexerRepresentation> representations;
  representations.reserve(tree_collection.trees_.size());
  for (const auto &tree : tree_collection.trees_) {
    representations.push_back(dag.IndexerRepresentationOf(indexer, tree.Topology(), out_of_sample_index));
  }
  return representations;
}

TEST_CASE("TopologySampler: SubsplitDAG::BuildUniformOnTopologicalSupportPrior tree sampling") {
  Driver driver;
  auto tree_collection =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile("data/five_taxon_rooted.nwk"));
  SubsplitDAG dag(tree_collection);
  TestSubsplitDAG(TopologySampler{}, dag, MakeIndexerRepresentations(dag, tree_collection),
    dag.BuildUniformOnTopologicalSupportPrior(), {
      {{1, 6, 7, 8}, 0.305},
      {{2, 9, 11, 13}, 0.195},
      {{2, 10, 12, 14}, 0.195},
      {{0, 3, 4, 5}, 0.305}
    });
}

TEST_CASE("TopologySampler: SubsplitDAG::BuildUniformOnAllTopologiesPrior tree sampling") {
  Driver driver;
  auto tree_collection =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile("data/five_taxon_rooted.nwk"));
  SubsplitDAG dag(tree_collection);
  TestSubsplitDAG(TopologySampler{}, dag, MakeIndexerRepresentations(dag, tree_collection),
    dag.BuildUniformOnAllTopologiesPrior(), {
      {{1, 6, 7, 8}, 0.345},
      {{2, 9, 11, 13}, 0.173},
      {{2, 10, 12, 14}, 0.173},
      {{0, 3, 4, 5}, 0.309}
    });
}

#endif  // DOCTEST_LIBRARY_INCLUDED