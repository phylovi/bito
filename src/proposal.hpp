#include <vector>
#include "beagle.hpp"

class PhyloModel {
 public:
  PhyloModel(std::unique_ptr<SubstitutionModel> substitution_model,
             std::unique_ptr<SiteModel> site_model,
             std::unique_ptr<ClockModel> clock_model);
};

using PhyloModelPtr = std::shared_ptr<PhyloModel>;

// @M: The most "controversial" part of this design may be that I have connected
// the eigenstuff and the beagle instance. I think this is a good idea because
// the beagle instance stores the eigenstuff in its own way, so it's nice to be
// able to only have to update it in one place and have things get automatically
// taken care of on the beagle side of things.

class FatBeagle {
 public:
  // This constructor makes the beagle_instance_;
  FatBeagle(PhyloModelPtr phylo_model, const SitePattern &site_pattern);
  ~FatBeagle() {
    // TODO Finalize beagle_instance_
  }
  // TODO disable copy and move assignment and constructors;

  // This function updates the eigenstuff in the object and in Beagle.
  // Question for Mathieu-- does this make sense for the clock model too?
  void SetParameters(SubstitutionModel::Parameterization parameterization);
  // TODO Rescaling should probably be set automatically.
  void SetRescaling();

  double LogLikelihood(const Tree &tree);
  std::pair<double, std::vector<double>> BranchGradient(const Tree &tree);

 private:
  void SetTipStates(const SitePattern &site_pattern);

  BeagleInstance beagle_instance_;
  PhyloModelPtr phylo_model_;
  bool rescaling_;
  Eigen::VectorXd frequencies_;
  EigenMatrixXd evec_;
  EigenMatrixXd ivec_;
  Eigen::VectorXd eval_;
  EigenMatrixXd Q_;
  // Mathieu-- I think that we will need some representation of the
  // parameterization of the site and clock model here too?
};

// This should be just like the existing parallelization routine, except now the
// tree_number indexes both the tree and the parameterization.
//
// So, f is expected to call FatBeagle.SetParameters and then LogLikelihood or
// whatever.
template <typename T>
std::vector<T> Parallelize(
    std::function<T(FatBeagle *, const SubstitutionModel::Parameterization &,
                    const Tree &, bool)>
        f,
    const std::vector<FatBeagle> &fat_beagles,
    const std::vector<SubstitutionModel::Parameterization> &parameterizations,
    const TreeCollection &tree_collection);

// This is the sort of _function_ (not method) that we anticipate passing to
// Parallelize.
double LogLikelihood(
    FatBeagle *fat_beagle,
    const SubstitutionModel::Parameterization &parameterization,
    const Tree &in_tree);

class NewEngine {
 public:
  NewEngine(PhyloModelPtr phylo_model, const SitePattern &site_pattern,
            size_t thread_count);
  // Make thread_count fat_beagles and store in fat_beagles_;

  // Re LogLikelihoods below...
  // It's a little ugly to have two vectors that have to be "in sync": the
  // parameterizations and the trees. This version would look like:
  // - sample trees (happens in libsbn)
  // - hand the branch length vector (with indexing info) for those trees to
  // Python
  // - Python sets those branch lengths and also gives us dicts that contain the
  // parameterizations.
  //
  // Another design would look like:
  // - we have a single ParameterizedTree object that has a tree and a
  // parameterization.
  // - we sample trees, and get a collection of such objects
  // - Python sets branch lengths in this object, and assigns to the
  // parameterization
  //
  // This seems like a cleaner design, and also expresses the fact that we can
  // only calculate likelihoods of trees with all their associated parameters.
  std::vector<double> LogLikelihoods(
      const std::vector<SubstitutionModel::Parameterization> &parameterizations,
      const TreeCollection &tree_collection);

  std::vector<std::pair<double, std::vector<double>>> BranchGradients(
      const std::vector<SubstitutionModel::Parameterization> &parameterizations,
      const TreeCollection &tree_collection);

 private:
  std::vector<FatBeagle> fat_beagles_;
};
