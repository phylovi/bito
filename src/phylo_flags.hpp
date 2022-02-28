// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// PhyloFlags are used for adding optional arguments to functions that are specified by
// the user for functions such as SBNInstance::PhyloGradients.
// PhyloMapkeys contains the keys used for accessing members of the
// GradientMap, the output returned by SBNInstance::PhyloGradients.
//
// For convenience, this should match the output mapkey if that mapkey
// directly corresponds to a function option for the computation of underlying data.
// (e.g. in FatBeagle::Gradient, we have an option to compute
// `substitution_model_rates`. If flag is set, then the return map will contain a
// `substitution_model_rates` key).
//

#pragma once

#include "sugar.hpp"

// ** Phylo Mapkey
// This is the base mapkey, used for enumerating possible keys for a given map.
class PhyloMapkey {
 public:
  PhyloMapkey(const std::string &name, const std::string &key);

  // Comparators
  static int Compare(const PhyloMapkey &mapkey_a, const PhyloMapkey &mapkey_b);
  // General compare.
  bool operator==(const PhyloMapkey &other) const;
  friend bool operator==(const PhyloMapkey &lhs, const PhyloMapkey &rhs);
  bool operator<(const PhyloMapkey &other) const;
  friend bool operator<(const PhyloMapkey &lhs, const PhyloMapkey &rhs);
  // Compare against String
  bool operator==(const std::string &other_name) const;
  bool operator<(const std::string &other_name) const;

  // Getters
  std::string GetName() const { return name_; };
  std::string GetKey() const { return key_; };

 private:
  // This is a descriptive name of the mapkey that will be visible to the user in
  // bito pybind interface.
  std::string name_;
  // This is the uniquely identifiable key that is used for accessing a map location.
  std::string key_;
};

// ** Phylo Mapkey Set
// Contains all possible options for function.
class PhyloMapkeySet {
 public:
  using MapkeyMap = std::map<std::string, PhyloMapkey>;

  explicit PhyloMapkeySet(const std::string &name) : name_(name){};
  PhyloMapkeySet(const std::string &name, const std::vector<PhyloMapkey> &mapkeys);

  // Insert individual mapkey.
  void AddMapkey(const PhyloMapkey &mapkey, const bool overwrite = false);
  // Does mapkey already exist in set?
  bool ContainsMapkey(const PhyloMapkey &mapkey);
  const MapkeyMap &GetAllNames() const { return all_mapkeys_; };
  std::string ToString() const;

 private:
  // Name for mapkey set.
  std::string name_;
  // List of all possible keys.
  MapkeyMap all_mapkeys_;
};

// ** Phylo FlagOption
// This is the base option flag type.  Also specifies default behaviour for flags.
class PhyloFlagOption {
 public:
  enum class FlagType { None, Boolean, SetValue, RunAll };
  enum class DataType { None, Double };

  PhyloFlagOption();
  PhyloFlagOption(const std::string &name, const std::string &flag,
                  const FlagType flag_type, const DataType data_type,
                  const bool is_set_when_running_defaults,
                  const bool is_set_when_not_running_defaults);
  // PhyloFlagOption FlagType-specific constructors.
  static PhyloFlagOption BooleanOption(
      const std::string &name, const std::string &flag,
      const bool is_set_when_running_defaults = true,
      const bool is_set_when_not_running_defaults = false);
  static PhyloFlagOption SetValueOption(const std::string &name,
                                        const std::string &flag,
                                        const DataType data_type);

  // Add Child Flags (these flags are set when Parent is set).
  void AddChild(const PhyloFlagOption &child);
  void AddChild(const std::string child_flag);
  // Output to String.
  std::string ToString() const;
  // Comparators
  static int Compare(const PhyloFlagOption &flag_a, const PhyloFlagOption &flag_b);
  // General compare.
  bool operator==(const PhyloFlagOption &other);
  friend bool operator==(const PhyloFlagOption &lhs, const PhyloFlagOption &rhs);
  bool operator<(const PhyloFlagOption &other);
  friend bool operator<(const PhyloFlagOption &lhs, const PhyloFlagOption &rhs);
  // Compare against String
  bool operator==(const std::string &other_name);
  bool operator<(const std::string &other_name);

  // Getters
  std::string GetName() const { return name_; };
  std::string GetFlag() const { return flag_; };
  std::string operator()() const { return GetFlag(); };
  bool IsSetWhenRunningDefaults() const { return is_set_when_running_defaults_; };
  bool IsSetWhenNotRunningDefaults() const {
    return is_set_when_not_running_defaults_;
  };
  FlagType GetFlagType() const { return flag_type_; };
  DataType GetDataType() const { return data_type_; };
  const StringVector &GetChildFlags() const { return child_flags_; };

 private:
  // This is a descriptive name of the flag option that will be visible to the user in
  // bito pybind interface.
  std::string name_;
  // This is the uniquely identifiable string that is used for setting/adding flag
  // options.
  std::string flag_;
  // Determines default behavior (whether to consider this option set or unset) when
  // `is_run_defaults_` flag is set. This behavior is overridden when by explicit flags.
  bool is_set_when_running_defaults_;
  bool is_set_when_not_running_defaults_;
  // This gives the type of flag. There are:
  // - Boolean: these options are either set or unset.
  // - SetValue: these options have an associated value.
  FlagType flag_type_;
  // This gives the underlying datatype of the flag.
  // Datatype is None if anything other than a SetValue.
  DataType data_type_;
  // These allow for subflags, corresponding to subroutines of given superflag routine.
  // (e.g. in FatBeagle::Gradient, `substitution_model` flag has two subflags,
  // `substitution_model_rates` and `substitution_model_frequencies`. If both subflags
  // are set, we should consider the superflag set as well.)
  StringVector child_flags_;
};

// ** Phylo FlagOption Set
// Contains all possible options for function.
class PhyloFlagOptionSet {
 public:
  using FlagOptionMap = std::map<std::string, PhyloFlagOption>;
  using SubFlagOptionSetMap = std::map<std::string, PhyloFlagOptionSet *>;

  explicit PhyloFlagOptionSet(const std::string &name);

  PhyloFlagOptionSet(const std::string &name,
                     const std::vector<PhyloFlagOption> &options);
  PhyloFlagOptionSet(const std::string &name,
                     const std::vector<PhyloFlagOption> &options,
                     PhyloFlagOptionSet &parent_optionset);
  // Add Flag Option.
  void AddFlagOption(const PhyloFlagOption &option, const bool overwrite = false);
  // Find Flag by name.
  bool ContainsFlagOption(const PhyloFlagOption &option);
  std::optional<const PhyloFlagOption> FindFlagOptionByName(
      const std::string &name) const;
  // Add Option Set for Subroutines.
  void AddSubPhyloFlagOptionSet(PhyloFlagOptionSet &sub_option_set,
                                const bool overwrite = false);
  std::optional<PhyloFlagOptionSet *> FindSubPhyloFlagOptionSet(
      const std::string name) const;

  // Getters
  std::string GetName() const { return name_; };
  const FlagOptionMap &GetOptions() const { return all_options_; };
  const SubFlagOptionSetMap &GetSubOptionsets() const { return sub_optionsets_; };
  // Get all FlagOption name, flag strings.
  StringPairVector GetAllNames(
      std::optional<StringPairVector> vec_to_append = std::nullopt) const;
  std::string ToString() const;

 private:
  // Name for option set.
  std::string name_;
  // List of all possible options user can set.
  // Map of each flag's name to the flag.
  FlagOptionMap all_options_;
  // Option Sets for Subroutines.
  SubFlagOptionSetMap sub_optionsets_;
};

namespace MasterFlagOptions {
// This determines whehter function will run its default behavior.
inline static auto run_defaults_ =
    PhyloFlagOption("RUN_DEFAULTS", "run_defaults", PhyloFlagOption::FlagType::RunAll,
                    PhyloFlagOption::DataType::None, false, false);

inline static auto set_ = PhyloFlagOptionSet("GLOBAL", {run_defaults_});
};  // namespace MasterFlagOptions

// ** Phylo Flags
// User-facing object.  Sets and stores flags for user and resolves flag value when
// function is called.
class PhyloFlags {
 public:
  using FlagMap = std::map<std::string, std::pair<bool, std::optional<double>>>;

  PhyloFlags(bool is_run_defaults = true,
             PhyloFlagOptionSet &optionset = MasterFlagOptions::set_)
      : explicit_flags_(), is_run_defaults_(is_run_defaults), optionset_(&optionset){};

  template <class T>
  PhyloFlags(const std::vector<T> &key_vec, bool is_run_defaults = true,
             PhyloFlagOptionSet &optionset = MasterFlagOptions::set_)
      : explicit_flags_(), is_run_defaults_(is_run_defaults), optionset_(&optionset) {
    for (auto &key : key_vec) {
      SetFlag(key);
    }
  };

  // ** Flag Setter
  // FlagSet functions add or return a boolean and associated value to/from the map.

  // Final SetFlag.
  void SetFlag(const PhyloFlagOption &flag, const bool set = true,
               const double value = 1.0);
  void SetFlag(const PhyloFlagOption &flag, const double value);

  // If passed SetFlag with flag_name string, look up associated PhyloFlagOption flag
  // and forward.
  template <typename... ArgTypes>
  void SetFlag(const std::string &flag_name, ArgTypes... args) {
    // Find Phyloflag.
    std::optional<PhyloFlagOption> flag = optionset_->FindFlagOptionByName(flag_name);
    Assert(flag.has_value(),
           "Attempted to set a option flag by name that does not exist: \"" +
               flag_name + "\"");
    SetFlag(flag.value(), args...);
  }

  // If passed SetFlag with tuples or pairs, unbind and forward.
  template <typename... PairTypes>
  void SetFlag(const std::pair<PairTypes...> pair) {
    std::apply([this](auto &&...args) { return SetFlag(args...); }, pair);
  };
  template <typename... TupleTypes>
  void SetFlag(const std::tuple<TupleTypes...> tuple) {
    std::apply([this](auto &&...args) { return SetFlag(args...); }, tuple);
  };

  // Add in all flags from a vector.
  void SetAllFlags(const StringVector &key_vec);
  // Add in all flags from another PhyloFlags.
  void AddPhyloFlags(const std::optional<PhyloFlags> phylo_flags,
                     const bool overwrite = true);
  // Clear all set flags and values.
  void ClearFlags();

  // ** Flag Checker
  // Determine whether the associated flag will be evaluated as true or false.
  // - (1) Returns the flag's value if it has been explicitly set.
  // - (2) If not, checks whether the `is_run_defaults` flag has been set, in which case
  // the flag's default behavior is returned.
  // - (3) If not, returns false.
  bool IsFlagSet(const PhyloFlagOption &flag) const;
  bool IsFlagNotSet(const PhyloFlagOption &flag) const;
  // Checks if a flag if user may or may not have passed any options.
  // If options have not been passed, uses flag's default behavior.
  static bool IsFlagSet(const std::optional<PhyloFlags> phylo_flags,
                        const PhyloFlagOption &flag);
  static bool IsFlagNotSet(const std::optional<PhyloFlags> phylo_flags,
                           const PhyloFlagOption &flag);

  // ** Flag Value Getter
  // Returns the value associated with the flag.
  std::optional<double> GetFlagValue(const std::string &flag_name) const;
  std::optional<double> GetFlagValue(const PhyloFlagOption &flag) const;
  // Returns the flag's value if set, otherwise returns default value.
  double GetFlagValueIfSet(const std::string &flag_name, double default_value) const;
  double GetFlagValueIfSet(const PhyloFlagOption &flag, double default_value) const;
  static double GetFlagValueIfSet(const std::optional<PhyloFlags> phylo_flags,
                                  const PhyloFlagOption &flag, double default_value);

  // ** "Run Defaults" Flag
  // Special flag that triggers all other flags' default behavior.
  void SetRunDefaultsFlag(bool is_set);
  bool IsRunDefaultsSet() const;

  // ** Optionset

  const PhyloFlagOptionSet &GetOptionSet() const { return *optionset_; }

  // ** Miscellaneous

  // Get Map of all Set Flags.
  const FlagMap &GetFlagMap() const;
  // Get PhyloFlagOptionSet in use.
  const PhyloFlagOptionSet &GetFlagOptionSet() const;
  // Interprets flags as a string.
  std::string ToString() const;

 private:
  // Check if flag option has been explicitly set.
  bool IsFlagInMap(const PhyloFlagOption &flag) const;
  bool IsFlagInMap(const std::string &flag) const;
  // Explictly set flag option by adding to map.
  void AddFlagToMap(const PhyloFlagOption &flag, const bool set = true,
                    const double value = 1.0f);

  // Stores all option flags that have been manually modified, with a bool whether the
  // flag has been set, and an associated data value.
  FlagMap explicit_flags_;
  // This is a special flag that determines behavior if option is not explicitly set.
  // If is_run_defaults_ is false, all flags are treated as if unset.
  // Otherwise, all flags are treated as their default.
  bool is_run_defaults_;
  // Current Options
  PhyloFlagOptionSet *optionset_ = &MasterFlagOptions::set_;
};

// ** FlagOption Sets

// Flag Options for requesting gradients via FatBeagle::Gradient
namespace PhyloGradientFlagOptions {
inline static const auto site_model_ =
    PhyloFlagOption::BooleanOption("SITE_MODEL", "site_model", true);
inline static const auto clock_model_ =
    PhyloFlagOption::BooleanOption("CLOCK_MODEL", "clock_model", true);
inline static const auto ratios_root_height_ =
    PhyloFlagOption::BooleanOption("RATIOS_ROOT_HEIGHT", "ratios_root_height", true);
inline static const auto substitution_model_ =
    PhyloFlagOption::BooleanOption("SUBSTITUTION_MODEL", "substitution_model", true);
inline static const auto include_log_det_jacobian_gradient_ =
    PhyloFlagOption::BooleanOption("INCLUDE_LOG_DET_JACOBIAN_GRADIENT",
                                   "include_log_det_jacobian_gradient", true, true);
inline static const auto use_stickbreaking_transform_ = PhyloFlagOption::BooleanOption(
    "USE_STICKBREAKING_TRANSFORM", "use_stickbreaking_transform", true, true);
inline static const auto set_gradient_delta_ = PhyloFlagOption::SetValueOption(
    "SET_GRADIENT_DELTA", "set_gradient_delta", PhyloFlagOption::DataType::Double);

inline static auto set_ = PhyloFlagOptionSet(
    "SBNInstance::Gradient",
    {site_model_, clock_model_, ratios_root_height_, site_model_, substitution_model_,
     include_log_det_jacobian_gradient_, set_gradient_delta_},
    MasterFlagOptions::set_);
};  // namespace PhyloGradientFlagOptions

// Flag Options for FatBeagle::LogLikelihood
namespace LogLikelihoodFlagOptions {
inline static const auto include_log_det_jacobian_likelihood_ =
    PhyloFlagOption::BooleanOption("INCLUDE_LOG_DET_JACOBIAN_LIKELIHOOD",
                                   "include_log_det_jacobian_likelihood", true, true);

inline static const PhyloFlagOptionSet set_ =
    PhyloFlagOptionSet("SBNInstance::LogLikelihood",
                       {include_log_det_jacobian_likelihood_}, MasterFlagOptions::set_);
};  // namespace LogLikelihoodFlagOptions
