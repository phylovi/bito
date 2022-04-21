// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "phylo_flags.hpp"

// ** Phylo Mapkey

PhyloMapkey::PhyloMapkey(const std::string &name, const std::string &key)
    : name_(name), key_(key){};

int PhyloMapkey::Compare(const PhyloMapkey &mapkey_a, const PhyloMapkey &mapkey_b) {
  // If mapkey_names are equal, so are the flag.
  if (mapkey_a.key_ == mapkey_b.key_) {
    return 0;
  }
  return (mapkey_a.key_ > mapkey_b.key_) ? 1 : -1;
};

bool PhyloMapkey::operator==(const PhyloMapkey &other) const {
  return Compare(*this, other) == 0;
};
bool operator==(const PhyloMapkey &lhs, const PhyloMapkey &rhs) {
  return PhyloMapkey::Compare(lhs, rhs) == 0;
};
bool PhyloMapkey::operator<(const PhyloMapkey &other) const {
  return Compare(*this, other) < 0;
};
bool operator<(const PhyloMapkey &lhs, const PhyloMapkey &rhs) {
  return PhyloMapkey::Compare(lhs, rhs) < 0;
};
// Compare against String
bool PhyloMapkey::operator==(const std::string &other_name) const {
  return key_ == other_name;
};
bool PhyloMapkey::operator<(const std::string &other_name) const {
  return key_ < other_name;
};

// ** Phylo Mapkey Set

PhyloMapkeySet::PhyloMapkeySet(const std::string &name,
                               const std::vector<PhyloMapkey> &mapkeys)
    : name_(name) {
  for (const auto &mapkey : mapkeys) {
    AddMapkey(mapkey);
  }
};

void PhyloMapkeySet::AddMapkey(const PhyloMapkey &mapkey, const bool overwrite) {
  if (!overwrite) {
    Assert(!ContainsMapkey(mapkey),
           "Attempted to insert Mapkey that already exists in MapkeySet, or non-unique "
           "flag name: " +
               mapkey.GetKey());
  }
  all_mapkeys_.insert(std::make_pair(mapkey.GetName(), mapkey));
};

bool PhyloMapkeySet::ContainsMapkey(const PhyloMapkey &mapkey) {
  return all_mapkeys_.find(mapkey.GetName()) != all_mapkeys_.end();
};

std::string PhyloMapkeySet::ToString() const {
  std::stringstream str;
  for (const auto &[name, mapkey] : all_mapkeys_) {
    std::ignore = name;
    str << mapkey.GetName() << " | " << mapkey.GetKey() << std::endl;
  }
  return str.str();
};

// ** Phylo FlagOption

PhyloFlagOption::PhyloFlagOption()
    : is_set_when_running_defaults_(false),
      is_set_when_not_running_defaults_(false),
      flag_type_(FlagType::None),
      data_type_(DataType::None) {}

PhyloFlagOption::PhyloFlagOption(const std::string &name, const std::string &flag,
                                 const FlagType flag_type, const DataType data_type,
                                 const bool is_set_when_running_defaults,
                                 const bool is_set_when_not_running_defaults)
    : name_(name),
      flag_(flag),
      is_set_when_running_defaults_(is_set_when_running_defaults),
      is_set_when_not_running_defaults_(is_set_when_not_running_defaults),
      flag_type_(flag_type),
      data_type_(data_type),
      child_flags_() {}

PhyloFlagOption PhyloFlagOption::BooleanOption(
    const std::string &name, const std::string &flag,
    const bool is_set_when_running_defaults,
    const bool is_set_when_not_running_defaults) {
  return {name,
          flag,
          FlagType::Boolean,
          DataType::None,
          is_set_when_running_defaults,
          is_set_when_not_running_defaults};
}

PhyloFlagOption PhyloFlagOption::SetValueOption(const std::string &name,
                                                const std::string &flag,
                                                const DataType data_type) {
  return {name, flag, FlagType::SetValue, data_type, false, false};
}

void PhyloFlagOption::AddChild(const PhyloFlagOption &child) { AddChild(child.flag_); }

void PhyloFlagOption::AddChild(const std::string child_flag) {
  child_flags_.push_back(child_flag);
}

std::string PhyloFlagOption::ToString() const {
  std::stringstream str;
  str << "{ " << name_ << ": " << flag_ << " }";
  return str.str();
}

int PhyloFlagOption::Compare(const PhyloFlagOption &flag_a,
                             const PhyloFlagOption &flag_b) {
  // If flag_names are equal, so are the flag.
  if (flag_a.flag_ == flag_b.flag_) {
    return 0;
  }
  return (flag_a.flag_ > flag_b.flag_) ? 1 : -1;
}

bool PhyloFlagOption::operator==(const PhyloFlagOption &other) {
  return Compare(*this, other) == 0;
}

bool operator==(const PhyloFlagOption &lhs, const PhyloFlagOption &rhs) {
  return PhyloFlagOption::Compare(lhs, rhs) == 0;
}

bool PhyloFlagOption::operator<(const PhyloFlagOption &other) {
  return Compare(*this, other) < 0;
}

bool operator<(const PhyloFlagOption &lhs, const PhyloFlagOption &rhs) {
  return PhyloFlagOption::Compare(lhs, rhs) < 0;
}

bool PhyloFlagOption::operator==(const std::string &other_name) {
  return flag_ == other_name;
}

bool PhyloFlagOption::operator<(const std::string &other_name) {
  return flag_ < other_name;
}

// ** Phylo FlagOption Set

PhyloFlagOptionSet::PhyloFlagOptionSet(const std::string &name) : name_(name) {
  AddFlagOption(MasterFlagOptions::run_defaults_);
}

PhyloFlagOptionSet::PhyloFlagOptionSet(const std::string &name,
                                       const std::vector<PhyloFlagOption> &options)
    : name_(name) {
  for (const auto &option : options) {
    AddFlagOption(option);
  }
  AddFlagOption(MasterFlagOptions::run_defaults_);
}

PhyloFlagOptionSet::PhyloFlagOptionSet(const std::string &name,
                                       const std::vector<PhyloFlagOption> &options,
                                       PhyloFlagOptionSet &parent_optionset)
    : name_(name) {
  for (const auto &option : options) {
    AddFlagOption(option);
  }
  AddFlagOption(MasterFlagOptions::run_defaults_);
  parent_optionset.AddSubPhyloFlagOptionSet(*this);
}

void PhyloFlagOptionSet::AddFlagOption(const PhyloFlagOption &option,
                                       const bool overwrite) {
  if (!overwrite) {
    Assert(!ContainsFlagOption(option),
           "Attempted to add FlagOption that already exists in FlagOptionSet, or "
           "non-unique flag name: " +
               option.GetFlag());
  }
  all_options_.insert(std::make_pair(option.GetFlag(), option));
}

bool PhyloFlagOptionSet::ContainsFlagOption(const PhyloFlagOption &option) {
  return all_options_.find(option.GetName()) != all_options_.end();
}

std::optional<const PhyloFlagOption> PhyloFlagOptionSet::FindFlagOptionByName(
    const std::string &flag_name) const {
  // Find if exists in given optionset.
  if (all_options_.find(flag_name) != all_options_.end()) {
    return all_options_.at(flag_name);
  }
  // Find if exists in any child optionsets.
  for (const auto &[name, sub_optionset] : sub_optionsets_) {
    std::ignore = name;
    auto sub_option = sub_optionset->FindFlagOptionByName(flag_name);
    if (sub_option.has_value()) {
      return sub_option;
    }
  }
  return std::nullopt;
}

void PhyloFlagOptionSet::AddSubPhyloFlagOptionSet(PhyloFlagOptionSet &sub_optionset,
                                                  const bool overwrite) {
  if (!overwrite) {
    Assert(sub_optionsets_.find(sub_optionset.GetName()) == sub_optionsets_.end(),
           "Attempted to add a PhyloFlagOptionSet under a pre-existing name: " +
               sub_optionset.GetName());
  }
  sub_optionsets_.insert(std::make_pair(sub_optionset.GetName(), &sub_optionset));
}

std::optional<PhyloFlagOptionSet *> PhyloFlagOptionSet::FindSubPhyloFlagOptionSet(
    const std::string name) const {
  auto sub_optionset = sub_optionsets_.find(name);
  if (sub_optionset == sub_optionsets_.end()) {
    return std::nullopt;
  }
  return sub_optionsets_.at(name);
}

StringPairVector PhyloFlagOptionSet::GetAllNames(
    std::optional<StringPairVector> flag_vec) const {
  if (!flag_vec.has_value()) {
    flag_vec = StringPairVector();
  }
  for (const auto &[name, flag] : GetOptions()) {
    std::ignore = name;
    flag_vec->push_back({flag.GetName(), flag.GetFlag()});
  }
  for (const auto &[name, sub_optionset] : GetSubOptionsets()) {
    std::ignore = name;
    sub_optionset->GetAllNames(flag_vec);
  }
  return *flag_vec;
}

std::string PhyloFlagOptionSet::ToString() const {
  std::stringstream str;
  str << "NAME:" << GetName() << std::endl;
  str << "FLAGS:" << std::endl;
  for (const auto &[name, option] : all_options_) {
    std::ignore = name;
    str << option.GetName() << " | " << option.GetFlag() << " | "
        << option.GetChildFlags() << std::endl;
  }
  return str.str();
}

// ** Phylo Flags

void PhyloFlags::ClearFlags() { explicit_flags_.clear(); }

void PhyloFlags::AddPhyloFlags(const std::optional<PhyloFlags> other_flags,
                               const bool overwrite) {
  if (other_flags.has_value()) {
    for (const auto &[name, bool_data] : other_flags->GetFlagMap()) {
      const auto &[set, data] = bool_data;
      // determines if other_flags will not overwrite flags, just supplements it.
      if (overwrite || (!IsFlagInMap(name))) {
        if (data.has_value()) {
          SetFlag(name, set, *data);
        } else {
          SetFlag(name, set);
        }
      }
    }
  }
}

void PhyloFlags::SetFlag(const PhyloFlagOption &flag, const bool is_set,
                         const double value) {
  // Add given flag.
  AddFlagToMap(flag, is_set, value);
  // Add all child flags of given flag.
  for (const auto &child_flag : flag.GetChildFlags()) {
    SetFlag(child_flag, value);
  }
  // If flag being set is the special run_defaults_ flag.
  if (MasterFlagOptions::run_defaults_.GetName() == flag.GetName()) {
    SetRunDefaultsFlag(true);
  }
}

void PhyloFlags::SetFlag(const PhyloFlagOption &flag, const double value) {
  SetFlag(flag, true, value);
}

void PhyloFlags::AddFlagToMap(const PhyloFlagOption &flag, const bool set,
                              const double value) {
  explicit_flags_.insert(std::make_pair(flag.GetFlag(), std::make_pair(set, value)));
}

void PhyloFlags::SetRunDefaultsFlag(bool is_set) { is_run_defaults_ = is_set; }

bool PhyloFlags::IsRunDefaultsSet() const { return is_run_defaults_; }

bool PhyloFlags::IsFlagInMap(const PhyloFlagOption &flag) const {
  return IsFlagInMap(flag.GetFlag());
}

bool PhyloFlags::IsFlagInMap(const std::string &flag_name) const {
  return (explicit_flags_.find(flag_name) != explicit_flags_.end());
}

std::optional<double> PhyloFlags::GetFlagValue(const PhyloFlagOption &flag) const {
  Assert(flag.GetFlagType() == PhyloFlagOption::FlagType::SetValue,
         "Requested FlagOption value from flag type that does not store associated "
         "value.");
  return GetFlagValue(flag.GetFlag());
}

std::optional<double> PhyloFlags::GetFlagValue(const std::string &flag_name) const {
  if (IsFlagInMap(flag_name)) {
    const auto &[set, value] = explicit_flags_.at(flag_name);
    std::ignore = set;
    return value;
  }
  return std::nullopt;
}

// Returns the value of the flag if set, otherwise returns default value.
double PhyloFlags::GetFlagValueIfSet(const std::string &flag_name,
                                     const double default_value) const {
  auto opt_value = GetFlagValue(flag_name);
  if (opt_value.has_value()) {
    return *opt_value;
  }
  return default_value;
}
double PhyloFlags::GetFlagValueIfSet(const PhyloFlagOption &flag,
                                     const double default_value) const {
  return GetFlagValueIfSet(flag.GetFlag(), default_value);
}
double PhyloFlags::GetFlagValueIfSet(const std::optional<PhyloFlags> phylo_flags,
                                     const PhyloFlagOption &flag,
                                     double default_value) {
  if (phylo_flags.has_value()) {
    return phylo_flags->GetFlagValueIfSet(flag, default_value);
  }
  return default_value;
}

const PhyloFlags::FlagMap &PhyloFlags::GetFlagMap() const { return explicit_flags_; }

const PhyloFlagOptionSet &PhyloFlags::GetFlagOptionSet() const { return *optionset_; }

std::string PhyloFlags::ToString() const {
  std::ostringstream rep;
  rep << "{ ";
  rep << "(DEFAULT: " << IsRunDefaultsSet() << "), ";
  for (const auto &[key, value] : explicit_flags_) {
    rep << "(" << key << ": " << value.first << "), ";
  }
  rep << "}";
  return rep.str();
}

bool PhyloFlags::IsFlagSet(const PhyloFlagOption &flag) const {
  // (1) Priority is given to explicitly set options.
  if (IsFlagInMap(flag)) {
    const auto &[set, value] = GetFlagMap().at(flag.GetFlag());
    std::ignore = value;
    return set;
  }
  // (2) If is_run_default_ option is set, use given individual flag's defined default
  // behavior.
  if (is_run_defaults_) {
    return flag.IsSetWhenRunningDefaults();
  }
  // (3) Otherwise, use flag type-based's default behavior.
  return flag.IsSetWhenNotRunningDefaults();
}

bool PhyloFlags::IsFlagNotSet(const PhyloFlagOption &flag) const {
  return !IsFlagSet(flag);
}

bool PhyloFlags::IsFlagSet(const std::optional<PhyloFlags> phylo_flags,
                           const PhyloFlagOption &flag) {
  // (1) If user has not passed any flags, then fall back to default behavior.
  if (!phylo_flags.has_value()) {
    return flag.IsSetWhenRunningDefaults();
  }
  // (2) If user passed flags, then check if option is set.
  return phylo_flags->IsFlagSet(flag);
}

bool PhyloFlags::IsFlagNotSet(const std::optional<PhyloFlags> phylo_flags,
                              const PhyloFlagOption &flag) {
  return !PhyloFlags::IsFlagSet(phylo_flags, flag);
}
