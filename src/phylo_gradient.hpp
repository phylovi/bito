// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <map>
#include <vector>
#include "phylo_flags.hpp"

using GradientMap = std::map<std::string, std::vector<double>>;

struct PhyloGradient {
  PhyloGradient() = default;
  PhyloGradient(double log_likelihood, GradientMap &gradient)
      : log_likelihood_(log_likelihood), gradient_(gradient){};

  std::vector<double> &operator[](const PhyloMapkey &key) {
    return gradient_[key.GetKey()];
  }

  double log_likelihood_;
  GradientMap gradient_;

  // Gradient mapkeys
  inline const static std::string site_model_key_ = "site_model";
  inline const static std::string clock_model_key_ = "clock_model";
  inline const static std::string substitution_model_key_ = "substitution_model";
  inline const static std::string substitution_model_rates_key_ =
      SubstitutionModel::rates_key_;
  inline const static std::string substitution_model_frequencies_key_ =
      SubstitutionModel::frequencies_key_;
  inline const static std::string branch_lengths_key_ = "branch_lengths";
  inline const static std::string ratios_root_height_key_ = "ratios_root_height";
};

// Mapkeys for GradientMap
namespace PhyloGradientMapkeys {
// Mapkeys
inline static const auto site_model_ =
    PhyloMapkey("SITE_MODEL", PhyloGradient::site_model_key_);
inline static const auto clock_model_ =
    PhyloMapkey("CLOCK_MODEL", PhyloGradient::clock_model_key_);
inline static const auto substitution_model_ =
    PhyloMapkey("SUBSTITUTION_MODEL", PhyloGradient::substitution_model_key_);
inline static const auto substitution_model_rates_ = PhyloMapkey(
    "SUBSTITUTION_MODEL_RATES", PhyloGradient::substitution_model_rates_key_);
inline static const auto substitution_model_frequencies_ =
    PhyloMapkey("SUBSTITUTION_MODEL_FREQUENCIES",
                PhyloGradient::substitution_model_frequencies_key_);
inline static const auto branch_lengths_ =
    PhyloMapkey("BRANCH_LENGTHS", PhyloGradient::branch_lengths_key_);
inline static const auto ratios_root_height_ =
    PhyloMapkey("RATIOS_ROOT_HEIGHT", PhyloGradient::ratios_root_height_key_);

inline static const auto set_ = PhyloMapkeySet(
    "PhyloModel",
    {site_model_, clock_model_, substitution_model_, substitution_model_rates_,
     substitution_model_frequencies_, branch_lengths_, ratios_root_height_});
};  // namespace PhyloGradientMapkeys
