// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This contains the meta information for the git commit that build bito.

#pragma once

#include "sugar.hpp"

class MetaData {
 public:
  static std::string version = "bito v0.1";
  static std::string git_commit_ = "";
  static std::string git_branch_ = "";
  static std::string git_tag = "";
  static std::string commit_date = "";
  static std::string build_date = "";

  static std::string GetVersion() { return version_(); }
  static std::string GetGitCommit() { return git_commit_; }
  static std::string GetGitTag() { return git_tag_; }
  static std::string GetGitBranch() { return git_branch_; }
}
