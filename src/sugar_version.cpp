// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Stores metadata like git commit of current build.

#include "sugar_version.hpp"

#define GIT_HASH "f045cf5"
#define GIT_BRANCH "483-nni-search-gp-branch-lengths-fix"
#define GIT_TAGS ""

Version::git_hash_ = #GIT_HASH;
Version::git_branch_ = #GIT_BRANCH;
Version::git_tags_ = #GIT_TAGS;
