// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BEAGLE_HPP_
#define SRC_BEAGLE_HPP_

#include <functional>
#include <queue>
#include <string>
#include <utility>
#include <vector>
#include "alignment.hpp"
#include "intpack.hpp"
#include "libhmsbeagle/beagle.h"
#include "site_pattern.hpp"
#include "sugar.hpp"
#include "task_processor.hpp"
#include "tree_collection.hpp"

using BeagleInstance = int;

namespace beagle {

int CreateInstance(int tip_count, int alignment_length,
                   BeagleInstanceDetails *return_info);
BeagleInstance CreateInstance(const SitePattern &site_pattern);

void PrepareBeagleInstance(const BeagleInstance beagle_instance,
                           const SitePattern &site_pattern);

void SetJCModel(BeagleInstance beagle_instance);
}  // namespace beagle

// The tests are in libsbn.hpp, where we have access to tree parsing.

#endif  // SRC_BEAGLE_HPP_
