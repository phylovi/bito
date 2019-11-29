// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_FILESYSTEM_HPP_
#define SRC_FILESYSTEM_HPP_

#include <filesystem>

bool is_legal_path(std::string file_path) {
  std::filesystem::path fs_path{file_path};
  if (std::filesystem::exists(fs_path)) {
    return true;
  }
  if (std::filesystem::exists(fs_path.parent_path())) {
    return true;
  }
  return false;
}

#endif  // SRC_FILESYSTEM_HPP_
