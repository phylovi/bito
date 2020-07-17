// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_operation.hpp"

GPOperations::PrepForMarginalization GPOperations::PrepForMarginalizationOfOperations(
    const GPOperationVector& operations) {
  return PrepForMarginalizationVisitor(operations).ToPrepForMarginalization();
}

std::ostream& operator<<(std::ostream& os, GPOperation const& operation) {
  std::visit(GPOperationOstream{os}, operation);
  return os;
}

std::ostream& operator<<(std::ostream& os, GPOperationVector const& operation_vector) {
  os << "[" << std::endl;
  for (const auto& operation : operation_vector) {
    os << "  " << operation << "," << std::endl;
  }
  os << "]" << std::endl;
  return os;
}
