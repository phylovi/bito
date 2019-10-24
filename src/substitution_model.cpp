#include "substitution_model.hpp"

std::unique_ptr<SubstitutionModel> SubstitutionModel::Create(
    SubstitutionModelType choice) {
  switch (choice) {
    case SubstitutionModelType::GTR:
      return std::make_unique<GTRModel>();
      break;
    default:
      throw std::runtime_error("oh no");
  }
}

