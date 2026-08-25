/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/usenearestneighbors/UseNearestNeighborsAlgorithmBase.h"

#include "eckit/exception/Exceptions.h"

namespace ufo {

UseNearestNeighborsAlgorithmBase::UseNearestNeighborsAlgorithmBase(
    const UseNearestNeighborsAlgorithmParametersBase& params,
    const ObsFilterData& data, const std::vector<bool>& apply,
    const Variables& filtervars, const ioda::ObsDataVector<int>& flags,
    std::vector<std::vector<bool>>& flagged)
    : data_(data),
      obsdb_(data.obsspace()),
      apply_(apply),
      filtervars_(filtervars),
      flags_(flags),
      flagged_(flagged) {}

void UseNearestNeighborsAlgorithmBase::runAlgorithm(
    const UseNearestNeighborsParameters& options) const {
  const UseNearestNeighborsAlgorithmParametersWrapper& wrapper =
      options.algorithmParameters.value();
  execute(wrapper.name.value(), options);
}

std::unique_ptr<UseNearestNeighborsAlgorithmBase>
UseNearestNeighborsAlgorithmFactory::create(
    const UseNearestNeighborsAlgorithmParametersBase& params,
    const ObsFilterData& data, const std::vector<bool>& apply,
    const Variables& filtervars, const ioda::ObsDataVector<int>& flags,
    std::vector<std::vector<bool>>& flagged) {
  auto& makers = getMakers();
  const std::string& name = params.algorithmName.value();
  auto it = makers.find(name);
  if (it == makers.end())
    throw eckit::BadParameter(
        "UseNearestNeighbors algorithm not found: " + name, Here());
  return it->second->make(params, data, apply, filtervars, flags, flagged);
}

std::unique_ptr<UseNearestNeighborsAlgorithmParametersBase>
UseNearestNeighborsAlgorithmFactory::createParameters(const std::string& name) {
  auto& makers = getMakers();
  auto it = makers.find(name);
  if (it == makers.end())
    throw eckit::BadParameter(
        "Parameters for UseNearestNeighbors algorithm not found: " + name,
        Here());
  return it->second->makeParameters();
}

UseNearestNeighborsAlgorithmFactory::UseNearestNeighborsAlgorithmFactory(
    const std::string& name) {
  auto& makers = getMakers();
  makers[name] = this;
}

void UseNearestNeighborsAlgorithmBase::verifyNearestNeighborIDTypesMatch(
    const Variable &idVariable,
    const std::vector<Variable> &nearestNeighborIDVars) const {
  auto idDtype = data_.obsspace().dtype(idVariable.group(), idVariable.variable());
  for (const auto &nearestNeighborIDVar : nearestNeighborIDVars) {
    auto nearestNeighborDtype = data_.obsspace().dtype(
        nearestNeighborIDVar.group(), nearestNeighborIDVar.variable());
    if (nearestNeighborDtype != idDtype) {
      throw eckit::BadParameter(
          "Nearest neighbor identifier variable '" + nearestNeighborIDVar.fullName() +
          "' has type " + std::to_string(static_cast<int>(nearestNeighborDtype)) +
          " but identifier variable '" + idVariable.fullName() +
          "' has type " + std::to_string(static_cast<int>(idDtype)) + ". " +
          "All nearest neighbor identifier variables must have the same type "
          "as the identifier variable.",
          Here());
    }
  }
}

}  // namespace ufo
