/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/UseNearestNeighbors.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"
#include "ufo/usenearestneighbors/UseNearestNeighborsAlgorithmBase.h"

namespace ufo {

UseNearestNeighbors::UseNearestNeighbors(
    ioda::ObsSpace& obsdb, const Parameters_& parameters,
    ioda::ObsDataVector<int> & flags,
    ioda::ObsDataVector<float> & obserr)
    : FilterBase(obsdb, parameters, flags, obserr), options_(parameters) {
  oops::Log::trace() << "UseNearestNeighbors constructor start" << std::endl;
  oops::Log::trace() << "UseNearestNeighbors constructor complete" << std::endl;
}

UseNearestNeighbors::~UseNearestNeighbors() {
  oops::Log::trace() << "UseNearestNeighbors destructor" << std::endl;
}

void UseNearestNeighbors::applyFilter(
    const std::vector<bool>& apply, const Variables& filtervars,
    std::vector<std::vector<bool>>& flagged) const {
  oops::Log::trace() << "UseNearestNeighbors applyFilter start" << std::endl;

  const UseNearestNeighborsAlgorithmParametersWrapper& algorithmParams =
      options_.algorithmParameters.value();

  std::unique_ptr<UseNearestNeighborsAlgorithmBase> alg =
      UseNearestNeighborsAlgorithmFactory::create(algorithmParams.name.value(),
                                                  data_, apply, filtervars,
                                                  flags_, flagged);

  alg->runAlgorithm(options_);

  oops::Log::trace() << "UseNearestNeighbors applyFilter complete" << std::endl;
}

void UseNearestNeighbors::print(std::ostream& os) const {
  os << "UseNearestNeighbors: config = " << options_ << std::endl;
}

}  // namespace ufo
