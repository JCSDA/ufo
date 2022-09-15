/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/OceanConversions/OceanDepthToPressure.h"

#include <cmath>

#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static TransformMaker<OceanDepthToPressure>
    makerOceanDepthToPressure_("OceanDepthToPressure");

OceanDepthToPressure::OceanDepthToPressure(
        const Parameters_ &options,
        const ObsFilterData &data,
        const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
        const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr),
      depthvariable_(options.DepthVariable),
      depthgroup_(options.DepthGroup),
      pressurevariable_(options.PressureVariable)
{}

// -----------------------------------------------------------------------------

void OceanDepthToPressure::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << "Retrieve Ocean Pressure from Depth and Latitude" << std::endl;

  // dimension
  const size_t nlocs = obsdb_.nlocs();

  // Get all required data
  std::vector<float> depth;
  std::vector<float> lats(nlocs);
  std::vector<float> pressure(nlocs, missingValueFloat);
  std::vector<float> pressureError(nlocs, 1.0);
  std::vector<float> pressurepge;
  std::vector<int> pressureflags;
  getObservation(depthgroup_,
                 depthvariable_,
                 depth, true);
  // Get latitude
  getObservation("MetaData",
                 "latitude",
                 lats, true);
  getObservation("GrossErrorProbability", depthvariable_,
                 pressurepge);
  getObservation("QCFlags", depthvariable_,
                 pressureflags);

  if (pressurepge.empty()) {
    pressurepge.assign(nlocs, missingValueFloat);
  }
  if (pressureflags.empty()) {
    pressureflags.assign(nlocs, 0);
  }

  // compute pressure as function of depth and latitude
  for (size_t loc = 0; loc < nlocs; ++loc) {
    if (!apply[loc]) continue;
    if (depth[loc] != missingValueFloat) {
      pressure[loc] = gsw_p_from_z_f90(-1.0*depth[loc],
                                       lats[loc]);
    }
  }
  obsdb_.put_db("DerivedObsValue", pressurevariable_, pressure);
  const size_t iv = obserr_.varnames().find(pressurevariable_);
  for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
    if (!apply[jobs])
      continue;
    obserr_[iv][jobs] = pressureError[jobs];
  }

  // copy PGEFinal and QCflags to new pressure
  obsdb_.put_db("GrossErrorProbability", pressurevariable_, pressurepge);
  obsdb_.put_db("QCFlags", pressurevariable_, pressureflags);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
