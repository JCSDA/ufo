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

  // Set default ObsError information
  {
    std::vector<float> pressureError(nlocs, 1.0);
    const size_t iv = obserr_.varnames().find(pressurevariable_);
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (!apply[jobs])
        continue;
      obserr_[iv][jobs] = pressureError[jobs];
    }
    putObservation(pressurevariable_, pressureError, "DerivedObsError");
  }

  // Copy across the GrossErrorProbability information to derived variable if present
  {
    std::vector<float> pressurepge;
    getObservation("GrossErrorProbability", depthvariable_,
                   pressurepge);
    if (!pressurepge.empty())
      putObservation(pressurevariable_, pressurepge, "GrossErrorProbability");
  }

  // Copy across QC information to derived variable if present
  {
    std::vector<int> pressureflags;
    getObservation("QCflagsData", depthvariable_,
                   pressureflags);
    if (pressureflags.empty()) {
      pressureflags.assign(nlocs, 0);
    }
    putObservation(pressurevariable_, pressureflags, "QCflagsData");
  }

  // compute pressure as function of depth and latitude
  {
    std::vector<float> depth;
    std::vector<float> lats(nlocs);
    std::vector<float> pressure(nlocs, missingValueFloat);
    getObservation(depthgroup_,
                   depthvariable_,
                   depth, true);
    getObservation("MetaData",
                   "latitude",
                   lats, true);
    for (size_t loc = 0; loc < nlocs; ++loc) {
      if (!apply[loc]) continue;
      if (depth[loc] != missingValueFloat) {
        pressure[loc] = gsw_p_from_z_f90(-1.0*depth[loc],
                                         lats[loc]);
      }
    }
    putObservation(pressurevariable_, pressure, "DerivedObsValue");
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
