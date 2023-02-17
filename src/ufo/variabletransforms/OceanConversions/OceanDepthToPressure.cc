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
{
    if (!obserr_.varnames().has(pressurevariable_) ||
        !flags_.varnames().has(pressurevariable_)) {
        throw eckit::BadValue("`" + pressurevariable_ +
            "` must be an observed or derived variable for the " +
            "`OceanDepthToPressure` variable transform.", Here());
    }
}

// -----------------------------------------------------------------------------

void OceanDepthToPressure::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << "Retrieve Ocean Pressure from Depth and Latitude" << std::endl;

  // dimension
  const size_t nlocs = obsdb_.nlocs();

  // Set ObsErrorData and QCflagsData
  {
    const size_t iErrPres = obserr_.varnames().find(pressurevariable_);
    const size_t iFlagPres = flags_.varnames().find(pressurevariable_);
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (!apply[jobs])
        continue;
      obserr_[iErrPres][jobs] = 1;
      flags_[iFlagPres][jobs] = 0;
    }
  }

  // Copy across the GrossErrorProbability information to derived variable if present
  {
    std::vector<float> pressurepge;
    getObservation("GrossErrorProbability", depthvariable_,
                   pressurepge);
    if (!pressurepge.empty())
      putObservation(pressurevariable_, pressurepge, "GrossErrorProbability");
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
