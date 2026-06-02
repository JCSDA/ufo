/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/Cal_IceThicknessFromFreeboard.h"

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include "eckit/log/CodeLocation.h"
#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static TransformMaker<Cal_IceThicknessFromFreeboard>
makerCal_IceThicknessFromFreeboard_("Calculate iceThickness from seaIceFreeboard");

Cal_IceThicknessFromFreeboard::Cal_IceThicknessFromFreeboard(
    const Parameters_ &options, const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr), options_(options) {}

// -----------------------------------------------------------------------------

void Cal_IceThicknessFromFreeboard::runTransform(
    const std::vector<bool> &apply) {
  oops::Log::trace()
    << classname_
    << "::runTransform: Calculate Ocean Ice Thickness and errors from ice "
       "freeboard, snow, and density observations "
    << std::endl;

  for (const Variable &var :
       {options_.IceFreeboardVariable.value(),
        options_.WaterDensityVariable.value(),
        options_.SnowDepthVariable.value(),
        options_.SnowDensityVariable.value(),
        options_.IceDensityVariable.value()}) {
    if (!obsdb_.has(var.group(), var.variable())) {
      throw eckit::BadValue(
          "" + classname_ + ":: `" + var.group() + '/' + var.variable() +
              "` must be a variable for this variable transform.",
          Here());
    }
  }
  if (options_.calculateErrors.value()) {
    for (const Variable &var : {options_.IceFreeboardESDVariable.value(),
                                options_.IceDensityESDVariable.value(),
                                options_.SnowDensityESDVariable.value(),
                                options_.SnowDepthESDVariable.value()}) {
      if (!obsdb_.has(var.group(), var.variable())) {
        throw eckit::BadValue("" + classname_ + ":: `" + var.fullName() +
                                  "` must be a variable for "
                                  "this variable transform if calculating errors.",
                              Here());
      }
    }
  }

  // dimension
  const size_t nlocs = obsdb_.nlocs();

  // Set default ObsErrorData and QCflagsData information

  const size_t iErrIceThickness =
      obserr_.varnames().find(options_.IceThicknessVariable.value().variable());
  const size_t iFlagIceThickness =
      flags_.varnames().find(options_.IceThicknessVariable.value().variable());
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (!apply[jobs]) continue;
    obserr_[iErrIceThickness][jobs] = missingValueFloat;
    flags_[iFlagIceThickness][jobs] = 0;
  }

  // compute ice thickness
  std::vector<float> iceFreeboard;
  std::vector<float> waterDensity;
  std::vector<float> snowDepth;
  std::vector<float> snowDensity;
  std::vector<float> iceDensity;

  getObservation(options_.IceFreeboardVariable.value().group(),
                 options_.IceFreeboardVariable.value().variable(), iceFreeboard,
                 true);
  getObservation(options_.WaterDensityVariable.value().group(),
                 options_.WaterDensityVariable.value().variable(), waterDensity,
                 true);
  getObservation(options_.SnowDepthVariable.value().group(),
                 options_.SnowDepthVariable.value().variable(), snowDepth,
                 true);
  getObservation(options_.SnowDensityVariable.value().group(),
                 options_.SnowDensityVariable.value().variable(), snowDensity,
                 true);
  getObservation(options_.IceDensityVariable.value().group(),
                 options_.IceDensityVariable.value().variable(), iceDensity,
                 true);

  std::vector<float> iceThickness(nlocs, missingValueFloat);
  for (size_t loc = 0; loc < nlocs; ++loc) {
    if (!apply[loc]) continue;
    if (iceFreeboard[loc] != missingValueFloat &&
        snowDepth[loc] != missingValueFloat &&
        snowDensity[loc] != missingValueFloat &&
        waterDensity[loc] != missingValueFloat &&
        iceDensity[loc] != missingValueFloat) {
      // fi * rhow + ds * rhos) / (rhow - rhoi)
      const float densityDifference = waterDensity[loc] - iceDensity[loc];
      if (densityDifference < std::numeric_limits<float>::min()) {
        iceThickness[loc] = missingValueFloat;
      } else {
        iceThickness[loc] = (iceFreeboard[loc] * waterDensity[loc] +
                              snowDepth[loc] * snowDensity[loc]) /
                            densityDifference;
      }
    }
  }
  putObservation(options_.IceThicknessVariable.value().variable(),
                  iceThickness, options_.IceThicknessVariable.value().group());

  if (!options_.calculateErrors.value()) {
    oops::Log::trace() << classname_ << "::runTransform: done." << std::endl;
    return;
  }

  std::vector<float> iceFreeboardESD;
  std::vector<float> iceDensityESD;
  std::vector<float> snowDensityESD;
  std::vector<float> snowDepthESD;

  getObservation(options_.IceFreeboardESDVariable.value().group(),
                 options_.IceFreeboardESDVariable.value().variable(),
                 iceFreeboardESD, true);
  getObservation(options_.IceDensityESDVariable.value().group(),
                 options_.IceDensityESDVariable.value().variable(),
                 iceDensityESD, true);
  getObservation(options_.SnowDensityESDVariable.value().group(),
                 options_.SnowDensityESDVariable.value().variable(),
                 snowDensityESD, true);
  getObservation(options_.SnowDepthESDVariable.value().group(),
                 options_.SnowDepthESDVariable.value().variable(), snowDepthESD,
                 true);

  std::vector<float> iceThicknessSystematicESD(nlocs, missingValueFloat);
  std::vector<float> iceThicknessRandomESD(nlocs, missingValueFloat);
  for (size_t loc = 0; loc < nlocs; ++loc) {
    if (!apply[loc]) continue;
    if (iceFreeboard[loc] != missingValueFloat &&
        snowDepth[loc] != missingValueFloat &&
        snowDensity[loc] != missingValueFloat &&
        waterDensity[loc] != missingValueFloat &&
        iceDensity[loc] != missingValueFloat &&
        iceFreeboardESD[loc] != missingValueFloat &&
        iceDensityESD[loc] != missingValueFloat &&
        snowDensityESD[loc] != missingValueFloat) {
      const float densityDifference = waterDensity[loc] - iceDensity[loc];
      float iceThicknessRandomVariance = missingValueFloat;
      float iceThicknessSystematicVariance = missingValueFloat;
      if (densityDifference > std::numeric_limits<float>::min()) {
        // random error
        // sqrt((fi_std * rhow / (rhow - rhoi))**2 +
        //       (rhoi_std * (fi * rhow + ds * rhos) /
        //        (rhow - rhoi)**2)**2)
        iceThicknessRandomVariance =
            std::pow(
                iceFreeboardESD[loc] * waterDensity[loc] / densityDifference,
                2) +
            std::pow(iceDensityESD[loc] *
                         (iceFreeboard[loc] * waterDensity[loc] +
                          snowDepth[loc] * snowDensity[loc]) /
                         (densityDifference * densityDifference),
                     2);
        iceThicknessRandomESD[loc] = std::sqrt(iceThicknessRandomVariance);
        // systematic error
        // sqrt((ds_std * rhos / (rhow - rhoi))**2 +
        //       (rhos_std * ds / (rhow - rhoi))**2)
        iceThicknessSystematicVariance =
            std::pow(snowDepthESD[loc] * snowDensity[loc] / densityDifference,
                     2) +
            std::pow(snowDensityESD[loc] * snowDepth[loc] / densityDifference,
                     2);
        iceThicknessSystematicESD[loc] =
            std::sqrt(iceThicknessSystematicVariance);
        obserr_[iErrIceThickness][loc] = std::sqrt(
            iceThicknessSystematicVariance + iceThicknessRandomVariance);
      }
    }
  }
  putObservation(options_.IceThicknessRandomESDVariable.value().variable(),
                 iceThicknessRandomESD,
                 options_.IceThicknessRandomESDVariable.value().group());
  putObservation(options_.IceThicknessSystematicESDVariable.value().variable(),
                 iceThicknessSystematicESD,
                 options_.IceThicknessSystematicESDVariable.value().group());
  oops::Log::trace() << classname_ << "::runTransform: done." << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
