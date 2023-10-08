/* * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/filters/ObsAccessor.h"
#include "ufo/filters/SpikeAndStepCheck.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// SpikeAndStepCheck: Flag observations that are likely erroneous because they constitute a spike
///  or step in e.g. profiles of ocean temperature vs. depth.

SpikeAndStepCheck::SpikeAndStepCheck(
        ioda::ObsSpace & obsdb,
        const Parameters_ & parameters,
        std::shared_ptr<ioda::ObsDataVector<int> > flags,
        std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "SpikeAndStepCheck constructor" << std::endl;
  // Detect incompatible parameter options:
  validateParameters();
}

// -----------------------------------------------------------------------------

SpikeAndStepCheck::~SpikeAndStepCheck() {
  oops::Log::trace() << "SpikeAndStepCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/// Detect incompatible parameter options.
void SpikeAndStepCheck::validateParameters() const {
  std::vector<float> tolBounds;
  if (((parameters_.toleranceOptions.value().toleranceBoundaries.value() != boost::none &&
        parameters_.toleranceOptions.value().toleranceFactors.value() != boost::none) &&
       (parameters_.toleranceOptions.value().toleranceBoundaries.value().value().size() !=
        parameters_.toleranceOptions.value().toleranceFactors.value().value().size())) ||
      (parameters_.toleranceOptions.value().toleranceBoundaries.value() != boost::none &&
       parameters_.toleranceOptions.value().toleranceFactors.value() == boost::none) ||
      (parameters_.toleranceOptions.value().toleranceBoundaries.value() == boost::none &&
       parameters_.toleranceOptions.value().toleranceFactors.value() != boost::none)) {
    throw eckit::UserError(
        ": 'tolerance factors' and 'tolerance boundaries' must be "
        "float arrays of equal size, or both none.", Here());
  }

  if (parameters_.toleranceOptions.value().toleranceBoundaries.value() != boost::none) {
    tolBounds = parameters_.toleranceOptions.value().toleranceBoundaries.value().value();
    if (!std::is_sorted(tolBounds.begin(), tolBounds.end())) {
      throw eckit::UserError(
        ": 'tolerance boundaries' must start from the lowest x upwards, "
        "regardless of how observations are sorted. Make sure the points "
        "it defines with 'tolerance factors' match up.", Here());
    }
  }

  if (parameters_.boundaryOptions.value() != boost::none) {
    if (parameters_.boundaryOptions.value().value().boundaryRange.value().size() != 2) {
      throw eckit::UserError(
        ": 'boundary layer range' must be a 2-element float array specifying "
        "the x-value at the top and bottom of the boundary layer.", Here());
    }
    if (parameters_.boundaryOptions.value().value().stepTolRange.value().size() != 2) {
      throw eckit::UserError(
        ": 'step tolerance range' must be a 2-element float array specifying "
        "the tolerance factor override [min, max] for x within the boundary layer.", Here());
    }
    if (parameters_.boundaryOptions.value().value().maxDx.value().size() != 2) {
      throw eckit::UserError(
        ": 'maximum x interval' must be a 2-element float array specifying "
        "the maximum dx when x [within, outside] the boundary layer, "
        "where if dx is greater, the filter doesn't check these y-values.", Here());
    }
  }
}

// -----------------------------------------------------------------------------
/// Apply the spike and step check filter.

void SpikeAndStepCheck::applyFilter(const std::vector<bool> & apply,
                                    const Variables & filtervars,
                                    std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "SpikeAndStepCheck filter" << std::endl;

  ObsAccessor obsAccessor =
    ObsAccessor::toObservationsSplitIntoIndependentGroupsByRecordId(obsdb_);
  const size_t totalNumObs = obsAccessor.totalNumObservations();

  const std::string yVarName = parameters_.yVar.value().variable();
  const std::string xVarName = parameters_.xVar.value().variable();
  const std::vector<float> y =
    obsAccessor.getFloatVariableFromObsSpace(parameters_.yVar.value().group(),
                                             yVarName);
  const std::vector<float> x =
    obsAccessor.getFloatVariableFromObsSpace(parameters_.xVar.value().group(),
                                             xVarName);

  std::vector<bool> isThinned(totalNumObs, false);
  std::vector<bool> spikeFlag(totalNumObs, false);
  std::vector<bool> stepFlag(totalNumObs, false);
  // get diagnostic flags from ObsSpace (can't create them in the code, must do in YAML...
  //  ...must be there when ObsSpace constructed)
  if (obsdb_.has("DiagnosticFlags/spike", yVarName)) {
    obsdb_.get_db("DiagnosticFlags/spike", yVarName, spikeFlag);
  } else {
    throw eckit::UserError("Variable 'DiagnosticFlags/spike/" + yVarName + "' does not exist yet. "
                           "It needs to be set up with the 'Create Diagnostic Flags' filter "
                           "prior to using the 'set' or 'unset' action.");
  }
  if (obsdb_.has("DiagnosticFlags/step", yVarName)) {
    obsdb_.get_db("DiagnosticFlags/step", yVarName, stepFlag);
  } else {
    throw eckit::UserError("Variable 'DiagnosticFlags/step/" + yVarName + "' does not exist yet. "
                           "It needs to be set up with the 'Create Diagnostic Flags' filter "
                           "prior to using the 'set' or 'unset' action.");
  }
  // Get vector of all record numbers:
  const std::vector<size_t> & record_numbers = obsdb_.recidx_all_recnums();
  // Loop over record numbers (i.e. profile to profile):
  for (size_t iProfile : record_numbers) {
    // Get non-missing where-included obs indices for this profile
    // (candidateForRetentionIfAnyFilterVariablesPassedQC false - i.e. unselect location
    //  if any filter variable fails QC):
    const std::vector<size_t> obs_indices = obsAccessor.getValidObsIdsInProfile(iProfile,
                                                                                apply,
                                                                                *flags_,
                                                                                filtervars,
                                                                                false);
    // Struct of y, x, dy, dx, dy/dx:
    xyStruct xy;

    // set ydiff and dydx for this record:
    set_xyrec(x,
              y,
              obs_indices,
              xy,
              parameters_);

    // set vector of tolerance values for each x-value:
    std::vector<float> tolerances = set_tolerances(xy.xrec,
                                                   obs_indices,
                                                   parameters_);

    // set flags for obs that fulfil spike/step conditions:
    identifyThinnedObservations(isThinned,
                                spikeFlag,
                                stepFlag,
                                xy,
                                tolerances,
                                obs_indices,
                                parameters_);
  }  // for each record
  obsAccessor.flagRejectedObservations(isThinned, flagged);
  obsdb_.put_db("DiagnosticFlags/spike", yVarName, spikeFlag);
  obsdb_.put_db("DiagnosticFlags/step", yVarName, stepFlag);
}

// -----------------------------------------------------------------------------
/// Given x (independent variable) and y (dependent variable), set x, y, dx, dy and dy/dx
///  for the given record (group).

void SpikeAndStepCheck::set_xyrec(const std::vector<float> &x,
                                  const std::vector<float> &y,
                                  const std::vector<size_t> &obs_indices,
                                  xyStruct &xy,
                                  const Parameters_ &parameters_) const {
  // return a or b, whichever has the larger absolute value, keeping the sign of a:
  auto max_abs = [](const float a, const float b) {
    if (std::abs(a) >= std::abs(b)) {
      return a;
    } else {
      if (a < 0) {
        return -std::abs(b);
      } else {
        return std::abs(b);
      }
    }
  };
  const float missingValueFloat = util::missingValue<float>();
  const float dxResolution = parameters_.toleranceOptions.value().dxResolution.value();
  std::vector<float> boundaryRange;
  std::vector<float> maxDx;
  if (parameters_.boundaryOptions.value() == boost::none) {
    boundaryRange = {0, 0};
    maxDx = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
  } else {
    boundaryRange = parameters_.boundaryOptions.value().value().boundaryRange.value();
    maxDx = parameters_.boundaryOptions.value().value().maxDx.value();
  }
  for (size_t ind = 0; ind < obs_indices.size(); ++ind) {
    xy.yrec.push_back(y[obs_indices[ind]]);
    xy.xrec.push_back(x[obs_indices[ind]]);
    if (ind > 0) {
      xy.ydiff.push_back(missingValueFloat);  // defaults
      xy.xdiff.push_back(missingValueFloat);
      xy.dydx.push_back(missingValueFloat);
      if (xy.xrec[ind] != missingValueFloat && xy.xrec[ind-1] != missingValueFloat &&
          xy.yrec[ind] != missingValueFloat && xy.yrec[ind-1] != missingValueFloat) {
        // Set ydiff, xdiff, dydx only if dx <= maxDx[0] if within boundary layer,
        //                    and only if dx <= maxDx[1] if outside the boundary layer:
        if ((xy.xrec[ind] >= std::min(boundaryRange[0], boundaryRange[1]) &&
            xy.xrec[ind] < std::max(boundaryRange[0], boundaryRange[1]) &&
            std::abs(xy.xrec[ind]-xy.xrec[ind-1]) <= maxDx[0]) ||
           ((xy.xrec[ind] < std::min(boundaryRange[0], boundaryRange[1]) ||
            xy.xrec[ind] >= std::max(boundaryRange[0], boundaryRange[1])) &&
            std::abs(xy.xrec[ind]-xy.xrec[ind-1]) <= maxDx[1])) {
          xy.ydiff[ind-1] = xy.yrec[ind]-xy.yrec[ind-1];
          xy.xdiff[ind-1] = xy.xrec[ind]-xy.xrec[ind-1];
          xy.dydx[ind-1] = xy.ydiff[ind-1]/max_abs(xy.xdiff[ind-1], dxResolution);
        }  // only set ydiff, xdiff, dydx, if dx not too large
      }  // x and y not missing
    }  // ind > 0
  }  // for each obs in the record
}

// -----------------------------------------------------------------------------
/// Given the tolerance factors and x-boundaries,
///  set the vector of tolerance values for every obs in the record (group).

std::vector<float> SpikeAndStepCheck::set_tolerances(const std::vector<float> &xrec,
                                    const std::vector<size_t> &obs_indices,
                                    const Parameters_ &parameters_) const {
  const float tolerance = parameters_.toleranceOptions.value().nominal.value();
  std::vector<float> toleranceFactors;
  std::vector<float> toleranceBoundaries;
  size_t nSections = 0;
  if (parameters_.toleranceOptions.value().toleranceBoundaries.value() != boost::none) {
    toleranceFactors =
      parameters_.toleranceOptions.value().toleranceFactors.value().value();
    toleranceBoundaries =
      parameters_.toleranceOptions.value().toleranceBoundaries.value().value();
    nSections = toleranceBoundaries.size();
  }
  std::vector<float> toleranceSlope;
  if (nSections == 1) {
    toleranceSlope.push_back(0);
  } else if (nSections > 1) {
    for (size_t iSec = 0; iSec < nSections-1; ++iSec) {
      if (std::abs(toleranceBoundaries[iSec+1]-toleranceBoundaries[iSec]) <
          std::numeric_limits<float>::epsilon()) {
        toleranceSlope.push_back(0);  // can't compute, denominator = 0
      } else {
        // slope in section bounded by toleranceBoundaries[iSec], toleranceBoundaries[iSec+1]
        toleranceSlope.push_back((toleranceFactors[iSec+1]-toleranceFactors[iSec])/
                                 (toleranceBoundaries[iSec+1]-toleranceBoundaries[iSec]));
      }
    }
    toleranceSlope.push_back(0);  // zero slope beyond final section
  }
  // tolerance should be set correctly whether x ascending or descending - AS LONG AS
  //   tolerance boundaries and factors are given in order of ascending x:
  std::vector<float> tolerances;
  if (nSections > 0) {
    for (size_t ind = 0; ind < obs_indices.size(); ++ind) {
      size_t iSec = 0;
      while (iSec < nSections) {
        if (xrec[ind] <= toleranceBoundaries[iSec]) {
          if (iSec == 0) {
            tolerances.push_back(tolerance);  // below lowest boundary
            break;  // next obs index
          } else {
            tolerances.push_back(tolerance*(toleranceFactors[iSec-1] + toleranceSlope[iSec-1]*
                                                    (xrec[ind]-toleranceBoundaries[iSec-1])));
            break;  // next obs index
          }
        } else {  // xrec[ind] > toleranceBoundaries[iSec]
          if (iSec == nSections-1) {
            tolerances.push_back(tolerance*(toleranceFactors[iSec]));  // above final boundary
            break;  // next obs index
          }
          iSec++;
        }
      }  // while there are more tolerance sections to populate
    }  // for obs in record
  } else {  // no toleranceBoundaries given
    for (size_t obsIndex : obs_indices) {
      tolerances.push_back(tolerance);
    }
  }  // nSections > 0 or not
  return tolerances;
}

// -----------------------------------------------------------------------------
/// Go through the obs in the record (group), finding which are spikes and steps,
///  according to given conditions, flag as appropriate, and increment respective counters.

void SpikeAndStepCheck::identifyThinnedObservations(std::vector<bool> &isThinned,
                                          std::vector<bool> &spikeFlag,
                                          std::vector<bool> &stepFlag,
                                          xyStruct &xy,
                                          const std::vector<float> &tolerances,
                                          const std::vector<size_t> &obs_indices,
                                          const Parameters_ &parameters_) const {
  if (xy.ydiff.size() < 2) {
    return;  // don't do anything if record is <=2 obs long (can't tell spike or step)
  }
  const float missingValueFloat = util::missingValue<float>();
  const bool yesSpikes = parameters_.yesSpikes.value();
  const bool yesSteps = parameters_.yesSteps.value();
  const float gradientTolerance =
    parameters_.toleranceOptions.value().gradientTolerance.value();
  const float threshold = parameters_.toleranceOptions.value().threshold.value();
  const float smallThreshold = parameters_.toleranceOptions.value().smallThreshold.value();
  std::vector<float> boundaryRange;
  std::vector<float> stepTolRange;
  if (parameters_.boundaryOptions.value() == boost::none) {
    boundaryRange = {0, 0};
    stepTolRange = {0, 0};
  } else {
    boundaryRange = parameters_.boundaryOptions.value().value().boundaryRange.value();
    stepTolRange = parameters_.boundaryOptions.value().value().stepTolRange.value();
  }
  oops::Log::debug() << "Identifying spikes and steps:" << std::endl;
  size_t lastChecked = 0;
  for (size_t ind = 0; ind < obs_indices.size(); ++ind) {
    if (yesSpikes &&
       (ind < xy.ydiff.size()-1) &&
       (ind > lastChecked) &&
       (xy.ydiff[ind] != missingValueFloat && xy.ydiff[ind+1] != missingValueFloat)) {
      if (std::abs(xy.ydiff[ind]) > tolerances[ind+1] ||
          std::abs(xy.ydiff[ind+1]) > tolerances[ind+1]) {
        if (std::abs(xy.ydiff[ind]+xy.ydiff[ind+1]) < threshold*tolerances[ind+1]) {
          isThinned[obs_indices[ind+1]] = true;  // spike
          spikeFlag[obs_indices[ind+1]] = true;
          lastChecked = ind+1;
          oops::Log::debug() << "Large spike at " << obs_indices[ind+1] << std::endl;
        } else {  // if (dx[i]+dx[i+1]) can be a denominator:
          if (std::abs(xy.xdiff[ind] + xy.xdiff[ind+1]) > std::numeric_limits<float>::epsilon()) {
            // if y[i+1] is within 0.5*tolerances[i] of the value interpolated
            //  between y[i] and y[i+2], it is not considered a spike after all,
            //  but part of a large steady gradient. y[i+1](interp'd) =
            //  y[i] + (x[i+1]-x[i])*(y[i+2]-y[i])/(x[i+2]-x[i]) =
            //  y[i] + xdiff[i]*(ydiff[i+1]+ydiff[i])/(xdiff[i+1]+xdiff[i])
            //  i.e. if abs diff between this and y[i+1], is < 0.5*tolerances[i]
            //  i.e. abs(ydiff[i] - xdiff[i]*(ydiff[i+1]+ydiff[i])/(xdiff[i+1]+xdiff[i]))
            //          < 0.5*tolerances[i]
            //  LHS = abs( (xdiff[i+1]*ydiff[i] - xdiff[i]*ydiff[i+1])/
            //             (xdiff[i] + xdiff[i+1]) )
            const float yInterpDiff = (xy.xdiff[ind+1]*xy.ydiff[ind] -
                                xy.xdiff[ind]*xy.ydiff[ind+1])/
                                (xy.xdiff[ind] + xy.xdiff[ind+1]);
            if (std::abs(yInterpDiff) < threshold*tolerances[ind+1]) {
              isThinned[obs_indices[ind+1]] = false;  // not a spike after all
              spikeFlag[obs_indices[ind+1]] = false;
              lastChecked = ind+1;  // OPS bug - should be lastChecked = ind;
              oops::Log::debug() << "Not spike at " << obs_indices[ind+1] << std::endl;
            }  // not a spike after all
          }  // can (dx[i]+dx[i+1]) be a denominator or not
        }  // suspect large spike
      } else if ((std::abs(xy.ydiff[ind]) > 0.5*tolerances[ind+1] ||
                  std::abs(xy.ydiff[ind+1]) > 0.5*tolerances[ind+1]) &&
                 (std::abs(xy.dydx[ind]) > gradientTolerance ||
                  std::abs(xy.dydx[ind+1]) > gradientTolerance) &&
                 (std::abs(xy.ydiff[ind]+xy.ydiff[ind+1]) <
                  smallThreshold*std::abs(xy.ydiff[ind]-xy.ydiff[ind+1]))) {
        isThinned[obs_indices[ind+1]] = true;  // smaller spike
        spikeFlag[obs_indices[ind+1]] = true;
        lastChecked = ind+1;
        oops::Log::debug() << "Small spike at " << obs_indices[ind+1] << std::endl;
      }  // spike or not
    }  // count spikes
    if (yesSteps &&
       (ind < xy.ydiff.size()) &&
       (ind > lastChecked) &&
       (xy.ydiff[ind] != missingValueFloat)) {
      if (std::abs(xy.ydiff[ind]) > tolerances[ind]) {
        isThinned[obs_indices[ind+1]] = true;  // end of step
        stepFlag[obs_indices[ind+1]] = true;
        lastChecked = ind+1;
        if (ind < xy.ydiff.size()-1) {  // only flag start of step if not final level
          isThinned[obs_indices[ind]] = true;  // start of step
        }
        oops::Log::debug() << "Step at " << obs_indices[ind] << " and " <<
                              obs_indices[ind+1] << std::endl;
        // unless thermocline/halocline/other expected sharp change:
        if (xy.xrec[ind] >= std::min(boundaryRange[0], boundaryRange[1]) &&
            xy.xrec[ind] <= std::max(boundaryRange[0], boundaryRange[1]) &&
            xy.ydiff[ind] > tolerances[ind]*
                            std::min(stepTolRange[0], stepTolRange[1]) &&
            xy.ydiff[ind] < tolerances[ind]*
                            std::max(stepTolRange[0], stepTolRange[1])) {
          isThinned[obs_indices[ind+1]] = false;  // not step after all
          stepFlag[obs_indices[ind+1]] = false;
          lastChecked = ind+1;
          if (ind < xy.ydiff.size()-1) {  // only unflag start of step if not final level
            isThinned[obs_indices[ind]] = false;  // not step after all
          }
          oops::Log::debug() << "Not step after all at " << obs_indices[ind] <<
                                " and " << obs_indices[ind+1] << std::endl;
        }  // thermocline/halocline/etc., not suspect step
      }  // step or not
    }  // count steps
  }  // for each level in record
}

// -----------------------------------------------------------------------------

void SpikeAndStepCheck::print(std::ostream & os) const {
  os << "SpikeAndStepCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
