/* * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"
#include "ufo/filters/FilterUtils.h"
#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/filters/ObsAccessor.h"
#include "ufo/filters/OceanVerticalStabilityCheck.h"
#include "ufo/filters/SpikeAndStepCheck.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// OceanVerticalStabilityCheck: Flag observations that are likely erroneous because they
///   constitute a density inversion.
OceanVerticalStabilityCheck::OceanVerticalStabilityCheck(
        ioda::ObsSpace & obsdb,
        const Parameters_ & parameters,
        std::shared_ptr<ioda::ObsDataVector<int> > flags,
        std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "OceanVerticalStabilityCheck constructor" << std::endl;
}

// -----------------------------------------------------------------------------

OceanVerticalStabilityCheck::~OceanVerticalStabilityCheck() {
  oops::Log::trace() << "OceanVerticalStabilityCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/// Apply the Ocean Vertical Stability check filter.

void OceanVerticalStabilityCheck::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "OceanVerticalStabilityCheck filter" << std::endl;

  ObsAccessor obsAccessor =
    ObsAccessor::toObservationsSplitIntoIndependentGroupsByRecordId(obsdb_);
  const size_t totalNumObs = obsAccessor.totalNumObservations();

  const std::vector<float> salinity = obsAccessor.getFloatVariableFromObsSpace
                                    (parameters_.varParams.value().salinity.value().group(),
                                     parameters_.varParams.value().salinity.value().variable());
  const std::vector<float> temperature = obsAccessor.getFloatVariableFromObsSpace
                                    (parameters_.varParams.value().temperature.value().group(),
                                     parameters_.varParams.value().temperature.value().variable());
  const std::vector<float> pressure = obsAccessor.getFloatVariableFromObsSpace
                                    (parameters_.varParams.value().pressure.value().group(),
                                     parameters_.varParams.value().pressure.value().variable());

  std::vector<bool> isThinned(totalNumObs, false);
  std::vector<bool> spikeFlag(totalNumObs, false);
  std::vector<bool> stepFlag(totalNumObs, false);
  // get diagnostic flags from ObsSpace (can't create them in the code, must do in YAML...
  //  ...must be there when ObsSpace constructed)
  for (size_t filterVarIndex = 0; filterVarIndex < filtervars.size(); ++filterVarIndex) {
    const std::string filterVarName = filtervars.variable(filterVarIndex).variable();
    if (obsdb_.has("DiagnosticFlags/DensitySpike", filterVarName)) {
      obsdb_.get_db("DiagnosticFlags/DensitySpike", filterVarName, spikeFlag);
    } else {
      throw eckit::UserError("Variable 'DiagnosticFlags/DensitySpike/" + filterVarName +
                            "' does not exist yet. "
                            "It needs to be set up with the 'Create Diagnostic Flags' filter "
                            "prior to using the 'set' or 'unset' action.");
    }
    if (obsdb_.has("DiagnosticFlags/DensityStep", filterVarName)) {
      obsdb_.get_db("DiagnosticFlags/DensityStep", filterVarName, stepFlag);
    } else {
      throw eckit::UserError("Variable 'DiagnosticFlags/DensityStep/" + filterVarName +
                            "' does not exist yet. "
                            "It needs to be set up with the 'Create Diagnostic Flags' filter "
                            "prior to using the 'set' or 'unset' action.");
    }
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

    // set level-to-level potential density difference for this profile:
    const std::vector<float> density_diff = calculateDensityDiff(salinity,
                                                                 temperature,
                                                                 pressure,
                                                                 obs_indices);

    // set flags for obs that fulfil density spike/step conditions:
    identifyDensityInversions(isThinned,
                              spikeFlag,
                              stepFlag,
                              density_diff,
                              obs_indices);
  }  // for each record
  obsAccessor.flagRejectedObservations(isThinned, flagged);
  for (size_t filterVarIndex = 0; filterVarIndex < filtervars.size(); ++filterVarIndex) {
    const std::string filterVarName = filtervars.variable(filterVarIndex).variable();
    obsdb_.put_db("DiagnosticFlags/DensitySpike", filterVarName, spikeFlag);
    obsdb_.put_db("DiagnosticFlags/DensityStep", filterVarName, stepFlag);
  }
}

// -----------------------------------------------------------------------------
/// Given salinity (g/kg), temperature (deg.C) and pressure (dbar),
///  set density difference, drho, for the given record (group).

std::vector<float> OceanVerticalStabilityCheck::calculateDensityDiff
                                    (const std::vector<float> &salinity,
                                    const std::vector<float> &temperature,
                                    const std::vector<float> &pressure,
                                    const std::vector<size_t> &obs_indices) const {
  if (obs_indices.size() < 3) {  // don't do anything if record is <=2 obs long
    std::vector<float> densityDiff;
    return densityDiff;
  }
  const float missingValueFloat = util::missingValue<float>();
  std::vector<float> densityDiff(obs_indices.size()-1, missingValueFloat);
  for (size_t ind = 1; ind < obs_indices.size(); ++ind) {
    if (salinity[obs_indices[ind]] != missingValueFloat &&
        salinity[obs_indices[ind-1]] != missingValueFloat &&
        temperature[obs_indices[ind]] != missingValueFloat &&
        temperature[obs_indices[ind-1]] != missingValueFloat &&
        pressure[obs_indices[ind]] != missingValueFloat) {
      const float density = gsw_rho_t_exact_f90(salinity[obs_indices[ind]],
                                    temperature[obs_indices[ind]],
                                    pressure[obs_indices[ind]]);
      const float densityPrevious = gsw_rho_t_exact_f90(salinity[obs_indices[ind-1]],
                                     temperature[obs_indices[ind-1]],
                                     pressure[obs_indices[ind]]);  // intentional use of [ind]
                                                                   //  rather than [ind-1]
      densityDiff[ind-1] = density - densityPrevious;
    }
  }  // for each obs in the record
  return densityDiff;
}

// -----------------------------------------------------------------------------
/// Go through the obs in the record (group), finding which are spikes and steps,
///  according to given conditions, flag as appropriate.

void OceanVerticalStabilityCheck::identifyDensityInversions(std::vector<bool> &isThinned,
                                          std::vector<bool> &spikeFlag,
                                          std::vector<bool> &stepFlag,
                                          const std::vector<float> &densityDiff,
                                          const std::vector<size_t> &obs_indices) const {
  if (densityDiff.size() < 2) {
    return;  // don't do anything if record is <=2 obs long (can't tell spike or step)
  }
  const float missingValueFloat = util::missingValue<float>();
  const bool yesSpikes = parameters_.yesSpikes.value();
  const bool yesSteps = parameters_.yesSteps.value();
  const float tolerance = parameters_.nominal.value();
  const float threshold = parameters_.threshold.value();
  oops::Log::debug() << "Identifying spikes and steps:" << std::endl;
  size_t lastChecked = 0;
  for (size_t ind = 0; ind < obs_indices.size(); ++ind) {
    if (yesSpikes &&
       (ind < densityDiff.size()-1) &&
       (ind > lastChecked) &&
       (densityDiff[ind] != missingValueFloat &&
        densityDiff[ind+1] != missingValueFloat)) {
      if ( (densityDiff[ind] < tolerance || densityDiff[ind+1] < tolerance) &&
           (std::abs(densityDiff[ind] + densityDiff[ind+1]) <
            threshold*std::abs(densityDiff[ind] - densityDiff[ind+1])) ) {
          isThinned[obs_indices[ind+1]] = true;  // spike
          spikeFlag[obs_indices[ind+1]] = true;
          lastChecked = ind+1;
          oops::Log::debug() << "Density spike at " << obs_indices[ind+1] << std::endl;
      }  // spike or not
    }  // count spikes
    if (yesSteps &&
       (ind < densityDiff.size()) &&
       (!spikeFlag[obs_indices[ind+1]]) &&
       (densityDiff[ind] != missingValueFloat)) {
      if (densityDiff[ind] < tolerance) {
        // only flag start of step if NOT final level:
        if (ind < densityDiff.size()-1) {
          isThinned[obs_indices[ind]] = true;  // start of step
          oops::Log::debug() << "Density step at " << obs_indices[ind] << " and " <<
                                obs_indices[ind+1] << std::endl;
        } else {
          oops::Log::debug() << "Density step at " << obs_indices[ind+1] <<
                                " (final level)." << std::endl;
        }
        isThinned[obs_indices[ind+1]] = true;  // end of step
        stepFlag[obs_indices[ind+1]] = true;
      }  // step or not
    }  // count steps
  }  // for each level in record
}

// -----------------------------------------------------------------------------

void OceanVerticalStabilityCheck::print(std::ostream & os) const {
  os << "OceanVerticalStabilityCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
