/* * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ProfileBackgroundCheck.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// ProfileBackgroundCheck: check observation closeness to background over a profile
///
/// This will use the record number in obsdb_ to identify which observations belong to a
/// given profile (all members of a profile must share the same record number).  For
/// each profile the RMS of observed-HofX is calculated, and this is compared against
/// the given threshold.  If "relative threshold" is used (rather than "absolute threshold")
/// then the RMS is normalised by observation error.
///
/// The threshold can be a floating-point constant or the name of a variable indicating
/// the threshold to use at each location.  In this case, the value of the threshold is
/// taken from the threshold given for the first observation in the profile.
///
/// This is related to BackgroundCheck, which checks each observation against a threshold.
/// There is also a group of other profile checks (ConventionalProfileProcessing) which are
/// mostly aimed to processing radiosondes.

ProfileBackgroundCheck::ProfileBackgroundCheck(
        ioda::ObsSpace & obsdb,
        const Parameters_ & parameters,
        std::shared_ptr<ioda::ObsDataVector<int> > flags,
        std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "ProfileBackgroundCheck constructor" << std::endl;

  ASSERT(parameters_.relativeThreshold.value() ||
         parameters_.absoluteThreshold.value());
}

// -----------------------------------------------------------------------------

ProfileBackgroundCheck::~ProfileBackgroundCheck() {
  oops::Log::trace() << "ProfileBackgroundCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/// Apply the profile background check filter.

void ProfileBackgroundCheck::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ProfileBackgroundCheck postFilter" << std::endl;
  const oops::ObsVariables observed = obsdb_.obsvariables();

  oops::Log::debug() << "ProfileBackgroundCheck obserr: " << *obserr_;

  ioda::ObsDataVector<float> obs(obsdb_, filtervars.toOopsObsVariables(), "ObsValue");
  ioda::ObsDataVector<float> bias(obsdb_, filtervars.toOopsObsVariables(), "ObsBias", false);

  // Get the record numbers from the observation data.  These will be used to identify
  // which observations belong to which profile.
  const std::vector<size_t> &unique = obsdb_.recidx_all_recnums();
  oops::Log::debug() << "Unique record numbers" << std::endl;
  for (size_t iUnique : unique)
    oops::Log::debug() << iUnique << ' ';
  oops::Log::debug() << std::endl;

  // Threshold for all variables
  std::vector<float> abs_thr(obsdb_.nlocs(), std::numeric_limits<float>::max());
  std::vector<float> rel_thr(obsdb_.nlocs(), std::numeric_limits<float>::max());
  if (parameters_.absoluteThreshold.value())
    abs_thr = getScalarOrFilterData(*parameters_.absoluteThreshold.value(), data_);
  if (parameters_.relativeThreshold.value())
    rel_thr = getScalarOrFilterData(*parameters_.relativeThreshold.value(), data_);

  Variables varhofx(filtervars_, "HofX");
  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    size_t iv = observed.find(filtervars.variable(jv).variable());
    // H(x)
    std::vector<float> hofx;
    data_.get(varhofx.variable(jv), hofx);

    // Loop over the unique profiles
    for (const auto& iprofile : unique) {
      std::vector<size_t> obs_numbers = obsdb_.recidx_vector(iprofile);
      // Initialise the accumulators for this profile
      double total_diff = 0;
      int total_nobs = 0;
      bool using_rel_thresh = false;
      bool first_ob = true;
      float profile_threshold = 0;

      // First run through all the observations to calculate the O-B statistics
      // (normalised by the observation error).
      // Determine whether this profile uses threshold or abs_threshold based
      // on the first value in the profile.
      for (size_t jobs : obs_numbers) {
        if (apply[jobs] && (*flags_)[iv][jobs] == QCflags::pass) {
          // style note: use missingValue<decltype(foo)> because obserr_ is declared in another file
          // and its underlying type could in principle change without it being obvious here.
          ASSERT((*obserr_)[iv][jobs] !=
              util::missingValue<std::remove_reference_t<decltype((*obserr_)[iv][jobs])>>());
          ASSERT(obs[jv][jobs] != util::missingValue<float>());
          ASSERT(hofx[jobs] != util::missingValue<float>());

          if (first_ob) {
            if (rel_thr[jobs] == std::numeric_limits<float>::max()) {
              using_rel_thresh = false;
              profile_threshold = abs_thr[jobs];
            } else {
              using_rel_thresh = true;
              profile_threshold = rel_thr[jobs];
            }
            first_ob = false;
            ASSERT(profile_threshold < std::numeric_limits<float>::max() &&
                   profile_threshold > 0.0f);
          }

          // Apply bias correction
          const float yy = obs[jv][jobs] - bias[jv][jobs];

          // Accumulate distance from background
          if ((*obserr_)[iv][jobs] > 0) {
            total_nobs++;
            if (using_rel_thresh) {
              // Normalise by observation error
              total_diff += std::pow((hofx[jobs] - yy) / (*obserr_)[iv][jobs], 2);
            } else {
              // Use the absolute differences
              total_diff += std::pow(hofx[jobs] - yy, 2);
            }
          }
        }
      }

      // Calculate average difference over profile
      if (total_nobs > 0) {
        total_diff = std::pow(total_diff / total_nobs, 0.5);
      }
      oops::Log::debug() << "ProfileBackgroundCheck: For profile " << iprofile <<
                            " rms diff is " << total_diff << std::endl;

      // Second run through the observations to apply the threshold to reject profiles
      for (size_t jobs : obs_numbers) {
        if (apply[jobs] && (*flags_)[iv][jobs] == QCflags::pass) {
          // If rejecting this profile, then flag all elements
          if (total_diff > static_cast<double>(profile_threshold)) flagged[jv][jobs] = true;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

void ProfileBackgroundCheck::print(std::ostream & os) const {
  os << "ProfileBackgroundCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
