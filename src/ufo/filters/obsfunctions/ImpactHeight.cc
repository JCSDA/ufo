/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------
 * Function to calculate the impact height for GNSS-RO.  This is the difference
 * between the impact parameter and the earth's radius of curvature.
 * -----------------------------------------------------------------------------
 */
#include "ufo/filters/obsfunctions/ImpactHeight.h"

#include <math.h>
#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

#include "eckit/exception/Exceptions.h"

namespace ufo {

static ObsFunctionMaker<ImpactHeight> makerImpactHeight_("ImpactHeight");

/* -----------------------------------------------------------------------------
 * Specify that impactParameterRO and earthRadiusCurvature need to be
 * provided to this function
 * -----------------------------------------------------------------------------
 */
ImpactHeight::ImpactHeight(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Get channels from options
  options_.deserialize(conf);
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));

  invars_ += Variable("MetaData/impactParameterRO", channels_);
  invars_ += Variable("MetaData/earthRadiusCurvature");
}

// -----------------------------------------------------------------------------

ImpactHeight::~ImpactHeight() {}

/* -----------------------------------------------------------------------------
 * Perform the computation.  Read in the required variables, and calculate
 * their difference, storing the difference in the output vector.
 * -----------------------------------------------------------------------------
 */
void ImpactHeight::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const size_t nchans = out.nvars();

  if (nchans == 1 && channels_.size() == 0) {
    // All good, do nothing
  } else if (nchans != channels_.size()) {
    std::string errString = "ImpactHeight: mismatch between channels (" +
      std::to_string(channels_.size()) + ") and number of variables (" +
      std::to_string(nchans) + ")";
    throw eckit::BadValue(errString);
  }

  std::vector<float> radius_curvature;
  in.get(Variable("MetaData/earthRadiusCurvature"), radius_curvature);
  for (size_t ichan=0; ichan < nchans ; ichan++) {
    std::vector<float> impact_parameter;
    in.get(Variable("MetaData/impactParameterRO", channels_)[ichan], impact_parameter);
    for (size_t jj = 0; jj < nlocs; ++jj) {
      if (impact_parameter[jj] == util::missingValue<float>() ||
          radius_curvature[jj] == util::missingValue<float>()) {
          out[ichan][jj] = util::missingValue<float>();
      } else {
          out[ichan][jj] = impact_parameter[jj] - radius_curvature[jj];
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ImpactHeight::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
