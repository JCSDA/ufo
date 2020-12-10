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
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<ImpactHeight> makerImpactHeight_("ImpactHeight");

/* -----------------------------------------------------------------------------
 * Specify that impact_parameter and earth_radius_of_curvature need to be
 * provided to this function
 * -----------------------------------------------------------------------------
 */
ImpactHeight::ImpactHeight(const eckit::LocalConfiguration & conf)
  : invars_() {
  invars_ += Variable("impact_parameter@MetaData");
  invars_ += Variable("earth_radius_of_curvature@MetaData");
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
  std::vector<float> impact_parameter;
  in.get(Variable("impact_parameter@MetaData"), impact_parameter);
  std::vector<float> radius_curvature;
  in.get(Variable("earth_radius_of_curvature@MetaData"), radius_curvature);
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (impact_parameter[jj] == util::missingValue(impact_parameter[jj]) ||
        radius_curvature[jj] == util::missingValue(radius_curvature[jj])) {
        out[0][jj] = util::missingValue(out[0][jj]);
    } else {
        out[0][jj] = impact_parameter[jj] - radius_curvature[jj];
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ImpactHeight::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
