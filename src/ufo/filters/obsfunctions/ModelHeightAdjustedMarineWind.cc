/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ModelHeightAdjustedMarineWind.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<ModelHeightAdjustedEastwardMarineWind>
        eastwardMaker_("ModelHeightAdjustedEastwardMarineWind");
static ObsFunctionMaker<ModelHeightAdjustedNorthwardMarineWind>
        northwardMaker_("ModelHeightAdjustedNorthwardMarineWind");

// -----------------------------------------------------------------------------

ModelHeightAdjustedMarineWindComponent::ModelHeightAdjustedMarineWindComponent(
        const eckit::LocalConfiguration & conf, const Variable &windComponent)
        : invars_(), wind_(windComponent) {
  // Required observation station height
  invars_ += Variable("MetaData/anemometerHeight");
  // Required wind component
  invars_ += wind_;
}

// -----------------------------------------------------------------------------

void ModelHeightAdjustedMarineWindComponent::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  std::vector<float> WindComponent(nlocs);
  std::vector<float> StationHeight(nlocs);

  in.get(Variable(wind_), WindComponent);
  in.get(Variable("MetaData/anemometerHeight"), StationHeight);

  const float missing = util::missingValue<float>();
  const float a = 1.0/0.0016;
  const float ref_height = std::log(10.0*a);
  // TODO(jwaller): When QC flaging is availiable add flags indicate corrected values
  // Compute adjusted marine wind.
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (StationHeight[jj] == missing || WindComponent[jj] == missing) {
      out[0][jj] = missing;
    } else {
      float ScaleFactor = ref_height/std::log(std::abs(StationHeight[jj])*a);
      out[0][jj] = WindComponent[jj]*ScaleFactor;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ModelHeightAdjustedMarineWindComponent::requiredVariables() const {
  return invars_;
}


// -----------------------------------------------------------------------------

}  // namespace ufo
