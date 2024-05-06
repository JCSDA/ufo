/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/BennartzScatIndex.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<BennartzScatIndex> makerBennartzScatIndex_("BennartzScatIndex");

BennartzScatIndex::BennartzScatIndex(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Get channels for computing scattering index from options
  channels_ = {options_.ch89.value(), options_.ch150.value()};

  // Include list of required data from ObsSpace
  invars_ += Variable("ObsValue/brightnessTemperature", channels_);
  if (options_.applyBias.value().size()) {
    invars_ += Variable(options_.applyBias.value()+"/brightnessTemperature", channels_);
  }
  invars_ += Variable("MetaData/sensorZenithAngle");
}

// -----------------------------------------------------------------------------

void BennartzScatIndex::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get dimension
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue<float>();

  ASSERT(out.nvars() == 1);

  // Get satellite zenith angle
  std::vector<float> satzen(nlocs);
  in.get(Variable("MetaData/sensorZenithAngle"), satzen);

  // Get observation values for required channels
  std::vector<float> bt89(nlocs), bt150(nlocs);
  in.get(Variable("ObsValue/brightnessTemperature", channels_)[0], bt89);
  in.get(Variable("ObsValue/brightnessTemperature", channels_)[1], bt150);

  // Apply bias correction if apply_bias is present in filter options
  std::vector<float> bias89(nlocs), bias150(nlocs);
  if (options_.applyBias.value().size()) {
    in.get(Variable(options_.applyBias.value()+"/brightnessTemperature", channels_)[0], bias89);
    in.get(Variable(options_.applyBias.value()+"/brightnessTemperature", channels_)[1], bias150);
    // Apply bias correction
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      bt89[iloc] -= bias89[iloc];
      bt150[iloc] -= bias150[iloc];
    }
  }

  // Compute Bennartz scattering index
  float Offset;
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    Offset = options_.coeff1.value() + options_.coeff2.value()*satzen[iloc];
    if (bt89[iloc] == missing || bt150[iloc] == missing || satzen[iloc] == missing) {
      out[0][iloc] = missing;
    } else {
      out[0][iloc] = bt89[iloc] - bt150[iloc] - Offset;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & BennartzScatIndex::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
