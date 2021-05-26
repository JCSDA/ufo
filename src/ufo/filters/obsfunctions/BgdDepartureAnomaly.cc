/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/obsfunctions/BgdDepartureAnomaly.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<BgdDepartureAnomaly> makerBgdDepartureAnomaly_("BgdDepartureAnomaly");

BgdDepartureAnomaly::BgdDepartureAnomaly(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Get channels for computing the background departure
  channels_ = {options_.obslow.value(), options_.obshigh.value()};

  // Use default or test version of HofX
  const std::string & hofxopt = options_.testHofX.value();

  // Include list of required data from ObsSpace
  invars_ += Variable("brightness_temperature@ObsValue", channels_);
  invars_ += Variable("brightness_temperature@"+hofxopt, channels_);

  if (options_.ObsBias.value().size()) {
    // Get optional bias correction
    const std::string & biasopt = options_.ObsBias.value();
    invars_ += Variable("brightness_temperature@"+biasopt, channels_);
  }
}

// -----------------------------------------------------------------------------

void BgdDepartureAnomaly::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get dimension
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue(missing);

  // Get obs space
  auto & obsdb = in.obsspace();

  ASSERT(channels_.size() == 2);

  // Get observation values for required channels
  const std::string & hofxopt = options_.testHofX.value();
  std::vector<float> obslow(nlocs), obshigh(nlocs);
  std::vector<float> bgdlow(nlocs), bgdhigh(nlocs);
  in.get(Variable("brightness_temperature@ObsValue", channels_)[0], obslow);
  in.get(Variable("brightness_temperature@ObsValue", channels_)[1], obshigh);
  in.get(Variable("brightness_temperature@"+hofxopt, channels_)[0], bgdlow);
  in.get(Variable("brightness_temperature@"+hofxopt, channels_)[1], bgdhigh);

  // Get bias correction if ObsBias is present in filter options
  if (options_.ObsBias.value().size()) {
    const std::string & biasopt = options_.ObsBias.value();
    std::vector<float> biaslow(nlocs), biashigh(nlocs);
    in.get(Variable("brightness_temperature@"+biasopt, channels_)[0], biaslow);
    in.get(Variable("brightness_temperature@"+biasopt, channels_)[1], biashigh);

    // Apply bias correction
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      bgdlow[iloc] -= biaslow[iloc];
      bgdhigh[iloc] -= biashigh[iloc];
    }
  }

  // Compute background departures
  std::vector<float> omblow(nlocs), ombhigh(nlocs);
  double MeanOmbLow = 0, MeanOmbHigh = 0;
  int missingVal = 0;
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (obslow[iloc] == missing || bgdlow[iloc] == missing ||
        obshigh[iloc] == missing || bgdhigh[iloc] == missing) {
      omblow[iloc] = missing;
      ombhigh[iloc] = missing;
      missingVal += 1;
    } else {
      omblow[iloc] = obslow[iloc] - bgdlow[iloc];
      ombhigh[iloc] = obshigh[iloc] - bgdhigh[iloc];
      MeanOmbLow += omblow[iloc];
      MeanOmbHigh += ombhigh[iloc];
    }
  }

  size_t nobslow = obsdb.distribution()->globalNumNonMissingObs(obslow);
  size_t nobshigh = obsdb.distribution()->globalNumNonMissingObs(obshigh);
  obsdb.distribution()->allReduceInPlace(MeanOmbLow, eckit::mpi::sum());
  obsdb.distribution()->allReduceInPlace(MeanOmbHigh, eckit::mpi::sum());

  MeanOmbLow = MeanOmbLow / nobslow;
  MeanOmbHigh = MeanOmbHigh / nobshigh;

  // Compute anomaly difference
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (obslow[iloc] == missing || bgdlow[iloc] == missing ||
        obshigh[iloc] == missing || bgdhigh[iloc] == missing) {
      out[0][iloc] = missing;
    } else {
      out[0][iloc] = std::abs((omblow[iloc]-MeanOmbLow) - (ombhigh[iloc]-MeanOmbHigh));
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & BgdDepartureAnomaly::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
