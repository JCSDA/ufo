/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <memory>

#include "ioda/distribution/Accumulator.h"
#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
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
  invars_ += Variable("ObsValue/brightnessTemperature", channels_);
  invars_ += Variable(hofxopt+"/brightnessTemperature", channels_);

  if (options_.ObsBias.value().size()) {
    // Get optional bias correction
    const std::string & biasopt = options_.ObsBias.value();
    invars_ += Variable(biasopt+"/brightnessTemperature", channels_);
  }
}

// -----------------------------------------------------------------------------

void BgdDepartureAnomaly::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get dimension
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue<float>();

  // Get obs space
  auto & obsdb = in.obsspace();

  ASSERT(channels_.size() == 2);

  // Get observation values for required channels
  const std::string & hofxopt = options_.testHofX.value();
  std::vector<float> obslow(nlocs), obshigh(nlocs);
  std::vector<float> bgdlow(nlocs), bgdhigh(nlocs);
  in.get(Variable("ObsValue/brightnessTemperature", channels_)[0], obslow);
  in.get(Variable("ObsValue/brightnessTemperature", channels_)[1], obshigh);
  in.get(Variable(hofxopt+"/brightnessTemperature", channels_)[0], bgdlow);
  in.get(Variable(hofxopt+"/brightnessTemperature", channels_)[1], bgdhigh);

  // Get bias correction if ObsBias is present in filter options
  if (options_.ObsBias.value().size()) {
    const std::string & biasopt = options_.ObsBias.value();
    std::vector<float> biaslow(nlocs), biashigh(nlocs);
    in.get(Variable(biasopt+"/brightnessTemperature", channels_)[0], biaslow);
    in.get(Variable(biasopt+"/brightnessTemperature", channels_)[1], biashigh);

    // Apply bias correction
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      bgdlow[iloc] -= biaslow[iloc];
      bgdhigh[iloc] -= biashigh[iloc];
    }
  }

  // Compute background departures

  // To calculate mean departures, we need to know the number of observations with non-missing
  // departures (taking into account all MPI ranks)...
  std::unique_ptr<ioda::Accumulator<size_t>> countAccumulator =
      obsdb.distribution()->createAccumulator<size_t>();

  // ... and the sum of all non-missing departures
  enum {OMB_LOW, OMB_HIGH, NUM_OMBS};
  std::unique_ptr<ioda::Accumulator<std::vector<double>>> totalsAccumulator =
      obsdb.distribution()->createAccumulator<double>(NUM_OMBS);

  // First, perform local reductions...
  std::vector<float> omblow(nlocs), ombhigh(nlocs);
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (obslow[iloc] == missing || bgdlow[iloc] == missing ||
        obshigh[iloc] == missing || bgdhigh[iloc] == missing) {
      omblow[iloc] = missing;
      ombhigh[iloc] = missing;
    } else {
      omblow[iloc] = obslow[iloc] - bgdlow[iloc];
      ombhigh[iloc] = obshigh[iloc] - bgdhigh[iloc];
      totalsAccumulator->addTerm(iloc, OMB_LOW, omblow[iloc]);
      totalsAccumulator->addTerm(iloc, OMB_HIGH, ombhigh[iloc]);
      countAccumulator->addTerm(iloc, 1);
    }
  }
  // ... and then global reductions.
  const std::size_t count = countAccumulator->computeResult();
  const std::vector<double> totals = totalsAccumulator->computeResult();

  // Calculate the means
  const double MeanOmbLow = totals[OMB_LOW] / count;
  const double MeanOmbHigh = totals[OMB_HIGH] / count;

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
