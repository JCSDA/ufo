/*
 * (C) Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <set>

#include "ufo/filters/obsfunctions/CloudFirstGuessMinimumResidual.h"

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace ufo {

static ObsFunctionMaker<CloudFirstGuessMinimumResidual>
       makerCloudFirstGuessMinimumResidual_("CloudFirstGuessMinimumResidual");

// -----------------------------------------------------------------------------

CloudFirstGuessMinimumResidual::CloudFirstGuessMinimumResidual(
        const eckit::LocalConfiguration & conf)
  : invars_(), channels_() {
  oops::Log::trace() << "CloudFirstGuessMinimumResidual constructor start" << std::endl;

  // Initialize options
  options_.validateAndDeserialize(conf);

  // Get channels used
  std::set<int> chanset = oops::parseIntSet(options_.channels.value());
  channels_.assign(chanset.begin(), chanset.end());
  ASSERT(channels_.size() > 0);

  // List of required data
  invars_ += Variable("ObsValue/brightness_temperature", channels_);
  invars_ += Variable("ObsError/brightness_temperature", channels_);
  invars_ += Variable(options_.obsBiasGroup.value() + "/brightness_temperature", channels_);
  invars_ += Variable("ObsDiag/brightness_temperature_assuming_clear_sky", channels_);
  invars_ += Variable("ObsDiag/brightness_temperature_from_atmosphere_layer_to_toa", channels_);
  invars_ += Variable("GeoVaLs/air_pressure");

  oops::Log::trace() << "CloudFirstGuessMinimumResidual constructor end" << std::endl;
}

// -----------------------------------------------------------------------------

CloudFirstGuessMinimumResidual::~CloudFirstGuessMinimumResidual() {}

// -----------------------------------------------------------------------------

void CloudFirstGuessMinimumResidual::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "CloudFirstGuessMinimumResidual compute start" << std::endl;
  ASSERT(out.nvars() == 1);

  /// Setup constant variables used throughout the routine
  const std::string btOvercastName = "ObsDiag/brightness_temperature_from_atmosphere_layer_to_toa";
  const std::vector<std::string> vars = {"brightness_temperature"};
  const std::vector<std::string> clearSkyName = {"brightness_temperature_assuming_clear_sky"};
  const float missing = util::missingValue(missing);
  const size_t nlocs = in.nlocs();
  const size_t nlevs = in.nlevs(Variable(btOvercastName, channels_)[0]);
  const size_t nchans = channels_.size();
  if (nlocs == 0) return;

  /// Create variables needed throughtout
  std::vector<std::vector<float>> cloudFraction(nlevs, std::vector<float>(nlocs, missing));
  std::vector<std::vector<float>> airPressure(nlevs, std::vector<float>(nlocs));
  std::vector<std::vector<float>> costFunction(nlocs, std::vector<float>(nlevs, 0.0f));
  ioda::ObsDataVector<float> obsVal(in.obsspace(), oops::Variables(vars, channels_));
  ioda::ObsDataVector<float> obsError(in.obsspace(), oops::Variables(vars, channels_));
  ioda::ObsDataVector<float> obsBias(in.obsspace(), oops::Variables(vars, channels_));
  ioda::ObsDataVector<float> obsClearVal(in.obsspace(), oops::Variables(clearSkyName, channels_));

  /// Read data from filterdata - for RTTOV ObsDiags are not bias corrected and therefore the
  /// ObsBias needs to be read in.
  const Variable clearVarName("ObsDiag/" + clearSkyName[0], channels_);
  const Variable obsBiasName(options_.obsBiasGroup.value() + "/" + vars[0], channels_);
  const Variable obsValName("ObsValue/" + vars[0], channels_);
  const Variable obsErrorName("ObsError/" + vars[0], channels_);
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(obsValName[ichan], obsVal);
    in.get(obsErrorName[ichan], obsError);
    in.get(clearVarName[ichan], obsClearVal);
    in.get(obsBiasName[ichan], obsBias);
  }

  /// Some checks before the loop
  ASSERT(obsClearVal.nlocs() == nlocs);
  ASSERT(obsVal.nlocs() == nlocs);
  ASSERT(obsVal.nvars() == nchans);

  /// Loop over the number of levels
  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
    in.get(Variable("GeoVaLs/air_pressure"), ilev, airPressure[ilev]);
    std::vector<float> cloudFractionNumerator(nlocs, 0.0f);
    std::vector<float> cloudFractionDenominator(nlocs, 0.0f);
    std::vector<std::vector<float>> obsCloudyVal(nchans, std::vector<float>(nlocs));

    /// This finds the optimum cloud fraction at each level.
    /// Calculates Equation 4 of Eyre and Menzel (1989) with error correction
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      in.get(Variable(btOvercastName, channels_)[ichan], ilev, obsCloudyVal[ichan]);

      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        float errorWeightedCloudyMinusClear =
                   (obsCloudyVal[ichan][iloc] - obsClearVal[ichan][iloc]) /
                      (obsError[ichan][iloc] * obsError[ichan][iloc]);
        float cloudFractionNumeratorIncrement =
                errorWeightedCloudyMinusClear *
                (obsVal[ichan][iloc] - obsBias[ichan][iloc] - obsClearVal[ichan][iloc]);
        float cloudFractionDenominatorIncrement =
                errorWeightedCloudyMinusClear *
                (obsCloudyVal[ichan][iloc] - obsClearVal[ichan][iloc]);
        cloudFractionNumerator[iloc] = cloudFractionNumerator[iloc] +
                                       cloudFractionNumeratorIncrement;
        cloudFractionDenominator[iloc] = cloudFractionDenominator[iloc] +
                                         cloudFractionDenominatorIncrement;
      }
    }

    /// This finds the cost function (residual) at each level
    /// Equation 5 of Eyra and Menzel (1989)
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (airPressure[ilev][iloc] < options_.minCloudPressure.value()) {
        costFunction[iloc][ilev] = 1.0e9;
        continue;
      }
      cloudFraction[ilev][iloc] = cloudFractionNumerator[iloc] / cloudFractionDenominator[iloc];
      cloudFraction[ilev][iloc] = std::min(cloudFraction[ilev][iloc], 1.0f);
      cloudFraction[ilev][iloc] = std::max(cloudFraction[ilev][iloc], 0.0f);
      for (size_t ichan = 0; ichan < nchans; ++ichan) {
        float LHS = obsVal[ichan][iloc] - obsBias[ichan][iloc] - obsClearVal[ichan][iloc];
        LHS *= LHS;
        float RHS = cloudFraction[ilev][iloc] *
                    (obsCloudyVal[ichan][iloc] - obsClearVal[ichan][iloc]);
        RHS *= RHS;
        float costFunctionIncrement = (LHS - RHS) /
                                      (obsError[ichan][iloc] * obsError[ichan][iloc]);
        costFunction[iloc][ilev] = costFunction[iloc][ilev] + costFunctionIncrement;
      }
    }
  }

  /// Use the minimum cost in a profile to evaluate the first guess cloud top pressure
  /// and effective cloud amount and save these to the ObsSpace.
  const std::vector<std::string> firstGuessNames = {"initial_cloud_top_pressure",
                                                    "initial_cloud_fraction"};
  ioda::ObsDataVector<float> firstGuessValues(in.obsspace(), oops::Variables(firstGuessNames));
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    size_t MinLevelIndex = std::min_element(costFunction[iloc].begin(), costFunction[iloc].end()) -
                           costFunction[iloc].begin();
    out[0][iloc] = costFunction[iloc][MinLevelIndex];
    firstGuessValues[firstGuessNames[0]][iloc] = airPressure[MinLevelIndex][iloc];
    firstGuessValues[firstGuessNames[1]][iloc] = cloudFraction[MinLevelIndex][iloc];
  }
  firstGuessValues.save(options_.outputGroup.value());
  oops::Log::trace() << "CloudFirstGuessMinimumResidual compute end" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
