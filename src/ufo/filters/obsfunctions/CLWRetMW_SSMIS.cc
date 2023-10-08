/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/CLWRetMW_SSMIS.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<CLWRetMW_SSMIS> makerCLWRetMW_SSMIS_("CLWRetMW_SSMIS");

CLWRetMW_SSMIS::CLWRetMW_SSMIS(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  channels_ = {options_.ch19h.value(), options_.ch19v.value(),
               options_.ch22v.value(), options_.ch37h.value(),
               options_.ch37v.value(), options_.ch91v.value(),
               options_.ch91h.value()};

  // Include list of required data from ObsSpace
  invars_ += Variable(options_.varGroup.value() + "/brightnessTemperature", channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/water_area_fraction");
}

// -----------------------------------------------------------------------------

CLWRetMW_SSMIS::~CLWRetMW_SSMIS() {}

// -----------------------------------------------------------------------------

void CLWRetMW_SSMIS::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get required parameters
  const std::string &vargrp = options_.varGroup.value();

  // Get dimension
  const size_t nlocs = in.nlocs();

  // Get area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  in.get(Variable("GeoVaLs/water_area_fraction"), water_frac);

  std::vector<float> bt19h(nlocs), bt19v(nlocs), bt22v(nlocs), bt37h(nlocs),
                     bt37v(nlocs), bt91v(nlocs), bt91h(nlocs);

  in.get(Variable(vargrp + "/brightnessTemperature", channels_)[0], bt19h);
  in.get(Variable(vargrp + "/brightnessTemperature", channels_)[1], bt19v);
  in.get(Variable(vargrp + "/brightnessTemperature", channels_)[2], bt22v);
  in.get(Variable(vargrp + "/brightnessTemperature", channels_)[3], bt37h);
  in.get(Variable(vargrp + "/brightnessTemperature", channels_)[4], bt37v);
  in.get(Variable(vargrp + "/brightnessTemperature", channels_)[5], bt91v);
  in.get(Variable(vargrp + "/brightnessTemperature", channels_)[6], bt91h);

  // Compute cloud liquid water amount
  cloudLiquidWater(bt19h, bt19v, bt22v, bt37h, bt37v, bt91v, bt91h, water_frac, out[0]);
}

// -----------------------------------------------------------------------------

void CLWRetMW_SSMIS::cloudLiquidWater(const std::vector<float> & bt19h,
                                      const std::vector<float> & bt19v,
                                      const std::vector<float> & bt22v,
                                      const std::vector<float> & bt37h,
                                      const std::vector<float> & bt37v,
                                      const std::vector<float> & bt91v,
                                      const std::vector<float> & bt91h,
                                      std::vector<float> & water_frac,
                                      std::vector<float> & clw) {
  ///
  /// \brief Retrieve cloud liquid water from channels 12-18 of SSMIS. Output is in
  ///        kg/m2 (mm) and bound between zero and 6.0.
  ///
  /// Reference:
  /// Weng, F., R. R. Ferraro, and N. C. Grody,2000: "Effects of AMSU cross-scan Symmetry of
  ///          brightness temperatures on  retrieval of atmospheric and surface parameters",
  ///          Ed. P. Pampaloni and S. Paloscia, VSP, Netherlands, 255-262, 2000.
  /// Yan B. and F. Weng, 'Intercalibration between Special Sensor Microwave Imager and
  ///          Sounder (SSMIS) and Special Sensor Microwave Imager (SSM/I)', TGARS Special
  ///          Issue on the DMSP SSMIS, 46, 984-995.

  const float missing = util::missingValue<float>();
  const size_t nchannels = 7;
  const size_t nlocs = bt19h.size();

  // Initialize output to missing.
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    clw[iloc] = missing;
  }

  // SSMIS information about channels and frequecies came from
  //  http://rain.atmos.colostate.edu/FCDR/ssmis.html
  //  Channels 12-18: 19.35h, 19.35v, 22.235, 37h, 37v, 91.655h, 91.655v GHz.
  //  Various coefficients and equations from GSI subroutine ret_ssmis.

  const std::vector<float> ap = {0.00424, -2.03627, -2.52875, 0.80170, -3.86053, -7.43913, 1.53650};
  const std::vector<float> bp = {1.00027, 1.00623, 0.99642, 0.99139,  1.00550, 1.03121, 0.99317};
  const std::vector<float> cp0 = {0.969, 0.969, 0.974, 0.986, 0.986, 0.988, 0.988};
  const std::vector<float> dp0 = {0.00415, 0.00473, 0.0107, 0.02612, 0.0217, 0.01383, 0.01947};

  std::vector<float> cp(nchannels), dp(nchannels);
  for (size_t i = 0; i < cp0.size(); ++i) {
    cp[i] = 1.0 / (cp0[i]*(1.0-dp0[i]));
    dp[i] = cp[i]*dp0[i];
  }

  // Setting the tax array to a relatively large number in the event that any missing data
  // will fail the IF-tests of brightness temperatures less than 285K, ensures the output
  // CLW is missing also.

  std::vector<float> bt_test(nchannels);
  std::vector<std::vector<float> > tax(nlocs, std::vector<float> (nchannels));
  std::vector<std::vector<float> > tay(nlocs, std::vector<float> (nchannels));
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    const std::vector<float> btest = {bt19h[iloc], bt19v[iloc], bt22v[iloc], bt37h[iloc],
                                      bt37v[iloc], bt91v[iloc], bt91h[iloc]};
    if (std::any_of(btest.begin(), btest.end(), [](float x){return (x <= 0.0f || x >= 300.0f);})) {
      for (size_t ich = 0; ich < nchannels; ++ich) {
        tax[iloc][ich] = 999.0f;
      }
    } else {
      tax[iloc][0] = (bt19h[iloc]*cp[1] + bt19v[iloc]*dp[0])/(cp[0]*cp[1] - dp[0]*dp[1]);
      tax[iloc][1] = (bt19h[iloc]*dp[1] + bt19v[iloc]*cp[0])/(cp[0]*cp[1] - dp[0]*dp[1]);
      tax[iloc][2] = 1.0/cp[2]*(bt22v[iloc] + dp[2]*(0.653*tax[iloc][1] + 96.6));
      tax[iloc][3] = (bt37h[iloc]*cp[4] + bt37v[iloc]*dp[3])/(cp[3]*cp[4] - dp[3]*dp[4]);
      tax[iloc][4] = (bt37h[iloc]*dp[4] + bt37v[iloc]*cp[3])/(cp[3]*cp[4] - dp[3]*dp[4]);
      tax[iloc][5] = (bt91v[iloc]*cp[6] + bt91h[iloc]*dp[5])/(cp[5]*cp[6] - dp[5]*dp[6]);
      tax[iloc][6] = (bt91v[iloc]*dp[6] + bt91h[iloc]*cp[5])/(cp[5]*cp[6] - dp[5]*dp[6]);
    }

    for (size_t ich = 0; ich < nchannels; ++ich) {
      tay[iloc][ich] =  ap[ich] + bp[ich]*tax[iloc][ich];
    }
  }

  float alg1, alg2, alg3;
  float tby1, tby3, tby4, tpwc;
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (water_frac[iloc] >= 0.99) {
      clw[iloc] = 0.0;
      alg1 = 0.0, alg2 = 0.0, alg3 = 0.0;
      tby1 = 0.0, tby3 = 0.0, tby4 = 0.0, tpwc = 0.0;
      // Try the quickest answer related to channels 2 and 3.
      alg1 = 0.0;
      if (tay[iloc][1] < 285.0 && tay[iloc][2] < 285.0) {
        alg1 = -3.20*(std::log(290.0f-tay[iloc][1]) - 2.80 - 0.42*std::log(290.0f-tay[iloc][2]));
      }
      if (alg1 > 0.70) {
        clw[iloc] = alg1;
        clw[iloc] = std::max(0.0f, std::min(alg1, 6.0f));
      } else {
        // Try the next quickest answer related to channels 3 and 5.
        alg2 = 0.0;
        if (tay[iloc][4] < 285.0 && tay[iloc][2] < 285.0) {
          alg2 = -1.66*(std::log(290.0f-tay[iloc][4]) - 2.90 - 0.349*std::log(290.0f-tay[iloc][2]));
        }
        if (alg2 > 0.28) {
          clw[iloc] = alg2;
          clw[iloc] = std::max(0.0f, std::min(alg2, 6.0f));
        } else {
          // Final test using channels 3 and 7, but we first need total precipitable water.
          tby1  =  cp[0]*tay[iloc][0] - dp[0]*tay[iloc][1];
          tby3  =  cp[2]*tay[iloc][2] - dp[2]*(0.653*tay[iloc][1] + 96.6);
          tby4  =  cp[3]*tay[iloc][3] - dp[3]*tay[iloc][4];
          tpwc = 232.89 - 0.1486*tby1 - 0.3695*tby4 - (1.8291 - 0.006193*tby3)*tby3;

          if (tpwc < 30.0) {
            alg3 = 0;
            if (tay[iloc][6] < 285.0 && tay[iloc][2] < 285.0) {
              alg3 = -0.44*(std::log(290.f-tay[iloc][6]) + 1.60
                     - 1.354*std::log(290.f-tay[iloc][2]));
              clw[iloc] = alg3;
              clw[iloc] = std::max(0.0f, std::min(alg3, 6.0f));
            }
          } else if (alg2 > 0.0) {
            clw[iloc] = alg2;
            clw[iloc] = std::max(0.0f, std::min(alg2, 6.0f));
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & CLWRetMW_SSMIS::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
