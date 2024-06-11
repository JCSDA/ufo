/*
 * (C) Copyright 2023 NASA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorPressureCheck.h"

#include <float.h>

#include <algorithm>
#include <cmath>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/Constants.h"

#include "ufo/GeoVaLs.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {

namespace {

float grdcrd1(const float & d, const std::vector<float> & gh,
      const int & nlevs, const int & iflag) {
  int ix;
  float result;
  ASSERT(gh.size() == nlevs);

/// To determine an observation’s vertical position relative
/// to a model's pressure or geometric height levels(i.e.,gh).
///
/// The input variable 'd' represents the reported pressure or geometric height
/// from the observation, while the input variable 'gh' represents the model's
/// pressure or geometric height levels. This function returns with the relative
/// position with respect to the model's levels(gh) in a unitless float value.

  if (iflag == 1) {
  //   Case in which gh is in increasing order
    if (d <= gh[0]) {
       ix = 0;
     } else {
       ix = nlevs - 1;
       for (size_t k = 0 ; k < nlevs-1 ; ++k) {
         if (d <= gh[k]) {
           ix = k-1;
           break;
         }
       }
     }
  } else if (iflag == -1) {
  //   Case in which gh is in decreasing order
    if (d >= gh[0]) {
       ix = 0;
     } else {
       ix = nlevs - 1;
       for (size_t k = 0 ; k < nlevs-1 ; ++k) {
         if (d >= gh[k]) {
           ix = k-1;
           break;
         }
       }
     }
  }
  result = 1.0f+static_cast<float>(ix) + (d-gh[ix])/(gh[ix+1]-gh[ix]);
  return result;
}

}   // namespace

static ObsFunctionMaker<ObsErrorFactorPressureCheck> makerSteps_("ObsErrorFactorPressureCheck");

// -----------------------------------------------------------------------------

ObsErrorFactorPressureCheck::ObsErrorFactorPressureCheck(const eckit::Configuration &config)
  : invars_() {
  const float missing = util::missingValue<float>();
  // Initialize options
  options_.reset(new ObsErrorFactorPressureCheckParameters());
  options_->deserialize(config);

  const std::string inflatevars = options_->inflatevars.value();
  const float infl_coeff = options_->infl_coeff.value();

  const std::string errgrp = options_->testObserr.value();
  const std::string flaggrp = options_->testQCflag.value();

  invars_ += Variable("ObsType/"+inflatevars);
  invars_ += Variable(errgrp+"/"+inflatevars);
  invars_ += Variable(flaggrp+"/"+inflatevars);

  // Include list of required data from MetaData
  invars_ += Variable("MetaData/height");
  invars_ += Variable("MetaData/stationElevation");
  invars_ += Variable("MetaData/latitude");
  invars_ += Variable("MetaData/pressure");

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/geopotential_height");
  invars_ += Variable("GeoVaLs/surface_pressure");
  invars_ += Variable("GeoVaLs/air_pressure");
  const std::string geovar_sfc_geomz = options_->geovar_sfc_geomz.value();
  invars_ += Variable("GeoVaLs/" + geovar_sfc_geomz);

  // Include list of optional data from GeoVaLs
  if (options_->requestQSat.value()) {
    invars_ += Variable("GeoVaLs/saturation_specific_humidity");
  }
}

// -----------------------------------------------------------------------------

ObsErrorFactorPressureCheck::~ObsErrorFactorPressureCheck() {}

// -----------------------------------------------------------------------------

void ObsErrorFactorPressureCheck::compute(const ObsFilterData & data,
                                     ioda::ObsDataVector<float> & obserr) const {
  const float missing = util::missingValue<float>();

  // Get output variable size
  int nvars = obserr.nvars();
  // Ensure that only one output variable is expected.
  ASSERT(nvars == 1);

  // Get dimensions
  size_t nlocs = data.nlocs();

  const std::string inflatevars = options_->inflatevars.value();
  const float infl_coeff = options_->infl_coeff.value();
  const std::string errgrp = options_->testObserr.value();
  const std::string flaggrp = options_->testQCflag.value();
  const std::string adjusterr_name = options_->adjusterr_name.value();
  const ufo::GeoVaLs * gvals = data.getGeoVaLs();

  std::vector<int> itype(nlocs);
  data.get(Variable("ObsType/"+inflatevars), itype);

  std::vector<float> currentObserr(nlocs);
  data.get(Variable(errgrp+"/"+inflatevars), currentObserr);

  std::vector<int> qcflagdata(nlocs);
  data.get(Variable(flaggrp+"/"+inflatevars), qcflagdata);

  size_t nlevs = data.nlevs(Variable("GeoVaLs/air_pressure"));
  // Get ObsValue of height.
  std::vector<float> obs_height(nlocs);
  data.get(Variable("MetaData/height"), obs_height);
  std::vector<float> dstn(nlocs);
  data.get(Variable("MetaData/stationElevation"), dstn);
  std::vector<float> lat(nlocs);
  data.get(Variable("MetaData/latitude"), lat);
  std::vector<float> obs_pressure(nlocs);
  data.get(Variable("MetaData/pressure"), obs_pressure);

  std::vector<float> adjustErr;
  data.get(Variable(adjusterr_name+"/"+inflatevars), adjustErr);

  std::vector<float> zsges(nlocs);
  const std::string geovar_sfc_geomz = options_->geovar_sfc_geomz.value();
  data.get(Variable("GeoVaLs/" + geovar_sfc_geomz), zsges);
  std::vector<float> model_pressure_sfc(nlocs);
  data.get(Variable("GeoVaLs/surface_pressure"), model_pressure_sfc);

  std::vector<std::vector<float>> zges(nlevs, std::vector<float>(nlocs));
  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
    const size_t level = nlevs - ilev - 1;
    data.get(Variable("GeoVaLs/geopotential_height"), level, zges[ilev]);
  }

  std::vector<std::vector<float>> prsl(nlevs, std::vector<float>(nlocs));
  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
    const size_t level = nlevs - ilev - 1;
    data.get(Variable("GeoVaLs/air_pressure"), level, prsl[ilev]);
  }

  int iflag;
  double sat_specific_humidity;
  const float grav = Constants::grav;
  const float deg2rad = Constants::deg2rad;
  const float grav_equator = Constants::grav_equator;
  const float somigliana = Constants::somigliana;
  const float eccentricity_sq = Constants::eccentricity_sq;
  const float semi_major_axis = Constants::semi_major_axis;
  const float flattening = Constants::flattening;
  const float grav_ratio = Constants::grav_ratio;
  float fact, slat, sin2, termg, termr, termrg;
  float dpres, sfcchk, logobspres, logsfcpres, rlow, rhgh, drpx;
  float obserror, new_error, error_factor;
  std::vector<float> zges_mh(nlevs);
  std::vector<float> logprsl(nlevs);
  std::vector<double> q_profile(nlevs);
  bool reported_height = false;
  bool iflag_print_one = true;
  bool iflag_print_negone = true;
  std::vector<double> logprsl_double(nlevs);
  double errorx;

  for (size_t iv = 0; iv < nvars; ++iv) {   // Variable loop
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
///   To determine the observation's relative location by pressure or
///   geometric height.  Default: pressure.

      reported_height = false;

///   For some wind observations, it is determined by geometric height.
///   Surface Marine, Surface Land, Atlas Buoy and Surface MESONET(280-299)
///   reported with geometric height.
///   PIBAL(221), WIND PROFILER(228) and WIND PROFILER DECODED FROM PILOT
///   (PIBAL)(229). If the reported geometric height is missing, then
///   the reported pressure is used instead.
      if (inflatevars.compare("windEastward") == 0 ||
          inflatevars.compare("windNorthward") == 0) {
        if (itype[iloc] >= 280 && itype[iloc] < 300) {
          reported_height = true;
        } else if ((itype[iloc] >= 221 && itype[iloc] <= 229) || (itype[iloc] == 261)) {
          if (abs(obs_height[iloc]) < 1.e10) {
            reported_height = true;
          }
        }
      }

      if (reported_height) {
        fact = 0.0f;
        if (obs_height[iloc]-dstn[iloc] > 10.0f) {
           if ( obs_height[iloc]-dstn[iloc] > 1000.0f ) {
              fact = 1.0f;
           } else {
              fact = (obs_height[iloc]-dstn[iloc])/990.0f;
           }
        }
        if (itype[iloc] == 261) {
          dpres = obs_height[iloc];
        } else {
          dpres = obs_height[iloc]-(dstn[iloc]+fact*(zsges[iloc]-dstn[iloc]));
        }
        for (size_t k = 0 ; k < nlevs ; ++k) {
          zges_mh[k] = zges[k][iloc];
        }

        if ((itype[iloc] >= 223 && itype[iloc] <= 228) ||
            (itype[iloc] >= 280 && itype[iloc] < 300)) {
          slat = lat[iloc]*deg2rad;
          sin2  = sin(slat)*sin(slat);
          termg = grav_equator *
             ((1.0f+somigliana*sin2)/sqrt(1.0f-eccentricity_sq*sin2));
          termr = semi_major_axis/(1.0f + flattening + grav_ratio -
                2.0f*flattening*sin2);
          termrg = (termg/grav)*termr;

          for (size_t k = 0 ; k < nlevs ; ++k) {
            zges_mh[k] = (termr*zges[k][iloc]) / (termrg-zges[k][iloc]);
          }
        }

        ASSERT(zges_mh[nlevs-1] > zges_mh[0]);
        iflag = 1;  // in increasing order
        if (iflag_print_one) {
          iflag_print_one = false;
        }
        dpres = grdcrd1(dpres, zges_mh, nlevs, iflag);

        drpx = 0.0f;
        if ((itype[iloc] >=280 && itype[iloc] < 300) || dpres < 1.0f)
          drpx = 0.005f*abs(dstn[iloc]-zsges[iloc])*(1.0f-fact);
        if (dpres > static_cast<float>(nlevs)) drpx = 1.e6f;

        sfcchk = 0.0f;

      } else {
        logobspres = std::log(obs_pressure[iloc]);
        logsfcpres = std::log(model_pressure_sfc[iloc]);
        for (size_t k = 0 ; k < nlevs ; ++k) {
          logprsl[k] = std::log(prsl[k][iloc]);
          logprsl_double[k] = std::log(prsl[k][iloc]);
        }

        ASSERT(logprsl[0] > logprsl[nlevs-1]);
        iflag = -1;    // in decreasing order
        if (iflag_print_negone) {
          iflag_print_negone = false;
        }
        dpres = grdcrd1(logobspres, logprsl, nlevs, iflag);
        sfcchk = grdcrd1(logsfcpres, logprsl, nlevs, iflag);

        // Apply this drpx correction only to surface or surface_ship data
        if ((itype[iloc] > 179 && itype[iloc] <= 190) ||
            (itype[iloc] >= 192 && itype[iloc] < 199)) {
          if (inflatevars.compare("airTemperature") == 0 ||
            inflatevars.compare("virtualTemperature") == 0) {
            drpx = abs(1.0f-pow(model_pressure_sfc[iloc]/obs_pressure[iloc],
                 ufo::Constants::rd_over_cp))*ufo::Constants::t0c;
            if (abs(dpres) > 4.0f) {
              drpx = 1.0e10f;
            }
          } else {
            drpx = abs(1.0f-(obs_pressure[iloc]/model_pressure_sfc[iloc])) * 10.0;
          }
        } else {
          drpx = 0.0f;
        }
        if (inflatevars.compare("specificHumidity") == 0) {
            if ((itype[iloc] > 179 && itype[iloc] < 186) ||
                (itype[iloc] == 199)) dpres = 1.0;
            gvals->getAtLocation(q_profile, oops::Variable{"saturation_specific_humidity"}, iloc);
            std::reverse(q_profile.begin(), q_profile.end());

            ufo::PiecewiseLinearInterpolation vert_interp_model(logprsl_double, q_profile);
            if ((itype[iloc] >= 180) && (itype[iloc] <= 184)) {
                sat_specific_humidity = q_profile[0];
            } else {
                sat_specific_humidity = vert_interp_model(logobspres);
            }
        }
      }  // reported pressure conditional statement bracket

      rlow = std::max(sfcchk-dpres, 0.0f);
      rhgh = std::max(dpres-0.001f- static_cast<float>(nlevs)-1.0f, 0.0f);
      obserr[iv][iloc] = 1.0;
      if (qcflagdata[iloc] == 0) {
          if (inflatevars.compare("specificHumidity") == 0) {
              errorx = (adjustErr[iloc]+drpx)*sat_specific_humidity;
              errorx = std::max(0.0001, errorx);
              obserr[iv][iloc] = (errorx + (1.e6*rhgh)+(infl_coeff*rlow)) /(currentObserr[iloc]);
          } else {
              obserr[iv][iloc] = (currentObserr[iloc]+drpx+1.e6*rhgh+infl_coeff*rlow)
                                /currentObserr[iloc];
          }
          if (dpres > nlevs) obserr[iv][iloc]=1.e20f;
          if ((itype[iloc] >= 221 && itype[iloc] <= 229) && dpres < 0.0f) obserr[iv][iloc]=1.e20f;
     }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorPressureCheck::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
