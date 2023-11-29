/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorConventional.h"

#include <float.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorConventional> makerSteps_("ObsErrorFactorConventional");

// -----------------------------------------------------------------------------

ObsErrorFactorConventional::ObsErrorFactorConventional(const eckit::Configuration &config)
  : invars_() {
  oops::Log::debug() << "ObsErrorFactorConventional: config = " << config << std::endl;

  // Initialize options
  options_.reset(new ObsErrorFactorConventionalParameters());
  options_->deserialize(config);

  // Variable to be inflated
  const std::vector<std::string> inflatevars = options_->inflatevars.value();

  // QC flags come from files or by default from filters
  const std::string qcgrp = options_->testQCflag.value();

  // Include list of required qc flags for test variable
  for (size_t ivar = 0; ivar < inflatevars.size(); ++ivar) {
    invars_ += Variable(qcgrp+"/"+inflatevars[ivar]);
  }

  // Include list of required data from MetaData
  invars_ += Variable(options_->pressureFullName);   // observed obs pressure
  invars_ += Variable("MetaData/stationIdentification");     // obs station ID
  invars_ += Variable("MetaData/dateTime");       // obs date and time

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/air_pressure");  // need model pressure at half
                                                // sigma levels (dimension=nsig in GSI)
}

// -----------------------------------------------------------------------------

ObsErrorFactorConventional::~ObsErrorFactorConventional() {
    oops::Log::debug() << "ObsErrorFactorCon: destructing "  << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorFactorConventional::compute(const ObsFilterData & data,
                                     ioda::ObsDataVector<float> & obserr) const {
  const int missing = util::missingValue<int>();
  static constexpr float con_g_rd = 500.0f*Constants::grav/(273.0f*Constants::rd);
  const float tiny_float = FLT_MIN;
  // TODO(HuiShao): replace the earth radius at equator after matching GSI
  // with the mean Earth radius in Constants.
  static constexpr double rearth_equator = 6.37813662e+06;  // equatorial earth radius (m)

  // Get dimensions
  const size_t nlocs = data.nlocs();
  const size_t nlevs = data.nlevs(Variable("GeoVaLs/air_pressure"));

  // Get obs space
  auto & obsdb = data.obsspace();

  // Get output variable size
  int varsize = obserr.nvars();

  // Get yaml options
  const std::vector<std::string> inflatevars = options_->inflatevars.value();
  const std::string qcgrp = options_->testQCflag.value();

  // QC flags can be either from PreQC or filter flags
  // qcthres is 0 if using flter flags
  // qcthres is 3 (default) or defined through yaml if using PreQC
  int qcthres = 0;
  if (qcgrp == "PreQC") qcthres = options_->qcthreshold.value().value_or(3);

  // threshold value for horizontal distance check (m)
  float distthres = options_->distthreshold.value().value_or(0.0f);

  // Get MetaData of obs
  std::vector<float> ob_pressure(nlocs);
  data.get(Variable(options_->pressureFullName), ob_pressure);
  std::vector<std::string> ob_stationID(nlocs);
  data.get(Variable("MetaData/stationIdentification"), ob_stationID);
  std::vector<util::DateTime> ob_datetime(nlocs);
  data.get(Variable("MetaData/dateTime"), ob_datetime);

  std::vector<float> ob_lat(nlocs);
  std::vector<float> ob_lon(nlocs);
  if (distthres > 0.) {
    data.get(Variable("MetaData/latitude"), ob_lat);
    data.get(Variable("MetaData/longitude"), ob_lon);
  }

  // Due to splitting across CPUs, it is possible we have zero obs, so just return nothing.
  if (ob_datetime.empty()) {
    return;
  }

  // Get GeoVaLs of air pressure [Pa] in vertical column
  std::vector<std::vector<float>> prsl(nlevs, std::vector<float>(nlocs));
  for (size_t geolev = 0; geolev < nlevs; ++geolev) {
    data.get(Variable("GeoVaLs/air_pressure"), geolev, prsl[geolev]);
  }

  for (size_t ivar = 0; ivar < varsize; ++ivar) {   // Variable loop
    // Get QC flags of test variable
    std::vector<int> ob_variable_QCflag(nlocs);
    std::vector<int> ob_QCflag(nlocs);
    data.get(Variable(qcgrp+"/"+inflatevars[ivar]), ob_variable_QCflag);

    // Only obs with QCflag <= qcthres will be checked in this obsFunc
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      (ob_variable_QCflag[iloc] <= qcthres && ob_variable_QCflag[iloc] != missing) ?
          ob_QCflag[iloc] = 0 : ob_QCflag[iloc] = 100;
    }

    // Define variables to be used
    std::vector<float> dprsl(nlevs-1);
    float error_factor;
    float pre1, pre2, maxpre, conpre, vmag, pdiffu, pdiffd, pdifftotal, tmp;
    float dist_x, dist_y, dist_z, dist, dist_chord, central_angle;
    float rlat_this, rlon_this, rlat_next, rlon_next;
    int nextPoint;
    int passCount = 0;
    int pointCount = 0;
    int profCount = 0;

    ioda::ObsSpace::RecIdxIter irec;   // Using obs grouping/sorting indices
    // record (profile) loop
    for ( irec = obsdb.recidx_begin(); irec != obsdb.recidx_end(); ++irec ) {
      std:: size_t rNum = obsdb.recidx_recnum(irec);
      std::vector<std::size_t> rSort = obsdb.recidx_vector(irec);

      profCount++;

      // loop over the vertical obs profile
      for (size_t thisPoint = 0; thisPoint < rSort.size(); ++thisPoint) {
        float thislevpress = 100000000.0f;
        float maxlev = 0.0f;
        int thislev = 0;
        pointCount++;

        for (size_t geolev = 0; geolev < nlevs-1; ++geolev) {
          // Background pressure level intervals [cb]
          dprsl[geolev] = \
             abs((prsl[geolev][rSort[thisPoint]]-prsl[geolev+1][rSort[thisPoint]]))*0.001f;
        }

        for (size_t geolev = 1; geolev < nlevs-1; ++geolev) {
          if (ob_pressure[rSort[thisPoint]] < prsl[geolev][rSort[thisPoint]] &&
            prsl[geolev][rSort[thisPoint]] < thislevpress) {
            thislev = geolev;
            thislevpress = prsl[thislev][rSort[thisPoint]];
          }
          if (prsl[geolev][rSort[thisPoint]] > prsl[maxlev][rSort[thisPoint]]) {
            maxlev = geolev;
          }
        }

        if (thislevpress > prsl[maxlev][rSort[thisPoint]]) {
          thislevpress = prsl[maxlev][rSort[thisPoint]];
          thislev = maxlev;
        }

        if (distthres > 0.) {
          rlat_this = ob_lat[rSort[thisPoint]]*Constants::deg2rad;
          rlon_this = ob_lon[rSort[thisPoint]]*Constants::deg2rad;
        }

        // pre1: half of the vertical pressure inverval
        // pre2: 2% of the value of pressure
        // conpre: the pressure interval (delt_p) for a height inverval (delt_h)=500m
        //         delt_p=rho*g*h=(p/rT)*g*delt_h where delt_h=500m, T=273.K
        // vmag: observations within a pressure interval <= vmag will be inflated.
        //       = pre1 or pre2, whichever is bigger, if not go beyond 500m in terms
        //         of layer depth in meter, or
        //       = the pressure interval with 500m in depth assuming T=273K
        // The pressures for this algorithm are assumed to be in centibars,
        // hence the factor of 0.001.
        // dprsl is already converted from Pa to centibar.
        pre1 = 0.5f*dprsl[thislev];
        pre2 = 0.02f*0.001f*thislevpress;
        maxpre = std::max(pre1, pre2);
        conpre = con_g_rd*0.001f*ob_pressure[rSort[thisPoint]];
        vmag = std::min(maxpre, conpre);

        pdiffu = vmag;
        pdiffd = vmag;

        if (ob_QCflag[rSort[thisPoint]] == 0) {
          passCount++;
          // Search obs from upper levels. If there are multiple observations
          // inside the same model interval, use the topmost one to compute
          // the pressure interval, pdiffu
          nextPoint = thisPoint+1;
          while (nextPoint < rSort.size()) {
            tmp = abs(ob_pressure[rSort[thisPoint]]-ob_pressure[rSort[nextPoint]])*0.001f;
            if (ob_QCflag[rSort[nextPoint]] == 0 && tmp < vmag) {
              pdiffu = tmp;
              if (distthres > 0.0f) {
                rlat_next = ob_lat[rSort[nextPoint]]*Constants::deg2rad;
                rlon_next = ob_lon[rSort[nextPoint]]*Constants::deg2rad;
                dist_x = cos(rlat_this)*cos(rlon_this)-cos(rlat_next)*cos(rlon_next);
                dist_y = cos(rlat_this)*sin(rlon_this)-cos(rlat_next)*sin(rlon_next);
                dist_z = sin(rlat_this)-sin(rlat_next);
                dist = std::min(1.0f, sqrt(dist_x*dist_x+dist_y*dist_y+dist_z*dist_z));
                central_angle = 2.0f*asin(dist/2.0f);
                dist_chord = rearth_equator*central_angle;
                if (dist_chord > distthres) pdiffu = vmag;
              }
              break;
            }
            nextPoint++;
          }

          // Search obs from lower levels. If there are multiple observations
          // inside the same model interval, use the lowest one to compute
          // the pressure interval, pdiffd
          nextPoint = thisPoint-1;
          while (nextPoint >= 0) {
            tmp = abs(ob_pressure[rSort[nextPoint]]-ob_pressure[rSort[thisPoint]])*0.001f;
            if (ob_QCflag[rSort[nextPoint]] == 0 && tmp < vmag) {
              pdiffd = tmp;
              if (distthres > 0.0f) {
                rlat_next = ob_lat[rSort[nextPoint]]*Constants::deg2rad;
                rlon_next = ob_lon[rSort[nextPoint]]*Constants::deg2rad;
                dist_x = cos(rlat_this)*cos(rlon_this)-cos(rlat_next)*cos(rlon_next);
                dist_y = cos(rlat_this)*sin(rlon_this)-cos(rlat_next)*sin(rlon_next);
                dist_z = sin(rlat_this)-sin(rlat_next);
                dist = std::min(1.0f, sqrt(dist_x*dist_x+dist_y*dist_y+dist_z*dist_z));
                central_angle = 2.0f*asin(dist/2.0f);
                dist_chord = rearth_equator*central_angle;
                if (dist_chord > distthres) pdiffd = vmag;
              }
              break;
            }
            nextPoint--;
          }
        }

        // When there are multiple observations inside the same model interval, the error_factor
        // will be bigger than 1 based on the spacing of the these observations
        pdifftotal = std::max(pdiffd+pdiffu, 5.0f * tiny_float);
        error_factor = sqrt(2.0f*vmag/pdifftotal);

       // Output
       obserr[ivar][rSort[thisPoint]] = error_factor;
      }
    }  // thisPoint (observations for single profile) loop
    if (pointCount != nlocs) {
      std::string errString = "The data should be sorted or total number of observations "
                              "after sorting is not consistent with nlocs: ";
      oops::Log::error() << errString << pointCount << ", "<< nlocs << std::endl;
      throw eckit::BadValue(errString);
    }  // irec (profile) loop
    oops::Log::debug() << "ObsErrorFactorCon: inflate var, # of profiles, total obs, "
            "filtered obs = " << inflatevars[ivar] << " " << profCount << " "<< nlocs
            << " " << passCount << std::endl;
  }  // ivar (variable) loop
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorConventional::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
