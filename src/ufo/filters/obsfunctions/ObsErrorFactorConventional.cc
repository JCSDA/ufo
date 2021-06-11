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
  for (size_t iv = 0; iv < inflatevars.size(); ++iv) {
    invars_ += Variable(inflatevars[iv]+"@"+qcgrp);
  }

  // Include list of required data from MetaData
  invars_ += Variable("air_pressure@MetaData");   // observed obs pressure
  invars_ += Variable("station_id@MetaData");     // obs station ID
  invars_ += Variable("datetime@MetaData");       // obs date and time

  // Include list of required data from GeoVaLs
  invars_ += Variable("air_pressure@GeoVaLs");  // need model pressure at half
                                                // sigma levels (dimension=nsig in GSI)
}

// -----------------------------------------------------------------------------

ObsErrorFactorConventional::~ObsErrorFactorConventional() {
    oops::Log::debug() << "ObsErrorFactorCon: destructing "  << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorFactorConventional::compute(const ObsFilterData & data,
                                     ioda::ObsDataVector<float> & obserr) const {
  const int missing = util::missingValue(missing);
  static constexpr float con_g_rd = 500.0f*Constants::grav/(273.0f*Constants::rd);
  const float tiny_float = FLT_MIN;

  // Get dimensions
  const size_t nlocs = data.nlocs();
  const size_t nlevs = data.nlevs(Variable("air_pressure@GeoVaLs"));

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

  // Get MetaData of obs
  std::vector<float> ob_pressure(nlocs);
  data.get(Variable("air_pressure@MetaData"), ob_pressure);
  std::vector<std::string> ob_stationID(nlocs);
  data.get(Variable("station_id@MetaData"), ob_stationID);
  std::vector<util::DateTime> ob_datetime(nlocs);
  data.get(Variable("datetime@MetaData"), ob_datetime);

  // Due to splitting across CPUs, it is possible we have zero obs, so just return nothing.
  if (ob_datetime.empty()) {
    return;
  }

  // Get GeoVals of air pressure [Pa] in vertical column
  std::vector<std::vector<float>> prsl(nlevs, std::vector<float>(nlocs));
  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
    data.get(Variable("air_pressure@GeoVaLs"), ilev+1, prsl[ilev]);
  }

  for (size_t iv = 0; iv < varsize; ++iv) {   // Variable loop
    // Get QC flags of test variable
    std::vector<int> ob_variable_QCflag(nlocs);
    std::vector<int> ob_QCflag(nlocs);
    data.get(Variable(inflatevars[iv]+"@"+qcgrp), ob_variable_QCflag);

    // Only obs with QCflag <= qcthres will be checked in this obsFunc
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      (ob_variable_QCflag[iloc] <= qcthres && ob_variable_QCflag[iloc] != missing) ?
          ob_QCflag[iloc] = 0 : ob_QCflag[iloc] = 100;
    }

    // Define variables to be used
    std::vector<float> dprsl(nlevs-1);
    float error_factor;
    float pre1, pre2, maxpre, conpre, vmag, pdiffu, pdiffd, pdifftotal;
    int ipass = 0;
    int icount = 0;
    int ireport = 0;

    ioda::ObsSpace::RecIdxIter irec;   // Using obs grouping/sorting indices
    for ( irec = obsdb.recidx_begin(); irec != obsdb.recidx_end(); ++irec ) {   // record loop
      std:: size_t rNum = obsdb.recidx_recnum(irec);
      std::vector<std::size_t> rSort = obsdb.recidx_vector(irec);

      ireport++;
      for (size_t iloc = 0; iloc < rSort.size(); ++iloc) {   // profile loop
        icount++;

        for (size_t ilev = 0; ilev < nlevs-1; ++ilev) {
          // Background pressure level intervals [cb]
          dprsl[ilev] = (prsl[ilev][rSort[iloc]]-prsl[ilev+1][rSort[iloc]])*0.001f;
        }

        int indexlev = 0;
        for (size_t ilev = 1; ilev < nlevs-1; ++ilev) {
          if (ob_pressure[rSort[iloc]] < prsl[ilev][rSort[iloc]]) {
            indexlev = ilev;
          }
        }

        pre1 = 0.5f*dprsl[indexlev];
        pre2 = 0.02f*0.001f*prsl[indexlev][rSort[iloc]];
        maxpre = std::max(pre1, pre2);
        conpre = con_g_rd*0.001f*ob_pressure[rSort[iloc]];
        vmag = std::min(maxpre, conpre);

        pdiffu = vmag;
        pdiffd = vmag;

        if (ob_QCflag[rSort[iloc]] == 0) {
          ipass++;
          // Search obs from upper levels. If there are multiple observations
          // inside the same model interval, use the toppest one to compute
          // the pressure interval, pdiffu
          int lloc = iloc+1;
          while (lloc < rSort.size()) {
            float tmp = abs(ob_pressure[rSort[iloc]]-ob_pressure[rSort[lloc]])*0.001f;
            if (ob_QCflag[rSort[lloc]] == 0 && tmp < vmag) {
              pdiffu = tmp;
              break;
            }
            lloc++;
          }

          // Search obs from lower levels. If there are multiple observations
          // inside the same model interval, use the lowest one to compute
          // the pressure interval, pdiffd
          lloc = iloc-1;
          while (lloc >= 0) {
            float tmp = abs(ob_pressure[rSort[lloc]]-ob_pressure[rSort[iloc]])*0.001f;
            if (ob_QCflag[rSort[lloc]] == 0 && tmp < vmag) {
              pdiffd = tmp;
              break;
            }
            lloc--;
          }
        }

        // When there are multiple observations inside the same model interval, the error_factor
        // will be bigger than 1 based on the spacing of the these observations
        pdifftotal = std::max(pdiffd+pdiffu, 5.0f * tiny_float);
        error_factor = sqrt(2.0f*vmag/pdifftotal);

        // Output
        obserr[iv][rSort[iloc]] = error_factor;
      }
    }  // iloc loop
    if (icount != nlocs) {
      std::string errString = "The data should be sorted or icount and nlocs are not consistent: ";
      oops::Log::error() << errString << icount << ", "<< nlocs << std::endl;
      throw eckit::BadValue(errString);
    }  // irec loop
    oops::Log::debug() << "ObsErrorFactorCon: inflate var, # of reports, total obs, filtered obs = "
                   << inflatevars[iv] << " " << ireport << " "<< nlocs << " " << ipass << std::endl;
  }  // iv loop
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorConventional::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
