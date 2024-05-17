/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/gnssro/QC/BackgroundCheckRONBAM.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"

#include "ufo/filters/QCflags.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// -----------------------------------------------------------------------------

BackgroundCheckRONBAM::BackgroundCheckRONBAM(ioda::ObsSpace & obsdb,
                                           const eckit::Configuration & config,
                                           std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                           std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  oops::Log::trace() << "BackgroundCheckRONBAM contructor: "
                     << "using NBAM style BackgroundCheck for GnssroBndNBAM" << std::endl;
  allvars_ += Variables(filtervars_, "HofX");
}

// -----------------------------------------------------------------------------

BackgroundCheckRONBAM::~BackgroundCheckRONBAM() {
  oops::Log::trace() << "BackgroundCheckRONBAM destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void BackgroundCheckRONBAM::applyFilter(const std::vector<bool> & apply,
                                        const Variables & filtervars,
                                        std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "BackgroundCheckRONBAM postFilter" << std::endl;

  const oops::ObsVariables observed = obsdb_.obsvariables();
  const float missing = util::missingValue<float>();

  ioda::ObsDataVector<float> obs(obsdb_, filtervars.toOopsObsVariables(), "ObsValue");
  ioda::ObsDataVector<float> bias(obsdb_, filtervars.toOopsObsVariables(), "ObsBias", false);
  ioda::ObsDataVector<float> impactparameter(obsdb_, "impactParameterRO", "MetaData");
  ioda::ObsDataVector<float> latitude(obsdb_, "latitude", "MetaData");
  ioda::ObsDataVector<float> earthradius(obsdb_, "earthRadiusCurvature", "MetaData");
  ioda::ObsDataVector<int> satelliteid(obsdb_, "satelliteIdentifier", "MetaData");
  // next line is actually background (model) virtual temperature at obs location
  ioda::ObsDataVector<float> temperature(obsdb_, "virtual_temperature", "MetaData");

  Variables varhofx(filtervars, "HofX");

  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    size_t iv = observed.find(filtervars.variable(jv).variable());

//  H(x)
    std::vector<float> hofx;
    data_.get(varhofx.variable(jv), hofx);

    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (apply[jobs] && (*flags_)[iv][jobs] == 0) {
        ASSERT(obs[jv][jobs] != missing);
        ASSERT(hofx[jobs] != missing);
        ASSERT(impactparameter[0][jobs] != missing);

        float imp = impactparameter[0][jobs]/1000.0 - earthradius[0][jobs]/1000.0;
        float lat = latitude[0][jobs]*Constants::deg2rad;
        float tmp = temperature[0][jobs];
        int satid = satelliteid[0][jobs];
//      Threshold for current observation
        float cutoff   = std::numeric_limits<float>::max();
        float cutoff1  = std::numeric_limits<float>::max();
        float cutoff2  = std::numeric_limits<float>::max();
        float cutoff3  = std::numeric_limits<float>::max();
        float cutoff4  = std::numeric_limits<float>::max();
        float cutoff12 = std::numeric_limits<float>::max();
        float cutoff23 = std::numeric_limits<float>::max();
        float cutoff34 = std::numeric_limits<float>::max();

        cutoff = 0.0;
        cutoff1  = (-4.725+0.045*imp+0.005*imp*imp);
        cutoff2  = 1.5+cos(lat);
        cutoff3  = 2.0;
        if (tmp > 240.0) cutoff3 = 0.005*tmp*tmp-2.3*tmp+266.0;
        cutoff4 = (4.0+8.0*cos(lat));
//      COSMIC-2 and Commercial providers
        if ( (satid >= 750 && satid <= 755) || (satid >= 265 && satid <= 269) ) {
           cutoff1 = cutoff1/2.0;
           cutoff3 = cutoff3/2.0;
           cutoff4 = cutoff4/2.0;
//      Default
        } else {
           cutoff1 = cutoff1*2.0/3.0;
           cutoff3 = cutoff3*2.0/3.0;
           cutoff4 = cutoff4*2.0/3.0;
        }
        cutoff12 = ( (36.0-imp)/2.0)*cutoff2 + ((imp-34.0)/2.0)*cutoff1;
        cutoff23 = ( (11.0-imp)/2.0)*cutoff3 + ((imp-9.0)/2.0)*cutoff2;
        cutoff34 = ( (6.0-imp)/2.0)*cutoff4 + ((imp-4.0)/2.0)*cutoff3;
        if (imp > 36.0) cutoff = cutoff1;
        if (imp <= 36.0 && imp > 34.0) cutoff = cutoff12;
        if (imp <= 34.0 && imp > 11.0) cutoff = cutoff2;
        if (imp <= 11.0 && imp > 9.0)  cutoff = cutoff23;
        if (imp <= 9.0  && imp > 6.0)  cutoff = cutoff3;
        if (imp <= 6.0  && imp > 4.0)  cutoff = cutoff34;
        if (imp <= 4.0) cutoff = cutoff4;
//      COSMIC-2 and Commercial providers
        if ( (satid >= 750 && satid <= 755) || (satid >= 265 && satid <= 269) ) {
          cutoff = 0.02*cutoff;
//      Default
        } else {
          cutoff = 0.03*cutoff;
        }
        ASSERT(cutoff < std::numeric_limits<float>::max() && cutoff > 0.0);

//      Apply bias correction
        float yy = obs[jv][jobs] + bias[jv][jobs];

//      NBAM style background check: if omb/o is greater than a cutoff
        if (std::abs(static_cast<float>(hofx[jobs])-yy) > yy*cutoff) {
           flagged[jv][jobs] = true; }
      }
    }
  }
}

// -----------------------------------------------------------------------------

void BackgroundCheckRONBAM::print(std::ostream & os) const {
  os << "BackgroundCheckRONBAM::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
