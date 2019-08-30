/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/MWCLWCheck.h"

#include <math.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/QCflags.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, MWCLWCheck> >
  mkMWCLWChk_("MWCLW Check");
// -----------------------------------------------------------------------------

MWCLWCheck::MWCLWCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                               boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                               boost::shared_ptr<ioda::ObsDataVector<float> >)
  : obsdb_(obsdb), data_(obsdb_), config_(config), geovars_(preProcessWhere(config_, "GeoVaLs")),
    diagvars_(), flags_(*flags)
{
  oops::Log::debug() << "MWCLWCheck: config = " << config_ << std::endl;
  oops::Log::debug() << "MWCLWCheck: geovars = " << geovars_ << std::endl;
}

// -----------------------------------------------------------------------------

MWCLWCheck::~MWCLWCheck() {}

// -----------------------------------------------------------------------------

void MWCLWCheck::priorFilter(const GeoVaLs & gv) const {
  const float missing = util::missingValue(missing);


  oops::Variables vars(config_);
  oops::Variables observed = obsdb_.obsvariables();


  ioda::ObsDataVector<float> obs(obsdb_, vars, "ObsValue");
// For now assume channels 1&2 are the same as those for AMSU-A and ATMS
  ioda::ObsDataVector<float> tobs1(obsdb_, "brightness_temperature_1", "ObsValue");
  ioda::ObsDataVector<float> tobs2(obsdb_, "brightness_temperature_2", "ObsValue");
  ioda::ObsDataVector<float> sza(obsdb_, "sensor_zenith_angle", "MetaData");
  ioda::ObsDataVector<float> clw(obsdb_, "cloud_liquid_water", "Diagnostic", false);

  const float clw_threshold = config_.getFloat("clw_threshold", missing);

  const float d1 = 0.754;
  const float d2 = -2.265;

// Get config
  std::cout << "MWCLWCheck: config = " << config_ << std::endl;

// Select where the bounds check will apply
  std::vector<bool> apply = processWhere(config_, data_);

// Loop over obs locations calculating CLW from observations
  for (size_t jobs = 0; jobs < obs.nlocs(); ++jobs) {
     float cossza = cos(M_PI * sza[0][jobs]/180.0);
     float d0 = 8.240 - (2.622 - 1.846*cossza)*cossza;
     if (tobs1[0][jobs] <= 284.0 && tobs2[0][jobs] <= 284.0 && tobs1[0][jobs] > 0.0
        && tobs1[0][jobs] > 0.0)
        {clw[0][jobs] = cossza*(d0 + d1*log(285.0-tobs1[0][jobs])) + d2*log(285.0-tobs2[0][jobs]);
        clw[0][jobs] = std::max(0.0f, clw[0][jobs]);
     } else {
        clw[0][jobs] = missing;
     }
// Apply CLW threshold to observations
     for (size_t jv = 0; jv < vars.size(); ++jv) {
        size_t iv = observed.find(vars[jv]);
        if (apply[jobs] && flags_[iv][jobs] == 0) {
           if (clw_threshold != missing && clw[0][jobs] != missing && clw[0][jobs] > clw_threshold)
            flags_[iv][jobs] = QCflags::clw;
      }
    }
  }
  clw.save("ObsValue");
}

// -----------------------------------------------------------------------------

void MWCLWCheck::print(std::ostream & os) const {
  os << "MWCLWCheck: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
