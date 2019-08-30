/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/MWSeaIceCheck.h"

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
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, MWSeaIceCheck> >
  mkSeaIceChk_("MW SeaIce Check");
// -----------------------------------------------------------------------------

MWSeaIceCheck::MWSeaIceCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                               boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                               boost::shared_ptr<ioda::ObsDataVector<float> >)
  : obsdb_(obsdb), data_(obsdb_), config_(config), geovars_(preProcessWhere(config_, "GeoVaLs")),
    diagvars_(), flags_(*flags)
{
  oops::Log::debug() << "MWSeaIceCheck: config = " << config_ << std::endl;
  oops::Log::debug() << "MWSeaIceCheck: geovars = " << geovars_ << std::endl;
}

// -----------------------------------------------------------------------------

MWSeaIceCheck::~MWSeaIceCheck() {}

// -----------------------------------------------------------------------------

void MWSeaIceCheck::priorFilter(const GeoVaLs & gv) const {
  data_.associate(gv);
  const float missing = util::missingValue(missing);

  oops::Variables vars(config_);
  oops::Variables observed = obsdb_.obsvariables();


  ioda::ObsDataVector<float> obs(obsdb_, vars, "ObsValue");
// this function was designed specifically for AMSU-A channels 1, 2, and 15
// may apply for ATMS but has not been tested
// will also flag majority of land points
  ioda::ObsDataVector<float> ObsCh1(obsdb_, "brightness_temperature_1", "ObsValue");
  ioda::ObsDataVector<float> ObsCh2(obsdb_, "brightness_temperature_2", "ObsValue");
  ioda::ObsDataVector<float> ObsCh15(obsdb_, "brightness_temperature_15", "ObsValue");

  const float SeaIce_threshold = config_.getFloat("limits.SeaIce_threshold", missing);

// Select which channels will have the sea ice check applied
  std::vector<bool> apply = processWhere(config_, data_);

  for (size_t jobs = 0; jobs < obs.nlocs(); ++jobs) {
    float siw = missing;
    if (ObsCh1[0][jobs] != missing && ObsCh2[0][jobs] != missing && ObsCh15[0][jobs] != missing) {
      siw = -113.2 + ((2.41 - (0.0049 * ObsCh1[0][jobs])) * ObsCh1[0][jobs]) +
            (0.454 * ObsCh2[0][jobs]) - ObsCh15[0][jobs];
    }
    for (size_t jv = 0; jv < vars.size(); ++jv) {
      size_t iv = observed.find(vars[jv]);
      if (apply[jobs] && flags_[iv][jobs] == 0) {
        if (SeaIce_threshold != missing && siw != missing && siw > SeaIce_threshold) {
            flags_[iv][jobs] = QCflags::seaice;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

void MWSeaIceCheck::print(std::ostream & os) const {
  os << "MWSeaIceCheck: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
