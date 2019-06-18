/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/ObsBoundsCheck.h"

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
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, ObsBoundsCheck> >
  mkBoundChk_("Bounds Check");
// -----------------------------------------------------------------------------

ObsBoundsCheck::ObsBoundsCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                               boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                               boost::shared_ptr<ioda::ObsDataVector<float> >)
  : obsdb_(obsdb), config_(config), geovars_(preProcessWhere(config_)), flags_(*flags)
{
  oops::Log::debug() << "ObsBoundsCheck: config = " << config_ << std::endl;
  oops::Log::debug() << "ObsBoundsCheck: geovars = " << geovars_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsBoundsCheck::~ObsBoundsCheck() {}

// -----------------------------------------------------------------------------

void ObsBoundsCheck::priorFilter(const GeoVaLs & gv) const {
  const float missing = util::missingValue(missing);


  oops::Variables vars(config_);
  oops::Variables observed = obsdb_.obsvariables();


  ioda::ObsDataVector<float> obs(obsdb_, vars, "ObsValue");

  const float vmin = config_.getFloat("minvalue", missing);
  const float vmax = config_.getFloat("maxvalue", missing);

// Select where the bounds check will apply
  std::vector<bool> apply = processWhere(obsdb_, gv, config_);

  for (size_t jv = 0; jv < vars.size(); ++jv) {
    size_t iv = observed.find(vars[jv]);
    for (size_t jobs = 0; jobs < obs.nlocs(); ++jobs) {
      if (apply[jobs] && flags_[iv][jobs] == 0) {
        ASSERT(obs[jv][jobs] != missing);
        if (vmin != missing && obs[jv][jobs] < vmin) flags_[iv][jobs] = QCflags::bounds;
        if (vmax != missing && obs[jv][jobs] > vmax) flags_[iv][jobs] = QCflags::bounds;
      }
    }
  }
}

// -----------------------------------------------------------------------------

void ObsBoundsCheck::print(std::ostream & os) const {
  os << "ObsBoundsCheck: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
