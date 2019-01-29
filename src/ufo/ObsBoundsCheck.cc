/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/ObsBoundsCheck.h"

#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/processWhere.h"
#include "ufo/QCflags.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, ObsBoundsCheck> >
  mkBoundChk_("Bounds Check");
// -----------------------------------------------------------------------------

ObsBoundsCheck::ObsBoundsCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config)
  : obsdb_(obsdb), config_(config), geovars_(preProcessWhere(config_))
{
  oops::Log::debug() << "ObsBoundsCheck: config = " << config_ << std::endl;
  oops::Log::debug() << "ObsBoundsCheck: geovars = " << geovars_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsBoundsCheck::~ObsBoundsCheck() {}

// -----------------------------------------------------------------------------

void ObsBoundsCheck::priorFilter(const GeoVaLs & gv) const {
  const std::string qcgrp = config_.getString("QCname");
  const std::string obgrp = "ObsValue";
  const float missing = util::missingValue(missing);

  std::vector<eckit::LocalConfiguration> bounds;
  config_.get("bounds", bounds);

  for (size_t jj = 0; jj < bounds.size(); ++jj) {
    const std::string var(bounds[jj].getString("variable"));
    const float vmin = bounds[jj].getFloat("minvalue", missing);
    const float vmax = bounds[jj].getFloat("maxvalue", missing);

    ioda::ObsDataVector<float> obs(obsdb_, var, obgrp);
    ioda::ObsDataVector<int> flags(obsdb_, var, qcgrp);

//  Select where the bounds check will apply
    std::vector<bool> apply = processWhere(obsdb_, gv, bounds[jj]);

    int ii = 0;
    for (size_t jobs = 0; jobs < obs.size(); ++jobs) {
      if (apply[jobs] && flags[jobs] == 0) {
        ASSERT(obs[jobs] != missing);
        if (vmin != missing && obs[jobs] < vmin) flags[jobs] = QCflags::bounds;
        if (vmax != missing && obs[jobs] > vmax) flags[jobs] = QCflags::bounds;
        if (flags[jobs] == QCflags::bounds) ++ii;
      }
    }

    flags.save();

    oops::Log::debug() << "ObsBoundsCheck: " << obsdb_.obsname() << " " << var
                       << " rejected " << ii << " obs." << std::endl;
  }
}

// -----------------------------------------------------------------------------

void ObsBoundsCheck::print(std::ostream & os) const {
  os << "ObsBoundsCheck: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
