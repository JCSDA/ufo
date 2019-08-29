/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/DifferenceCheck.h"

#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"

#include "ufo/filters/processWhere.h"
#include "ufo/filters/QCflags.h"
#include "ufo/GeoVaLs.h"
#include "ufo/UfoTrait.h"
#include "ufo/utils/SplitVarGroup.h"

namespace ufo {

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, DifferenceCheck> >
  makerDiff_("Difference Check");
// -----------------------------------------------------------------------------

DifferenceCheck::DifferenceCheck(ioda::ObsSpace & os, const eckit::Configuration & config,
                                 boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 boost::shared_ptr<ioda::ObsDataVector<float> >)
  : obsdb_(os), flags_(*flags), config_(config), geovars_(), diagvars_(),
    rvar_(), rgrp_(), vvar_(), vgrp_()
{
  oops::Log::trace() << "DifferenceCheck contructor starting" << std::endl;


// Reference setup
  const std::string sref = config_.getString("reference");
  splitVarGroup(sref, rvar_, rgrp_);
  if (rgrp_ == "GeoVaLs") geovars_.push_back(rvar_);

// Value to compare setup
  const std::string sval = config_.getString("value");
  splitVarGroup(sval, vvar_, vgrp_);
  if (vgrp_ == "GeoVaLs") geovars_.push_back(vvar_);
}

// -----------------------------------------------------------------------------

DifferenceCheck::~DifferenceCheck() {
  oops::Log::trace() << "DifferenceCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void DifferenceCheck::priorFilter(const GeoVaLs & gvals) const {
  oops::Log::trace() << "DifferenceCheck priorFilter" << std::endl;

  const float missing = util::missingValue(missing);
  const size_t nlocs = obsdb_.nlocs();

// min/max value setup
  float vmin = config_.getFloat("minvalue", missing);
  float vmax = config_.getFloat("maxvalue", missing);

// check for threshold and if exists, set vmin and vmax appropriately
  const float thresh = config_.getFloat("threshold", missing);
  if (thresh != missing) {
    vmin = -thresh;
    vmax = thresh;
  }

// Process "where" mask
  std::vector<bool> apply = processWhere(config_, obsdb_, &gvals);

// Get reference values (as floats)
  std::vector<float> ref(nlocs);
  if (rgrp_ == "GeoVaLs") {
    gvals.get(ref, rvar_);
  } else {
    ioda::ObsDataVector<float> tmp(obsdb_, rvar_, rgrp_);
    for (size_t jj = 0; jj < nlocs; ++jj) {
      ref[jj] = tmp[rvar_][jj];
    }
  }

// Get values to compare (as floats)
  std::vector<float> val(nlocs);
  if (vgrp_ == "GeoVaLs") {
    gvals.get(val, vvar_);
  } else {
    ioda::ObsDataVector<float> tmp(obsdb_, vvar_, vgrp_);
    for (size_t jj = 0; jj < nlocs; ++jj) {
      val[jj] = tmp[vvar_][jj];
    }
  }

// Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (apply[jobs]) {
      // check to see if one of the reference or value is missing
      if (val[jobs] == missing || ref[jobs] == missing) {
        for (size_t jv = 0; jv < flags_.nvars(); ++jv) {
          if (flags_[jv][jobs] == 0) flags_[jv][jobs] = QCflags::diffref;
        }
      } else {
// Check if difference is within min/max value range and set flag
        float diff = val[jobs] - ref[jobs];
        for (size_t jv = 0; jv < flags_.nvars(); ++jv) {
          if (vmin != missing && diff < vmin) flags_[jv][jobs] = QCflags::diffref;
          if (vmax != missing && diff > vmax) flags_[jv][jobs] = QCflags::diffref;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

void DifferenceCheck::print(std::ostream & os) const {
  os << "DifferenceCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
