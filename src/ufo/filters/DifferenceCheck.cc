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

DifferenceCheck::DifferenceCheck(ioda::ObsSpace & os, const eckit::Configuration & config,
                                 boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 boost::shared_ptr<ioda::ObsDataVector<float> >)
  : obsdb_(os), data_(obsdb_), flags_(*flags), config_(config), geovars_(), diagvars_(),
    ref_(config_.getString("reference")), val_(config_.getString("value"))
{
  oops::Log::trace() << "DifferenceCheck contructor starting" << std::endl;


  std::string var, grp;
// Reference setup
  splitVarGroup(ref_, var, grp);
  if (grp == "GeoVaLs") geovars_.push_back(var);

// Value to compare setup
  splitVarGroup(val_, var, grp);
  if (grp == "GeoVaLs") geovars_.push_back(var);
}

// -----------------------------------------------------------------------------

DifferenceCheck::~DifferenceCheck() {
  oops::Log::trace() << "DifferenceCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void DifferenceCheck::priorFilter(const GeoVaLs & gv) {
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
  data_.associate(gv);
  std::vector<bool> apply = processWhere(config_, data_);

// Get reference values and values to compare (as floats)
  std::vector<float> ref = data_.get(ref_);
  std::vector<float> val = data_.get(val_);
  ASSERT(ref.size() == val.size());

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
