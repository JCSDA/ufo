/*
 * (C) Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

/* 1-D var qc
 *   J(x) = (x-xb)T B-1 (x-xb) + (y-H(x))T R-1 (y-H(x))
 *   Code adapted from Met Office OPS System
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "ufo/GeoVaLs.h"
#include "ufo/rttovonedvarcheck/RTTOVOneDVarCheck.h"
#include "ufo/rttovonedvarcheck/RTTOVOneDVarCheck.interface.h"

#include "eckit/config/Configuration.h"

#include "oops/util/IntSetParser.h"

namespace ufo {

// -----------------------------------------------------------------------------

RTTOVOneDVarCheck::RTTOVOneDVarCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                 boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr), config_(config)
{
  oops::Log::debug() << "RTTOVOneDVarCheck contructor starting" << std::endl;

  // parse channels from the config and create variable names
  const oops::Variables & variables = obsdb.obsvariables();
  channels_ = variables.channels();

  // Choose when to apply filter - this is a temporary fix
  // to run as a post filter
  if (config_.has("applyfilter")) {
    std::vector<eckit::LocalConfiguration> testvarconf;
    config_.get("applyfilter", testvarconf);
    allvars_ += ufo::Variables(testvarconf);
  }

  // Setup fortran object
  const eckit::Configuration * conf = &config_;
  ufo_rttovonedvarcheck_create_f90(key_, obsdb, conf, channels_.size(), channels_[0],
                                   RTTOVOneDVarCheck::qcFlag());

  oops::Log::debug() << "RTTOVOneDVarCheck contructor complete. " << std::endl;
}

// -----------------------------------------------------------------------------

RTTOVOneDVarCheck::~RTTOVOneDVarCheck() {
  ufo_rttovonedvarcheck_delete_f90(key_);
  oops::Log::trace() << "RTTOVOneDVarCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void RTTOVOneDVarCheck::applyFilter(const std::vector<bool> & apply,
                               const Variables & filtervars,
                               std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "RTTOVOneDVarCheck Filter starting" << std::endl;

// Get GeoVaLs
  const ufo::GeoVaLs * gvals = data_.getGeoVaLs();

// Create oops variable
  oops::Variables variables = filtervars.toOopsVariables();

// Convert apply to char for passing to fortran
// needed for channel selection
  std::vector<char> apply_char(apply.size(), 'F');
  for (size_t i = 0; i < apply_char.size(); i++) {
    if (apply[i]) {apply_char[i]='T';}
  }

// Save qc flags to database for retrieval in fortran - needed for channel selection
  flags_->save("FortranQC");    // temporary measure as per gnss qc

// Pass it all to fortran
  const eckit::Configuration * conf = &config_;
  ufo_rttovonedvarcheck_apply_f90(key_, variables, gvals->toFortran(),
                                  apply_char.size(), apply_char[0]);

// Read qc flags from database
  flags_->read("FortranQC");    // temporary measure as per gnss qc

  oops::Log::trace() << "RTTOVOneDVarCheck Filter complete" << std::endl;
}

// -----------------------------------------------------------------------------

void RTTOVOneDVarCheck::print(std::ostream & os) const {
  os << "RTTOVOneDVarCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
