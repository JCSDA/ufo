/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "ufo/onedvarcheck/OneDVarCheck.h"
#include "ufo/onedvarcheck/OneDVarFortran.interface.h"
#include "ufo/GeoVaLs.h"

#include "eckit/config/Configuration.h"

#include "oops/util/IntSetParser.h"

namespace ufo {

// -----------------------------------------------------------------------------

OneDVarCheck::OneDVarCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                 boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr), config_(config),
    gv_(NULL), cost_converge_(0.01)
{
  oops::Log::debug() << "OneDVarCheck contructor starting" << std::endl;

  // parse channels from the config and create variable names
  const oops::Variables & variables = obsdb.obsvariables();
  if (config_.has("channels")) {
    channels_ = variables.channels();
  }

  // Choose when to apply filter - this is a temporary fudge
  if (config_.has("applyfilter")) {
    std::vector<eckit::LocalConfiguration> testvarconf;
    config_.get("applyfilter", testvarconf);
    allvars_ += ufo::Variables(testvarconf);
  }

  // Setup fortran object
  const eckit::Configuration * conf = &config_;
  ufo_onedvarfortran_create_f90(key_, obsdb, conf, channels_.size(), channels_[0]);

  oops::Log::debug() << "OneDVarCheck contructor complete. " << std::endl;

}

// -----------------------------------------------------------------------------

OneDVarCheck::~OneDVarCheck() {
  ufo_onedvarfortran_delete_f90(key_);
  oops::Log::trace() << "OneDVarCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void OneDVarCheck::applyFilter(const std::vector<bool> & apply,
                               const Variables & filtervars,
                               std::vector<std::vector<bool>> & flagged) const {

  oops::Log::trace() << "OneDVarCheck Filter starting" << std::endl;

  /* 1-D var qc
     J(x) = (x-xb)T B-1 (x-xb) + (y-H(x))T R-1 (y-H(x))
     Code adapted from Variational.h and IncrementalAssimilation.h
  */

// Get GeoVaLs
  const ufo::GeoVaLs * gvals = data_.getGeoVaLs();

// Create oops variable
  oops::Variables variables = filtervars.toOopsVariables();
  oops::Log::trace() << "OneDVarCheck variables = " << variables << std::endl;

// Save qc flags to database for retrieval in fortran - needed for channel selection
  flags_->save("FortranQC");    // should pass values to fortran properly

// Pass it all to fortran
  const eckit::Configuration * conf = &config_;
  ufo_onedvarfortran_post_f90(key_, variables, gvals->toFortran());

// Read qc flags from database
  flags_->read("FortranQC");    // should get values from fortran properly

// Print output flags_
  oops::Log::trace() << "OneDVarCheck flags_ = " << *flags_ << std::endl;

  oops::Log::trace() << "OneDVarCheck Filter complete" << std::endl;

}

// -----------------------------------------------------------------------------

void OneDVarCheck::print(std::ostream & os) const {
  os << "OneDVarCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
