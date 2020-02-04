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
  std::string chlist = config.getString("channels");
  std::set<int> channels = oops::parseIntSet(chlist);
  channels_.reserve(channels.size());
  for (const int jj : channels) {
    channels_.push_back(jj);
  }
  oops::Log::info() << "OneDVarCheck channels: " << channels_ << std::endl;

  // Read in cost convergence criteria to change default
  if (config_.has("variational.cost_convergence")) {
    config_.get("variational.cost_convergence", cost_converge_);
    oops::Log::info() << "OneDVarCheck updating cost convergence to: "
                      << cost_converge_ << std::endl;
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

// Get hofx
  Variables varhofx(filtervars_, "HofX");
  oops::Log::trace() << "OneDVarCheck Filter hofx vars = " << varhofx << std::endl;

  for (std::size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
    oops::Log::trace() << "OneDVarCheck Filter loop index = " << jobs << std::endl;

    std::vector<float> hofx;
    data_.get(varhofx.variable(jobs), hofx);
    oops::Log::trace() << "OneDVarCheck hofx = " << hofx << std::endl;

  }

// Pass it all to fortran
//  const eckit::Configuration * conf = &config_;
//  ufo_onedvarfortran_post_f90(key_, hofx->nvars(), hofx->nlocs(), hofx->toFortran(),
//                              hofx->varnames().toFortranBetter(), gvals->toFortran(),
//                              conf);

oops::Log::trace() << "OneDVarCheck Filter complete" << std::endl;

}

// -----------------------------------------------------------------------------

void OneDVarCheck::print(std::ostream & os) const {
  os << "OneDVarCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
