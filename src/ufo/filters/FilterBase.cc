/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/FilterBase.h"

#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"

#include "ufo/filters/actions/FilterAction.h"
#include "ufo/filters/processWhere.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------

FilterBase::FilterBase(ioda::ObsSpace & os, const eckit::Configuration & config,
                       boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                       boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : obsdb_(os), config_(config), flags_(flags), obserr_(obserr),
    allvars_(getAllWhereVariables(config_)),
    filtervars_(), data_(obsdb_), prior_(false), post_(false),
    defer_to_post_(false)
{
  oops::Log::trace() << "FilterBase contructor" << std::endl;
  ASSERT(flags);
  ASSERT(obserr);
  data_.associate(*flags_, "QCflagsData");
  data_.associate(*obserr_, "ObsErrorData");
  if (config_.has("filter variables")) {
  // read filter variables
    filtervars_ += Variables(config_.getSubConfigurations("filter variables"));
  } else {
  // if no filter variables explicitly specified, filter out all variables
    filtervars_ += Variables(obsdb_.obsvariables());
  }
  // Defer filter to run as a post filter even if hofx not needed.
  defer_to_post_ = config_.getBool("defer to post", false);
  eckit::LocalConfiguration aconf;
  config_.get("action", aconf);
  FilterAction action(aconf);
  allvars_ += action.requiredVariables();
}

// -----------------------------------------------------------------------------

FilterBase::~FilterBase() {
  oops::Log::trace() << "FilterBase destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void FilterBase::preProcess() {
  oops::Log::trace() << "FilterBase preProcess begin" << std::endl;
// Cannot determine earlier when to apply filter because subclass
// constructors add to allvars
  if (allvars_.hasGroup("HofX") || allvars_.hasGroup("ObsDiag") ||
      defer_to_post_ ) {
    post_ = true;
  } else {
    if (allvars_.hasGroup("GeoVaLs")) {
      prior_ = true;
    } else {
      this->doFilter();
    }
  }
  oops::Log::trace() << "FilterBase preProcess end" << std::endl;
}

// -----------------------------------------------------------------------------

void FilterBase::priorFilter(const GeoVaLs & gv) {
  oops::Log::trace() << "FilterBase priorFilter begin" << std::endl;
  if (prior_ || post_) data_.associate(gv);
  if (prior_) this->doFilter();
  oops::Log::trace() << "FilterBase priorFilter end" << std::endl;
}

// -----------------------------------------------------------------------------

void FilterBase::postFilter(const ioda::ObsVector & hofx, const ObsDiagnostics & diags) {
  oops::Log::trace() << "FilterBase postFilter begin" << std::endl;
  if (post_) {
    data_.associate(hofx, "HofX");
    data_.associate(diags);
    this->doFilter();
  }
  oops::Log::trace() << "FilterBase postFilter end" << std::endl;
}

// -----------------------------------------------------------------------------

void FilterBase::doFilter() const {
  oops::Log::trace() << "FilterBase doFilter begin" << std::endl;

// Select where the background check will apply
  std::vector<bool> apply = processWhere(config_, data_);

// Allocate flagged obs indicator (false by default)
  const size_t nvars = filtervars_.nvars();
  std::vector<std::vector<bool>> flagged(nvars);
  for (size_t jv = 0; jv < flagged.size(); ++jv) flagged[jv].resize(obsdb_.nlocs());

// Apply filter
  this->applyFilter(apply, filtervars_, flagged);

// Take action
  eckit::LocalConfiguration aconf;
  config_.get("action", aconf);
  aconf.set("flag", this->qcFlag());
  FilterAction action(aconf);
  action.apply(filtervars_, flagged, data_, *flags_, *obserr_);

// Done
  oops::Log::trace() << "FilterBase doFilter end" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
