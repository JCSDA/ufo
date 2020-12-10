/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/FilterBase.h"

#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"

#include "ufo/filters/actions/FilterAction.h"
#include "ufo/filters/GenericFilterParameters.h"
#include "ufo/filters/processWhere.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------

FilterBase::FilterBase(ioda::ObsSpace & os, const FilterParametersBase & parameters,
                       std::shared_ptr<ioda::ObsDataVector<int> > flags,
                       std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : obsdb_(os), config_(parameters.toConfiguration()),
    flags_(flags), obserr_(obserr),
    allvars_(getAllWhereVariables(parameters.where)),
    filtervars_(), data_(obsdb_), prior_(false), post_(false),
    deferToPost_(parameters.deferToPost),
    whereConfig_(parameters.where),
    actionConfig_(parameters.action)
{
  oops::Log::trace() << "FilterBase contructor" << std::endl;
  ASSERT(flags);
  ASSERT(obserr);
  data_.associate(*flags_, "QCflagsData");
  data_.associate(*obserr_, "ObsErrorData");
  if (parameters.filterVariables.value() != boost::none) {
  // read filter variables
    for (const Variable &var : *parameters.filterVariables.value())
      filtervars_ += var;
  } else {
  // if no filter variables explicitly specified, filter out all variables
    filtervars_ += Variables(obsdb_.obsvariables());
  }
  FilterAction action(parameters.action);
  allvars_ += action.requiredVariables();
}

// -----------------------------------------------------------------------------

FilterBase::FilterBase(ioda::ObsSpace & os, const eckit::Configuration & config,
                       std::shared_ptr<ioda::ObsDataVector<int> > flags,
                       std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(os,
               oops::validateAndDeserialize<GenericFilterParameters>(config),
               std::move(flags),
               std::move(obserr))
{}

// -----------------------------------------------------------------------------

FilterBase::~FilterBase() {
  oops::Log::trace() << "FilterBase destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void FilterBase::preProcess() {
  oops::Log::trace() << "FilterBase preProcess begin" << std::endl;
// Cannot determine earlier when to apply filter because subclass
// constructors add to allvars
  if (allvars_.hasGroup("HofX") || allvars_.hasGroup("ObsDiag") || deferToPost_) {
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
  std::vector<bool> apply = processWhere(whereConfig_, data_);

// Allocate flagged obs indicator (false by default)
  const size_t nvars = filtervars_.nvars();
  std::vector<std::vector<bool>> flagged(nvars);
  for (size_t jv = 0; jv < flagged.size(); ++jv) flagged[jv].resize(obsdb_.nlocs());

// Apply filter
  this->applyFilter(apply, filtervars_, flagged);

// Take action
  eckit::LocalConfiguration aconf = actionConfig_;
  aconf.set("flag", this->qcFlag());
  FilterAction action(aconf);
  action.apply(filtervars_, flagged, data_, *flags_, *obserr_);

// Done
  oops::Log::trace() << "FilterBase doFilter end" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
