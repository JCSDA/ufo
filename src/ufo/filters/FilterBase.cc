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

#include "oops/util/Logger.h"

#include "ufo/filters/actions/FilterAction.h"
#include "ufo/filters/GenericFilterParameters.h"
#include "ufo/filters/processWhere.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------

FilterBase::FilterBase(ioda::ObsSpace & os,
                       const FilterParametersBaseWithAbstractActions & parameters,
                       std::shared_ptr<ioda::ObsDataVector<int> > flags,
                       std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : ObsProcessorBase(os, parameters.deferToPost, std::move(flags), std::move(obserr)),
    filtervars_(),
    whereParameters_(parameters.where),
    actionsParameters_(parameters.actions())
{
  oops::Log::trace() << "FilterBase constructor" << std::endl;

  // Identify filter variables
  if (parameters.filterVariables.value() != boost::none) {
  // read filter variables
    for (const Variable &var : *parameters.filterVariables.value())
      filtervars_ += var;
  } else {
  // if no filter variables explicitly specified, filter out all variables
    filtervars_ += Variables(obsdb_.obsvariables());
  }

  // Identify input variables required by the filter and notify user if any action except the last
  // modifies QC flags
  allvars_ += getAllWhereVariables(whereParameters_);

  const size_t numActions = actionsParameters_.size();
  for (size_t i = 0; i < numActions; ++i) {
    const std::unique_ptr<FilterActionParametersBase> &actionParameters = actionsParameters_[i];
    FilterAction action(*actionParameters);
    if (i < numActions - 1 && action.modifiesQCFlags()) {
      throw eckit::UserError("Actions modifying QC flags, such as '" +
                             actionParameters->name.value().value() + "', must not be followed by "
                             "any other actions performed by the same filter", Here());
    }
    allvars_ += action.requiredVariables();
  }
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

void FilterBase::doFilter() const {
  oops::Log::trace() << "FilterBase doFilter begin" << std::endl;

// Select locations to which the filter will be applied
  std::vector<bool> apply = processWhere(whereParameters_, data_);

// Allocate flagged obs indicator (false by default)
  const size_t nvars = filtervars_.nvars();
  std::vector<std::vector<bool>> flagged(nvars);
  for (size_t jv = 0; jv < flagged.size(); ++jv) flagged[jv].resize(obsdb_.nlocs());

// Apply filter
  this->applyFilter(apply, filtervars_, flagged);

// Take actions
  for (const std::unique_ptr<FilterActionParametersBase> &actionParameters : actionsParameters_) {
    FilterAction action(*actionParameters);
    action.apply(filtervars_, flagged, data_, this->qcFlag(), *flags_, *obserr_);
  }

// Done
  oops::Log::trace() << "FilterBase doFilter end" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
