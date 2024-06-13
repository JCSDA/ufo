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
                       std::shared_ptr<ioda::ObsDataVector<float> > obserr,
                       const VariableNameMap & nameMap)
  : ObsProcessorBase(os, parameters.deferToPost, std::move(flags), std::move(obserr)),
    filtervars_(),
    nameMap_(nameMap),
    whereParameters_(parameters.where),
    whereOperator_(parameters.whereOperator),
    actionsParameters_(parameters.actions())
{
  oops::Log::trace() << "FilterBase constructor" << std::endl;

  // Identify filter variables
  if (parameters.filterVariables.value() != boost::none) {
  // read filter variables
    for (const Variable &var : *parameters.filterVariables.value()) {
      filtervars_ += var;
      filtersimvars_ += var;
    }
  } else {
  // if no filter variables explicitly specified, filter out all variables
    filtervars_ += Variables(obsdb_.obsvariables());
    filtersimvars_ += Variables(obsdb_.assimvariables());
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
                       std::shared_ptr<ioda::ObsDataVector<float> > obserr,
                       const VariableNameMap & nameMap)
  : FilterBase(os,
               oops::validateAndDeserialize<GenericFilterParameters>(config),
               std::move(flags),
               std::move(obserr),
               nameMap)
{}

// -----------------------------------------------------------------------------

FilterBase::~FilterBase() {
  oops::Log::trace() << "FilterBase destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void FilterBase::doFilter() {
  oops::Log::trace() << "FilterBase doFilter begin" << std::endl;

// Select locations to which the filter will be applied
  std::vector<bool> apply = processWhere(whereParameters_, data_, whereOperator_);

  ufo::Variables vars;
  if (post_) {
    vars += filtersimvars_;
    if (allvars_.hasGroup("HofX")) {
      for (size_t jv = 0; jv < filtersimvars_.toOopsObsVariables().size(); ++jv) {
        if (!obsdb_.assimvariables().has(filtersimvars_.toOopsObsVariables()[jv])) {
          throw eckit::UserError("Filter variable '"
                                 + filtersimvars_.toOopsObsVariables()[jv] +
                                 "' is not a simulated variable,"
                                 " but an HofX is required", Here());
        }
      }
    }
  } else {
    vars += filtervars_;
  }

// Allocate flagged obs indicator (false by default)
  const size_t nvars = vars.nvars();
  std::vector<std::vector<bool>> flagged(nvars);
  for (size_t jv = 0; jv < flagged.size(); ++jv) flagged[jv].resize(obsdb_.nlocs());

// Apply filter
  this->applyFilter(apply, vars, flagged);

// Take actions
  for (const std::unique_ptr<FilterActionParametersBase> &actionParameters : actionsParameters_) {
    FilterAction action(*actionParameters);
    action.apply(vars, flagged, data_, this->qcFlag(), *flags_, *obserr_);
  }

// Done
  oops::Log::trace() << "FilterBase doFilter end" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
