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

#include "ioda/distribution/Accumulator.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"

#include "ufo/filters/actions/FilterAction.h"
#include "ufo/filters/GenericFilterParameters.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/QCflags.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------

FilterBase::FilterBase(ioda::ObsSpace & os,
                       const FilterParametersBaseWithAbstractActions & parameters,
                       ioda::ObsDataVector<int> & flags,
                       ioda::ObsDataVector<float> & obserr,
                       const VariableNameMap & nameMap)
  : ObsProcessorBase(os, parameters.deferToPost, flags, obserr),
    filtervars_(),
    nameMap_(nameMap),
    whereParameters_(parameters.where),
    whereOperator_(parameters.whereOperator),
    actionsParameters_(parameters.actions()),
    identifier_(parameters.identifier.value())
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
                       ioda::ObsDataVector<int> & flags,
                       ioda::ObsDataVector<float> & obserr,
                       const VariableNameMap & nameMap)
  : FilterBase(os,
               oops::validateAndDeserialize<GenericFilterParameters>(config),
               flags,
               obserr,
               nameMap)
{}

// -----------------------------------------------------------------------------

FilterBase::~FilterBase() {
  oops::Log::trace() << "FilterBase destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void FilterBase::doFilter() {
  oops::Log::trace() << "FilterBase doFilter start" << std::endl;

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

// Log flagged count if identifier is specified and logging enabled
  if (identifier_ != boost::none && identifier_->logging) {
    this->outputIdentifierLogging(vars, nvars, flagged);
  }

// Write DiagnosticFlags if identifier is specified and diagnostic flag enabled
  if (identifier_ != boost::none && identifier_->diagnosticFlag) {
    this->writeIdentifierDiagnosticFlags(vars, nvars, flagged, false);
  }

// Write DiagnosticFlags (newly flagged only) if enabled
  if (identifier_ != boost::none && identifier_->diagnosticFlagNew) {
    this->writeIdentifierDiagnosticFlags(vars, nvars, flagged, true);
  }

// Take actions
  for (const std::unique_ptr<FilterActionParametersBase> &actionParameters : actionsParameters_) {
    FilterAction action(*actionParameters);
    action.apply(vars, flagged, data_, this->qcFlag(), flags_, obserr_);
  }

// Done
  oops::Log::trace() << "FilterBase doFilter complete" << std::endl;
}

// -----------------------------------------------------------------------------

void FilterBase::outputIdentifierLogging(
    const Variables & vars, size_t nvars,
    const std::vector<std::vector<bool>> & flagged) const {
  std::unique_ptr<ioda::Accumulator<std::vector<size_t>>> accumulator =
      obsdb_.distribution()->createAccumulator<size_t>(nvars * 3);
  for (size_t jv = 0; jv < nvars; ++jv) {
    const size_t iv = flags_.varnames().find(vars.variable(jv).variable());
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      accumulator->addTerm(jobs, jv * 3 + 2, 1);  // total count
      if (flagged[jv][jobs]) {
        accumulator->addTerm(jobs, jv * 3, 1);    // flagged count
        if (flags_[iv][jobs] == QCflags::pass)
          accumulator->addTerm(jobs, jv * 3 + 1, 1);  // newly flagged count
      }
    }
  }
  const std::vector<size_t> counts = accumulator->computeResult();
  for (size_t jv = 0; jv < nvars; ++jv) {
    oops::Log::info() << "FilterID [" << identifier_->name.value()
                      << "] " << obsdb_.obsname()
                      << " outerLoop" << getIteration()
                      << " " << vars.variable(jv).fullName()
                      << ": flagged " << counts[jv * 3]
                      << " (newly " << counts[jv * 3 + 1] << ")"
                      << " out of "
                      << counts[jv * 3 + 2] << " obs" << std::endl;
  }
}

// -----------------------------------------------------------------------------

void FilterBase::writeIdentifierDiagnosticFlags(
    const Variables & vars, size_t nvars,
    const std::vector<std::vector<bool>> & flagged, bool newlyFlaggedOnly) const {
  const std::string & filterId = identifier_->name.value();
  const std::string suffix = newlyFlaggedOnly ? "_new" : "";
  const std::vector<std::string> dimList{"Location"};
  for (size_t jv = 0; jv < nvars; ++jv) {
    const size_t iv = flags_.varnames().find(vars.variable(jv).variable());
    std::vector<bool> diagFlag(obsdb_.nlocs(), false);
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (flagged[jv][jobs]) {
        if (!newlyFlaggedOnly || flags_[iv][jobs] == QCflags::pass) {
          diagFlag[jobs] = true;
        }
      }
    }
    obsdb_.put_db("DiagnosticFlags/" + filterId + std::to_string(getIteration()) + suffix,
                  vars.variable(jv).variable(),
                  diagFlag, dimList);
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
