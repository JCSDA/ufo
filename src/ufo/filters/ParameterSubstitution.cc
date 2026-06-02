/*
 * (C) Crown Copyright 2024 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ParameterSubstitution.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

namespace ufo {

ParameterSubstitution::ParameterSubstitution(ioda::ObsSpace & obsdb,
                                             const Parameters_ & parameters,
                                             std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                             std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : parameters_(parameters) {
  oops::Log::trace() << "ParameterSubstitution constructor" << std::endl;

  const eckit::LocalConfiguration sectionToRepeat = parameters.sectionToRepeat;
  const std::vector<eckit::LocalConfiguration> repetitions = parameters.repetitions.value();

  // Find out how many repetitions there are.
  // This has to be the same for each block to repeat.
  const std::size_t nreps =
    repetitions.front().getSubConfigurations(repetitions.front().keys().front()).size();
  for (const auto & variable : repetitions) {
    const std::string description = variable.keys().front();
    const std::vector<eckit::LocalConfiguration> repeatedParams =
      variable.getSubConfigurations(description);
    if (repeatedParams.size() != nreps) {
      throw eckit::UserError("Invalid number of repetitions for `" + description +
                             "`", Here());
    }
  }

  // Loop over blocks to repeat and use the nth entry in each case.
  for (size_t jrep = 0; jrep < nreps; ++jrep) {
    // Filter configuration to be modified.
    eckit::LocalConfiguration filterConf(sectionToRepeat);
    for (const auto & variable : repetitions) {
      const std::string description = variable.keys().front();
      const std::vector<eckit::LocalConfiguration> repeatedParams =
        variable.getSubConfigurations(description);
      if (!filterConf.has(description)) {
        throw eckit::UserError("The parameter `" + description +
                               "` does not exist for this filter", Here());
      }
      // Transfer repeated parameter values to the filter configuration.
      const eckit::LocalConfiguration values = repeatedParams[jrep];
      filterConf.set(description, values);
    }

    // Instantiate the component filter for this repetition.
    ObsFilterParametersWrapper params;
    params.validateAndDeserialize(filterConf);
    filters_.emplace_back(FilterFactory::create(obsdb,
                                                params.filterParameters,
                                                flags,
                                                obserr));
  }
}

void ParameterSubstitution::preProcess() {
  oops::Log::trace() << "ParameterSubstitution preProcess start" << std::endl;
  for (int jrep = 0; jrep < filters_.size(); ++jrep) {
    oops::Log::trace() << "preProcess repetition " << jrep << std::endl;
    filters_[jrep]->preProcess();
  }
  oops::Log::trace() << "ParameterSubstitution preProcess complete" << std::endl;
}

void ParameterSubstitution::priorFilter(const GeoVaLs & gv) {
  oops::Log::trace() << "ParameterSubstitution priorFilter start" << std::endl;
  for (int jrep = 0; jrep < filters_.size(); ++jrep) {
    oops::Log::trace() << "priorFilter repetition " << jrep << std::endl;
    filters_[jrep]->priorFilter(gv);
  }
  oops::Log::trace() << "ParameterSubstitution priorFilter complete" << std::endl;
}

void ParameterSubstitution::postFilter(const GeoVaLs & gv,
                                       const ioda::ObsVector & ov,
                                       const ioda::ObsVector & bv,
                                       const ObsDiagnostics & dv) {
  oops::Log::trace() << "ParameterSubstitution postFilter start" << std::endl;
  for (int jrep = 0; jrep < filters_.size(); ++jrep) {
    oops::Log::trace() << "postFilter repetition " << jrep << std::endl;
    filters_[jrep]->postFilter(gv, ov, bv, dv);
  }
  oops::Log::trace() << "ParameterSubstitution postFilter complete" << std::endl;
}

void ParameterSubstitution::checkFilterData(const FilterStage filterStage) {
  for (int jrep = 0; jrep < filters_.size(); ++jrep) {
    filters_[jrep]->checkFilterData(filterStage);
  }
}

oops::Variables ParameterSubstitution::requiredVars() const {
  return filters_.front()->requiredVars();
}

oops::ObsVariables ParameterSubstitution::requiredHdiagnostics() const {
  return filters_.front()->requiredHdiagnostics();
}

void ParameterSubstitution::print(std::ostream & os) const {
  os << "Parameter substitution: " << parameters_;
}

}  // namespace ufo
