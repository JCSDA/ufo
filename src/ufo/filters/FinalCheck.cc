/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/FinalCheck.h"

#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsSpace.h"
#include "oops/base/ObsVariables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

FinalCheck::FinalCheck(ioda::ObsSpace & obsdb, const Parameters_ &,
                       std::shared_ptr<ioda::ObsDataVector<int>> qcflags,
                       std::shared_ptr<ioda::ObsDataVector<float>> obserr)
  : ObsProcessorBase(obsdb, true /*deferToPost?*/, std::move(qcflags), std::move(obserr))
{
  oops::Log::trace() << "FinalCheck constructed" << std::endl;
}

// -----------------------------------------------------------------------------

FinalCheck::~FinalCheck() {
  oops::Log::trace() << "FinalCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void FinalCheck::doFilter() {
  oops::Log::trace() << "FinalCheck doFilter starts" << std::endl;

  const oops::ObsVariables &derived = obsdb_.derived_obsvariables();
  for (size_t jv = 0; jv < derived.size(); ++jv) {
    if (!obsdb_.has("ObsValue", derived[jv]))
      throw eckit::UnexpectedState(
          "All filters have been run, but the derived simulated variable " + derived[jv] +
          " can't be found either in the ObsValue or the DerivedObsValue group", Here());
  }

  // Set the QC flag to missing for any observations that haven't been rejected yet,
  // but have missing error estimates.
  const float missing = util::missingValue<float>();
  for (size_t jv = 0; jv < obsdb_.obsvariables().size(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if ((*flags_)[jv][jobs] == QCflags::pass && (*obserr_)[jv][jobs] == missing) {
        (*flags_)[jv][jobs] = QCflags::missing;
      }
    }
  }

  // Set the QC flag to reject observations that have been processed but are not to
  // be assimilated.
  for (size_t jv = 0; jv < obsdb_.obsvariables().size(); ++jv) {
    if (!obsdb_.assimvariables().has(obsdb_.obsvariables()[jv])) {
      for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
        (*flags_)[jv][jobs] = QCflags::processed;
      }
    }
  }

  oops::Log::trace() << "FinalCheck doFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

void FinalCheck::print(std::ostream & os) const {
  os << "Final Check" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
