/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/QCmanager.h"

#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/distribution/Accumulator.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// Presets for QC filters could be performed in a function outside of any class.
// We keep them as a filter for now. The main reason for this is to be able to use
// the factory for models not in UFO/IODA.

namespace {

/// At each location, set the QC flag to QCflags::missing if the current QC flag is invalid
/// or if the ObsValue is missing.
void updateQCFlags(const std::vector<float> *obsValues, std::vector<int>& qcflags) {
  const float rmiss = util::missingValue(rmiss);
  const int imiss = util::missingValue(imiss);
  for (size_t jobs = 0; jobs < qcflags.size(); ++jobs) {
    if (qcflags[jobs] == imiss || !obsValues || (*obsValues)[jobs] == rmiss) {
      qcflags[jobs] = QCflags::missing;
    }
  }
}

}  // namespace

// -----------------------------------------------------------------------------

QCmanager::QCmanager(ioda::ObsSpace & obsdb, const eckit::Configuration & /*config*/,
                     std::shared_ptr<ioda::ObsDataVector<int> > qcflags,
                     std::shared_ptr<ioda::ObsDataVector<float> > /*obserr*/)
  : obsdb_(obsdb), nogeovals_(), nodiags_(), flags_(qcflags)
{
  oops::Log::trace() << "QCmanager::QCmanager starting" << std::endl;

  ASSERT(qcflags);

  const oops::Variables &allSimulatedVars = obsdb.obsvariables();
  const oops::Variables &initialSimulatedVars = obsdb.initial_obsvariables();
  const oops::Variables &derivedSimulatedVars = obsdb.derived_obsvariables();

  ASSERT(allSimulatedVars.size() == initialSimulatedVars.size() + derivedSimulatedVars.size());
  ASSERT(flags_->nvars() == allSimulatedVars.size());
  ASSERT(flags_->nlocs() == obsdb_.nlocs());

  const float rmiss = util::missingValue(rmiss);
  const int imiss = util::missingValue(imiss);

  const ioda::ObsDataVector<float> obs(obsdb, initialSimulatedVars, "ObsValue");

  // Iterate over initial simulated variables
  for (size_t jv = 0; jv < initialSimulatedVars.size(); ++jv) {
    const ioda::ObsDataRow<float> &currentObsValues = obs[jv];
    ioda::ObsDataRow<int> &currentQCFlags = (*qcflags)[obs.varnames()[jv]];
    updateQCFlags(&currentObsValues, currentQCFlags);
  }

  // Iterate over derived simulated variables and if they don't exist yet, set their QC flags to
  // 'missing'.
  for (size_t jv = 0; jv < derivedSimulatedVars.size(); ++jv) {
    ioda::ObsDataRow<int> &currentQCFlags = (*qcflags)[derivedSimulatedVars[jv]];
    if (!obsdb.has("ObsValue", derivedSimulatedVars[jv])) {
      updateQCFlags(nullptr, currentQCFlags);
    } else {
      std::vector<float> currentObsValues(obsdb_.nlocs());
      obsdb_.get_db("ObsValue", derivedSimulatedVars[jv], currentObsValues);
      updateQCFlags(&currentObsValues, currentQCFlags);
    }
  }

  oops::Log::trace() << "QCmanager::QCmanager done" << std::endl;
}

// -----------------------------------------------------------------------------

void QCmanager::postFilter(const ioda::ObsVector & hofx,
                           const ioda::ObsVector & /*bias*/,
                           const ObsDiagnostics & /*diags*/) {
  oops::Log::trace() << "QCmanager postFilter" << std::endl;

  const double missing = util::missingValue(missing);
  const oops::Variables &allSimulatedVars = obsdb_.obsvariables();

  for (size_t jv = 0; jv < allSimulatedVars.size(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      size_t iobs = allSimulatedVars.size() * jobs + jv;
      if ((*flags_)[jv][jobs] == 0 && hofx[iobs] == missing) {
        (*flags_)[jv][jobs] = QCflags::Hfailed;
      }
    }
  }
  oops::Log::trace() << "QCmanager postFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

QCmanager::~QCmanager() {
  oops::Log::trace() << "QCmanager::~QCmanager starting" << std::endl;
  oops::Log::info() << *this;
  oops::Log::trace() << "QCmanager::~QCmanager done" << std::endl;
}

// -----------------------------------------------------------------------------

void QCmanager::print(std::ostream & os) const {
  const std::vector<std::pair<int, const char*>> cases{
    // Special cases reported using dedicated code
    {QCflags::pass, nullptr},
    {76, nullptr},  // } The numbers of observations with these two flags
    {77, nullptr},  // } will be added up and reported together

    // "Normal" cases reported in a uniform way
    {QCflags::missing,       "missing values"},
    {QCflags::preQC,         "rejected by pre QC"},
    {QCflags::bounds,        "out of bounds"},
    {QCflags::domain,        "out of domain of use"},
    {QCflags::black,         "black-listed"},
    {QCflags::Hfailed,       "H(x) failed"},
    {QCflags::thinned,       "removed by thinning"},
    {QCflags::derivative,    "dy/dx out of valid range"},
    {QCflags::clw,           "removed by cloud liquid water check"},
    {QCflags::profile,       "removed by profile consistency check"},
    {QCflags::fguess,        "rejected by first-guess check"},
    {QCflags::diffref,       "rejected by difference check"},
    {QCflags::seaice,        "removed by sea ice check"},
    {QCflags::track,         "removed by track check"},
    {QCflags::buddy,         "removed by buddy check"},
    {QCflags::onedvar,       "removed by 1D Var check"},
    {QCflags::bayesianQC,    "removed by Bayesian background check"},
    {QCflags::modelobthresh, "removed by ModelOb threshold"}
  };
  const size_t numSpecialCases = 3;

  const size_t nlocs = obsdb_.nlocs();
  const size_t gnlocs = obsdb_.globalNumLocs();

  const oops::Variables &allSimulatedVars = obsdb_.obsvariables();

  for (size_t jvar = 0; jvar < allSimulatedVars.size(); ++jvar) {
    std::unique_ptr<ioda::Accumulator<std::vector<size_t>>> accumulator =
        obsdb_.distribution()->createAccumulator<size_t>(cases.size());

    for (size_t jobs = 0; jobs < nlocs; ++jobs) {
      const int actualFlag = (*flags_)[jvar][jobs];
      for (size_t jcase = 0; jcase < cases.size(); ++jcase)
        if (actualFlag == cases[jcase].first)
          accumulator->addTerm(jobs, jcase, 1);
    }
    const std::vector<std::size_t> counts = accumulator->computeResult();

    if (obsdb_.comm().rank() == 0) {
      const std::string info = "QC " + flags_->obstype() + " " + allSimulatedVars[jvar] + ": ";

      // Normal cases
      for (size_t i = numSpecialCases; i < counts.size(); ++i)
        if (counts[i] > 0)
          os << info << counts[i] << " " << cases[i].second << "." << std::endl;

      // Special cases: the GNSSRO check...
      const size_t nGNSSRO = counts[1] + counts[2];
      if (nGNSSRO > 0)
        os << info << nGNSSRO << " rejected by GNSSRO reality check." << std::endl;

      // ... the number of passed observations and the total number of observations.
      const size_t npass = counts[0];
      os << info << npass << " passed out of " << gnlocs << " observations." << std::endl;
    }

    const size_t numRecognizedFlags = std::accumulate(counts.begin(), counts.end(), 0);
    ASSERT(numRecognizedFlags == gnlocs);
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
