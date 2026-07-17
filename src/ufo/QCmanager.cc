/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/QCmanager.h"

#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "eckit/utils/StringTools.h"
#include "ioda/distribution/Accumulator.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/ObsVariables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

namespace {

/// At each location, set the QC flag to QCflags::missing if the current QC flag is invalid
/// or if the ObsValue is missing.
void updateQCFlags(const std::vector<float> *obsValues, std::vector<int>& qcflags) {
  const float rmiss = util::missingValue<float>();
  const int imiss = util::missingValue<int>();
  for (size_t jobs = 0; jobs < qcflags.size(); ++jobs) {
    if (qcflags[jobs] == imiss || !obsValues || (*obsValues)[jobs] == rmiss) {
      qcflags[jobs] = QCflags::missing;
    }
  }
}

}  // namespace

// -----------------------------------------------------------------------------

QCmanager::QCmanager(ioda::ObsSpace & obsdb,
                     ioda::ObsDataVector<int> & flags,
                     ioda::ObsDataVector<float> & obserr)
  : obsdb_(obsdb), flags_(flags), obserr_(obserr)
{
  oops::Log::trace() << "QCmanager constructor" << std::endl;
}

// -----------------------------------------------------------------------------

QCmanager::~QCmanager() {
  oops::Log::trace() << "QCmanager destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void QCmanager::preSetQc() {
  oops::Log::trace() << "QCmanager preSetQc start" << std::endl;
  const oops::ObsVariables &allObservedVars = obsdb_.obsvariables();
  const oops::ObsVariables &initialObservedVars = obsdb_.initial_obsvariables();
  const oops::ObsVariables &derivedObservedVars = obsdb_.derived_obsvariables();

  ASSERT(allObservedVars.size() == initialObservedVars.size() + derivedObservedVars.size());
  ASSERT(flags_.nvars() == allObservedVars.size());
  ASSERT(flags_.nlocs() == obsdb_.nlocs());

  const ioda::ObsDataVector<float> obs(obsdb_, initialObservedVars, "ObsValue");

  // Iterate over initial observed variables
  for (size_t jv = 0; jv < initialObservedVars.size(); ++jv) {
    const ioda::ObsDataRow<float> &currentObsValues = obs[jv];
    ioda::ObsDataRow<int> &currentQCFlags = flags_[obs.varnames()[jv]];
    updateQCFlags(&currentObsValues, currentQCFlags);
  }

  // Iterate over derived variables and if they don't exist yet, set their QC flags to
  // 'missing'.
  for (size_t jv = 0; jv < derivedObservedVars.size(); ++jv) {
    ioda::ObsDataRow<int> &currentQCFlags = flags_[derivedObservedVars[jv]];
    if (!obsdb_.has("ObsValue", derivedObservedVars[jv])) {
      updateQCFlags(nullptr, currentQCFlags);
    } else {
      std::vector<float> currentObsValues(obsdb_.nlocs());
      obsdb_.get_db("ObsValue", derivedObservedVars[jv], currentObsValues);
      updateQCFlags(&currentObsValues, currentQCFlags);
    }
  }
  oops::Log::trace() << "QCmanager preSetQc complete" << std::endl;
}

// -----------------------------------------------------------------------------

void QCmanager::postSetQc(const ioda::ObsVector & hofx) {
  oops::Log::trace() << "QCmanager postSetQc start" << std::endl;

  const double missing = util::missingValue<double>();
  const oops::ObsVariables &allObservedVars = obsdb_.assimvariables();

  for (size_t jv = 0; jv < allObservedVars.size(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      const size_t iobs = allObservedVars.size() * jobs + jv;
      if (flags_[allObservedVars[jv]][jobs] == 0 && hofx[iobs] == missing) {
        flags_[allObservedVars[jv]][jobs] = QCflags::Hfailed;
      }
    }
  }
  oops::Log::trace() << "QCmanager postSetQc complete" << std::endl;
}

// -----------------------------------------------------------------------------

void QCmanager::finalSetQc() {
  oops::Log::trace() << "QCmanager finalSetQc start" << std::endl;

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
      if (flags_[jv][jobs] == QCflags::pass && obserr_[jv][jobs] == missing) {
        flags_[jv][jobs] = QCflags::missing;
      }
    }
  }

  // Set the QC flag to observations that are processed but not to
  // be assimilated.
  for (size_t jv = 0; jv < obsdb_.obsvariables().size(); ++jv) {
    if (!obsdb_.assimvariables().has(obsdb_.obsvariables()[jv])) {
      for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
        flags_[jv][jobs] = QCflags::processed;
      }
    }
  }

  oops::Log::trace() << "QCmanager finalSetQc complete" << std::endl;
  // Print QC statistics
  oops::Log::info() << *this;
}

// -----------------------------------------------------------------------------

void QCmanager::print(std::ostream & os) const {
  const std::vector<std::pair<int, const char*>> cases{
    // Special cases reported using dedicated code
    {QCflags::pass, nullptr},
    {76, nullptr},  // } The numbers of observations with these two flags
    {77, nullptr},  // } will be added up and reported together

    // "Normal" cases reported in a uniform way
    {QCflags::passive,       "passive observations"},
    {QCflags::missing,       "missing values"},
    {QCflags::preQC,         "rejected by pre QC"},
    {QCflags::bounds,        "out of bounds"},
    {QCflags::domain,        "out of domain of use"},
    {QCflags::black,         "black-listed"},
    {QCflags::Hfailed,       "H(x) failed"},
    {QCflags::thinned,       "removed by thinning"},
    {QCflags::derivative,    "dy/dx out of valid range"},
    {QCflags::clw,           "removed by cloud liquid water check"},
    {QCflags::profile,       "removed by conventional profile processing"},
    {QCflags::fguess,        "rejected by first-guess check"},
    {QCflags::diffref,       "rejected by difference check"},
    {QCflags::seaice,        "removed by sea ice check"},
    {QCflags::track,         "removed by track check"},
    {QCflags::buddy,         "removed by buddy check"},
    {QCflags::onedvar,       "removed by 1D Var check"},
    {QCflags::bayesianQC,    "removed by Bayesian background check"},
    {QCflags::modelobthresh, "removed by ModelOb threshold"},
    {QCflags::history,       "removed by history check"},
    {QCflags::processed,     "rejected as processed but not assimilated"},
    {QCflags::superrefraction, "rejected by GNSSRO super refraction QC"},
    {QCflags::superob,       "rejected by superobbing"},
    {QCflags::percentile,    "rejected by percentile filter"}
  };
  const size_t numSpecialCases = 3;

  const size_t nlocs = obsdb_.nlocs();
  const size_t gnlocs = obsdb_.globalNumLocs();

  const oops::ObsVariables &allObservedVars = obsdb_.obsvariables();

  const std::vector<std::string> obsSpaceVarList = obsdb_.listVariables();

  for (size_t jvar = 0; jvar < allObservedVars.size(); ++jvar) {
    const std::string varName = allObservedVars[jvar];
    std::unique_ptr<ioda::Accumulator<std::vector<size_t>>> accumulator =
        obsdb_.distribution()->createAccumulator<size_t>(cases.size());
    for (size_t jobs = 0; jobs < nlocs; ++jobs) {
      const int actualFlag = flags_[jvar][jobs];
      for (size_t jcase = 0; jcase < cases.size(); ++jcase)
        if (actualFlag == cases[jcase].first)
          accumulator->addTerm(jobs, jcase, 1);
    }
    const std::vector<std::size_t> counts = accumulator->computeResult();

    // Get list of diagnostic flags for this variable.
    // The first slash is guaranteed to appear at position 15 in the full variable name.
    const std::size_t slashFirst = 15;
    std::vector<std::string> diagFlagNames;
    for (const std::string & obsSpaceVarName : obsSpaceVarList) {
      const std::size_t slashLast = obsSpaceVarName.find_last_of("/");
      if (slashLast == std::string::npos) {
        continue;
      }
      // Check whether the full variable name starts with DiagnosticFlags and ends
      // with the name of the ObsSpace variable (potentially including a channel number).
      if (eckit::StringTools::startsWith(obsSpaceVarName, "DiagnosticFlags") &&
          varName.find(obsSpaceVarName.substr(slashLast + 1)) != std::string::npos) {
        // There must be at least two slashes in the full variable name.
        if (slashFirst == slashLast) {
          continue;
        }
        const std::string diagFlagName =
          obsSpaceVarName.substr(slashFirst + 1, slashLast - slashFirst - 1);
        // Ensure this flag is in the ObsSpace.
        if (obsdb_.has("DiagnosticFlags/" + diagFlagName, varName)) {
          diagFlagNames.push_back(diagFlagName);
        }
      }
    }

    // Create an accumulator to count the number of locations at which
    // the diagnostic flags have been set to true.
    std::unique_ptr<ioda::Accumulator<std::vector<size_t>>> accumulatorDiags =
        obsdb_.distribution()->createAccumulator<size_t>(diagFlagNames.size());

    // Retrieve all of the relevant diagnostic flags from the ObsSpace.
    std::vector<std::vector<bool>> diagFlagValues;
    for (const std::string & diagFlagName : diagFlagNames) {
      std::vector<bool> diagFlagBool(nlocs);
      obsdb_.get_db("DiagnosticFlags/" + diagFlagName, varName, diagFlagBool);
      diagFlagValues.push_back(diagFlagBool);
    }

    // Increment the accumulator for each diagnostic flag.
    for (size_t jdiag = 0; jdiag < diagFlagValues.size(); ++jdiag) {
      for (size_t jobs = 0; jobs < nlocs; ++jobs) {
        if (diagFlagValues[jdiag][jobs]) {
          accumulatorDiags->addTerm(jobs, jdiag, 1);
        }
      }
    }

    // Produce totals for each diagnostic flag.
    const std::vector<std::size_t> countDiags = accumulatorDiags->computeResult();

    if (obsdb_.comm().rank() == 0) {
      const std::string info = "QC " + flags_.obstype() + " " + varName + ": ";

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

      // Print diagnostic flags.
      for (size_t jdiag = 0; jdiag < diagFlagValues.size(); ++jdiag) {
        os << info << countDiags[jdiag] << " assigned the "
           << diagFlagNames[jdiag] << " diagnostic flag." << std::endl;
      }
    }
    const size_t numRecognizedFlags = std::accumulate(counts.begin(), counts.end(), 0);
    ASSERT(numRecognizedFlags == gnlocs);
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
