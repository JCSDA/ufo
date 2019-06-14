/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/PreQC.h"

#include <string>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/QCflags.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// Presets for QC filters could be performed in a function outside of any class.
// We keep them as a filter for now. The main reason for this is to be able to use
// the factory for models not in UFO/IODA.

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, PreQC>> mkPreQC_("PreQC");
// -----------------------------------------------------------------------------

PreQC::PreQC(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
             boost::shared_ptr<ioda::ObsDataVector<int> > qcflags,
             boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : nogeovals_()
{
  oops::Log::trace() << "PreQC::PreQC starting " << config << std::endl;
  const int missing = util::missingValue(missing);

// Basic arguments checks
  ASSERT(qcflags);
  ASSERT(obserr);

  const oops::Variables observed = obsdb.obsvariables();

  ASSERT(qcflags->nvars() == observed.size());
  ASSERT(qcflags->nlocs() == obsdb.nlocs());
  ASSERT(obserr->nvars() == observed.size());
  ASSERT(obserr->nlocs() == obsdb.nlocs());

// Read QC flags from pre-processing
  const std::string qcin = config.getString("inputQC", "PreQC");
  ioda::ObsDataVector<int> preqc(obsdb, observed, qcin);
  oops::Log::debug() << "PreQC::PreQC preqc: " << preqc;

// Get threshold and reject above threshold
  const int threshold = config.getInt("threshold", 0);

  for (size_t jv = 0; jv < observed.size(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb.nlocs(); ++jobs) {
      if (preqc[jv][jobs] == missing || preqc[jv][jobs] > threshold) {
        (*qcflags)[jv][jobs] = QCflags::preQC;
      }
    }
  }

  oops::Log::trace() << "PreQC::PreQC done" << std::endl;
}

// -----------------------------------------------------------------------------

void PreQC::print(std::ostream & os) const {
  os << "PreQC";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
