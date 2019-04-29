/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/ObsPreQC.h"

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/QCflags.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// Presets for QC filters could be performed in a function outside of any class.
// We keep them as a filter for now. The main reason for this is to be able to use
// the factory for models not in UFO/IODA.

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, ObsPreQC>> mkPreQC_("PreQC");
// -----------------------------------------------------------------------------

ObsPreQC::ObsPreQC(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                   boost::shared_ptr<ioda::ObsDataVector<int> > qcflags,
                   boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : obsdb_(obsdb), config_(config), nogeovals_(), flags_(*qcflags)
{
  oops::Log::trace() << "ObsPreQC::ObsPreQC starting " << config_ << std::endl;

  ASSERT(qcflags);
  ASSERT(obserr);

  const oops::Variables observed(config.getStringVector("observed"));

  ASSERT(flags_.nvars() == observed.size());
  ASSERT(flags_.nlocs() == obsdb_.nlocs());
  ASSERT(obserr->nvars() == observed.size());
  ASSERT(obserr->nlocs() == obsdb_.nlocs());

  const float rmiss = util::missingValue(rmiss);
  const int imiss = util::missingValue(imiss);

  const ioda::ObsDataVector<float> obs(obsdb, observed, "ObsValue");

  for (size_t jv = 0; jv < observed.size(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (flags_[jv][jobs] == imiss || obs[jv][jobs] == rmiss || (*obserr)[jv][jobs] == rmiss) {
        flags_[jv][jobs] = QCflags::missing;
      } else {
        if (flags_[jv][jobs] > 3) {
          flags_[jv][jobs] = QCflags::preQC;
        } else {
          flags_[jv][jobs] = 0;
        }
      }
    }
  }
  oops::Log::trace() << "ObsPreQC::ObsPreQC done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsPreQC::postFilter(const ioda::ObsVector & hofx) const {
  oops::Log::trace() << "ObsPreQC postFilter" << std::endl;

  const oops::Variables observed(config_.getStringVector("observed"));
  const double missing = util::missingValue(missing);

  for (size_t jv = 0; jv < observed.size(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      size_t iobs = observed.size() * jobs + jv;
      if (flags_[jv][jobs] == 0 && hofx[iobs] == missing) {
        flags_[jv][jobs] = QCflags::Hfailed;
      }
    }
  }
}

// -----------------------------------------------------------------------------

ObsPreQC::~ObsPreQC() {
  const oops::Variables observed(config_.getStringVector("observed"));

  for (size_t jj = 0; jj < observed.size(); ++jj) {
    size_t iobs = obsdb_.nlocs();
    size_t ipass = 0;
    size_t imiss = 0;
    size_t ipreq = 0;
    size_t ibnds = 0;
    size_t iwhit = 0;
    size_t iblck = 0;
    size_t iherr = 0;
    size_t ifgss = 0;
    size_t ignss = 0;
    size_t ithin = 0;

    for (size_t jobs = 0; jobs < iobs; ++jobs) {
      if (flags_[jj][jobs] == QCflags::pass)    ++ipass;
      if (flags_[jj][jobs] == QCflags::missing) ++imiss;
      if (flags_[jj][jobs] == QCflags::preQC)   ++ipreq;
      if (flags_[jj][jobs] == QCflags::bounds)  ++ibnds;
      if (flags_[jj][jobs] == QCflags::domain)  ++iwhit;
      if (flags_[jj][jobs] == QCflags::black)   ++iblck;
      if (flags_[jj][jobs] == QCflags::Hfailed) ++iherr;
      if (flags_[jj][jobs] == QCflags::fguess)  ++ifgss;
      if (flags_[jj][jobs] == QCflags::thinned) ++ithin;
      if (flags_[jj][jobs] == 76 || flags_[jj][jobs] == 77)  ++ignss;
    }

    obsdb_.comm().allReduceInPlace(iobs, eckit::mpi::sum());
    obsdb_.comm().allReduceInPlace(ipass, eckit::mpi::sum());
    obsdb_.comm().allReduceInPlace(imiss, eckit::mpi::sum());
    obsdb_.comm().allReduceInPlace(ipreq, eckit::mpi::sum());
    obsdb_.comm().allReduceInPlace(ibnds, eckit::mpi::sum());
    obsdb_.comm().allReduceInPlace(iwhit, eckit::mpi::sum());
    obsdb_.comm().allReduceInPlace(iblck, eckit::mpi::sum());
    obsdb_.comm().allReduceInPlace(iherr, eckit::mpi::sum());
    obsdb_.comm().allReduceInPlace(ifgss, eckit::mpi::sum());
    obsdb_.comm().allReduceInPlace(ignss, eckit::mpi::sum());
    obsdb_.comm().allReduceInPlace(ithin, eckit::mpi::sum());

    if (obsdb_.comm().rank() == 0) {
      const std::string info = "QC " + flags_.obstype() + " " + observed[jj] + ": ";
      if (imiss > 0) oops::Log::info() << info << imiss << " missing values." << std::endl;
      if (ipreq > 0) oops::Log::info() << info << ipreq << " rejected by pre QC." << std::endl;
      if (ibnds > 0) oops::Log::info() << info << ibnds << " out of bounds." << std::endl;
      if (iwhit > 0) oops::Log::info() << info << iwhit << " out of domain of use." << std::endl;
      if (iblck > 0) oops::Log::info() << info << iblck << " black-listed." << std::endl;
      if (iherr > 0) oops::Log::info() << info << iherr << " H(x) failed." << std::endl;
      if (ithin > 0) oops::Log::info() << info << ithin << " removed by thinning." << std::endl;
      if (ifgss > 0) oops::Log::info() << info << ifgss << " rejected by first-guess check."
                                       << std::endl;
      if (ignss > 0) oops::Log::info() << info << ignss << " rejected by GNSSRO reality check."
                                       << std::endl;
      oops::Log::info() << info << ipass << " passed out of "
                        << iobs << " observations." << std::endl;
    }

    ASSERT(ipass + imiss + ipreq + ibnds + iwhit + iblck + iherr + ithin + ifgss + ignss == iobs);
  }
}

// -----------------------------------------------------------------------------

void ObsPreQC::print(std::ostream & os) const {
  os << "ObsPreQC";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
