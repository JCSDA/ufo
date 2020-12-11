/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/QCmanager.h"

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// Presets for QC filters could be performed in a function outside of any class.
// We keep them as a filter for now. The main reason for this is to be able to use
// the factory for models not in UFO/IODA.

// -----------------------------------------------------------------------------

QCmanager::QCmanager(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                     std::shared_ptr<ioda::ObsDataVector<int> > qcflags,
                     std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : obsdb_(obsdb), config_(config), nogeovals_(), nodiags_(), flags_(qcflags),
    observed_(obsdb.obsvariables())
{
  oops::Log::trace() << "QCmanager::QCmanager starting " << config_ << std::endl;

  ASSERT(qcflags);
  ASSERT(obserr);

  ASSERT(flags_->nvars() == observed_.size());
  ASSERT(flags_->nlocs() == obsdb_.nlocs());
  ASSERT(obserr->nvars() == observed_.size());
  ASSERT(obserr->nlocs() == obsdb_.nlocs());

  const float rmiss = util::missingValue(rmiss);
  const int imiss = util::missingValue(imiss);

  const ioda::ObsDataVector<float> obs(obsdb, observed_, "ObsValue");

  for (size_t jv = 0; jv < observed_.size(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if ((*flags_)[jv][jobs] == imiss || obs[jv][jobs] == rmiss || (*obserr)[jv][jobs] == rmiss) {
        (*flags_)[jv][jobs] = QCflags::missing;
      }
    }
  }

  oops::Log::trace() << "QCmanager::QCmanager done" << std::endl;
}

// -----------------------------------------------------------------------------

void QCmanager::postFilter(const ioda::ObsVector & hofx, const ObsDiagnostics &) const {
  oops::Log::trace() << "QCmanager postFilter" << std::endl;

  const double missing = util::missingValue(missing);

  for (size_t jv = 0; jv < observed_.size(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      size_t iobs = observed_.size() * jobs + jv;
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
  for (size_t jj = 0; jj < observed_.size(); ++jj) {
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
    size_t idydx = 0;
    size_t iclw  = 0;
    size_t iprof = 0;
    size_t idiffref = 0;
    size_t iseaice  = 0;
    size_t itrack   = 0;
    size_t ibuddy   = 0;
    size_t iratioref = 0;
    size_t ionedvar  = 0;

    for (size_t jobs = 0; jobs < iobs; ++jobs) {
      if ((*flags_)[jj][jobs] == QCflags::pass)    ++ipass;
      if ((*flags_)[jj][jobs] == QCflags::missing) ++imiss;
      if ((*flags_)[jj][jobs] == QCflags::preQC)   ++ipreq;
      if ((*flags_)[jj][jobs] == QCflags::bounds)  ++ibnds;
      if ((*flags_)[jj][jobs] == QCflags::domain)  ++iwhit;
      if ((*flags_)[jj][jobs] == QCflags::black)   ++iblck;
      if ((*flags_)[jj][jobs] == QCflags::Hfailed) ++iherr;
      if ((*flags_)[jj][jobs] == QCflags::fguess)  ++ifgss;
      if ((*flags_)[jj][jobs] == QCflags::thinned) ++ithin;
      if ((*flags_)[jj][jobs] == QCflags::clw)     ++iclw;
      if ((*flags_)[jj][jobs] == QCflags::profile) ++iprof;
      if ((*flags_)[jj][jobs] == QCflags::diffref) ++idiffref;
      if ((*flags_)[jj][jobs] == QCflags::seaice)  ++iseaice;
      if ((*flags_)[jj][jobs] == 76 || (*flags_)[jj][jobs] == 77)  ++ignss;
      if ((*flags_)[jj][jobs] == QCflags::track)  ++itrack;
      if ((*flags_)[jj][jobs] == QCflags::buddy)  ++ibuddy;
      if ((*flags_)[jj][jobs] == QCflags::derivative) ++idydx;
      if ((*flags_)[jj][jobs] == QCflags::ratioref) ++iratioref;
      if ((*flags_)[jj][jobs] == QCflags::onedvar) ++ionedvar;
    }

    const ioda::Distribution & distribution = obsdb_.distribution();
    distribution.sum(iobs);
    distribution.sum(ipass);
    distribution.sum(imiss);
    distribution.sum(ipreq);
    distribution.sum(ibnds);
    distribution.sum(iwhit);
    distribution.sum(iblck);
    distribution.sum(iherr);
    distribution.sum(ifgss);
    distribution.sum(iclw);
    distribution.sum(iprof);
    distribution.sum(ignss);
    distribution.sum(ithin);
    distribution.sum(idiffref);
    distribution.sum(iseaice);
    distribution.sum(itrack);
    distribution.sum(ibuddy);
    distribution.sum(idydx);
    distribution.sum(iratioref);
    distribution.sum(ionedvar);

    if (obsdb_.comm().rank() == 0) {
      const std::string info = "QC " + flags_->obstype() + " " + observed_[jj] + ": ";
      if (imiss > 0) os << info << imiss << " missing values." << std::endl;
      if (ipreq > 0) os << info << ipreq << " rejected by pre QC." << std::endl;
      if (ibnds > 0) os << info << ibnds << " out of bounds." << std::endl;
      if (iwhit > 0) os << info << iwhit << " out of domain of use." << std::endl;
      if (iblck > 0) os << info << iblck << " black-listed." << std::endl;
      if (iherr > 0) os << info << iherr << " H(x) failed." << std::endl;
      if (ithin > 0) os << info << ithin << " removed by thinning." << std::endl;
      if (idydx > 0) os << info << idydx << " dy/dx out of valid range." << std::endl;
      if (iclw  > 0) os << info << iclw  << " removed by cloud liquid water check." << std::endl;
      if (iprof > 0) os << info << iprof  << " removed by profile consistency check." << std::endl;
      if (ifgss > 0) os << info << ifgss << " rejected by first-guess check." << std::endl;
      if (ignss > 0) os << info << ignss << " rejected by GNSSRO reality check." << std::endl;
      if (idiffref > 0) os << info << idiffref << " rejected by difference check." << std::endl;
      if (iseaice  > 0) os << info << iseaice  << " removed by sea ice check." << std::endl;
      if (itrack   > 0) os << info << itrack  << " removed by track check." << std::endl;
      if (ibuddy   > 0) os << info << ibuddy  << " removed by buddy check." << std::endl;
      if (iratioref > 0) os << info << iratioref << " rejected by ratio check." << std::endl;
      if (ionedvar  > 0) os << info << ionedvar  << " removed by 1D Var check." << std::endl;

      os << info << ipass << " passed out of " << iobs << " observations." << std::endl;
    }

    ASSERT(ipass + imiss + ipreq + ibnds + iwhit + iblck + iherr + ithin + iclw + iprof + ifgss + \
           ignss + idiffref + iseaice + itrack + ibuddy + idydx + iratioref + ionedvar == iobs);
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
