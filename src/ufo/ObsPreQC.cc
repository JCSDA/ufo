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

ObsPreQC::ObsPreQC(ioda::ObsSpace & obsdb, const eckit::Configuration & config)
  : obsdb_(obsdb), config_(config), nogeovals_()
{
  const float rmiss = util::missingValue(rmiss);
  const int imiss = util::missingValue(imiss);
  const oops::Variables vars(config.getStringVector("observed"));

  ioda::ObsDataVector<float> obs(obsdb, vars, "ObsValue");
  ioda::ObsDataVector<float> err(obsdb, vars, "ObsError");
  ioda::ObsDataVector<int> flags(obsdb, vars);
  ioda::ObsDataVector<int> preqc(obsdb, vars);
  if (obsdb.has(vars[0], "PreQC")) flags.read("PreQC");

  for (size_t jv = 0; jv < vars.size(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (flags[jv][jobs] == imiss || obs[jv][jobs] == rmiss || err[jv][jobs] == rmiss) {
        preqc[jv][jobs] = QCflags::missing;
      } else {
        if (flags[jv][jobs] != 0) preqc[jv][jobs] = QCflags::preQC;
      }
    }
    preqc.save(config.getString("QCname"));
  }
}

// -----------------------------------------------------------------------------

void ObsPreQC::postFilter(const ioda::ObsVector & hofx) const {
  oops::Log::trace() << "ObsPreQC postFilter" << std::endl;

  const oops::Variables vars(config_.getStringVector("observed"));
  const std::string qcgrp = config_.getString("QCname");
  const double missing = util::missingValue(missing);

  ioda::ObsDataVector<int> flags(obsdb_, vars, qcgrp);

  for (size_t jv = 0; jv < vars.size(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      size_t iobs = vars.size() * jobs + jv;
      if (flags[jv][jobs] == 0 && hofx[iobs] == missing) {
        flags[jv][jobs] = QCflags::Hfailed;
      }
    }
  }
  flags.save(qcgrp);
}

// -----------------------------------------------------------------------------

ObsPreQC::~ObsPreQC() {
  const std::string qcgrp = config_.getString("QCname");
  const oops::Variables vars(config_.getStringVector("observed"));
  ioda::ObsDataVector<int> flags(obsdb_, vars, qcgrp);

  for (size_t jj = 0; jj < vars.size(); ++jj) {
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
      if (flags[jj][jobs] == QCflags::pass)    ++ipass;
      if (flags[jj][jobs] == QCflags::missing) ++imiss;
      if (flags[jj][jobs] == QCflags::preQC)   ++ipreq;
      if (flags[jj][jobs] == QCflags::bounds)  ++ibnds;
      if (flags[jj][jobs] == QCflags::domain)  ++iwhit;
      if (flags[jj][jobs] == QCflags::black)   ++iblck;
      if (flags[jj][jobs] == QCflags::Hfailed) ++iherr;
      if (flags[jj][jobs] == QCflags::fguess)  ++ifgss;
      if (flags[jj][jobs] == QCflags::thinned) ++ithin;
      if (flags[jj][jobs] == 76 || flags[jj][jobs] == 77)  ++ignss;
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
      const std::string info = "QC " + flags.obstype() + " " + vars[jj] + ": ";
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
