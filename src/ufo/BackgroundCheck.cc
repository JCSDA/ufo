/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/BackgroundCheck.h"

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, BackgroundCheck> >
  makerBgChk_("Background Check");
// -----------------------------------------------------------------------------

BackgroundCheck::BackgroundCheck(const ioda::ObsSpace & os,
                                 const eckit::Configuration & config) {
  oops::Log::debug() << "BackgroundCheck contructor starting" << std::endl;
  const eckit::Configuration * conf = &config;
  ufo_bgcheck_create_f90(key_, os, conf);
  oops::Log::debug() << "BackgroundCheck contructor key = " << key_ << std::endl;
}

// -----------------------------------------------------------------------------

BackgroundCheck::~BackgroundCheck() {
  oops::Log::debug() << "BackgroundCheck destructor key = " << key_ << std::endl;
  ufo_bgcheck_delete_f90(key_);
}

// -----------------------------------------------------------------------------

void BackgroundCheck::priorFilter(const GeoVaLs & gv) const {}

// -----------------------------------------------------------------------------

void BackgroundCheck::postFilter(const ioda::ObsVector & hofx) const {
  oops::Log::debug() << "BackgroundCheck postFilter" << std::endl;
  ufo_bgcheck_post_f90(key_, hofx.nvars(), hofx.nlocs(), hofx.toFortran(),
                       hofx.varnames().toFortranBetter());
}

// -----------------------------------------------------------------------------

void BackgroundCheck::print(std::ostream & os) const {
  os << "BackgroundCheck::print not yet implemented " << key_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
