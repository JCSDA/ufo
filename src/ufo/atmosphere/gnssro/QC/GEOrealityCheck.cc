/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "GEOrealityCheck.h"
#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/interface/ObsFilter.h"
#include "oops/base/ObsFilterBase.h"

#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, GEOrealityCheck> >
     makerROgeorealityChk_("GEOreality Check");
// -----------------------------------------------------------------------------

GEOrealityCheck::GEOrealityCheck(const ioda::ObsSpace & os,
                                 const eckit::Configuration & config) {
  oops::Log::debug() << "GEOrealityCheck contructor starting" << std::endl;
  const eckit::Configuration * conf = &config;
  ufo_georealitycheck_create_f90(key_, os, conf);
  oops::Log::debug() << "GEOrealityCheck contructor key = " << key_ << std::endl;
}

// -----------------------------------------------------------------------------

GEOrealityCheck::~GEOrealityCheck() {
  oops::Log::debug() << "GEOrealityCheck destructor key = " << key_ << std::endl;
  ufo_georealitycheck_delete_f90(key_);
}

// -----------------------------------------------------------------------------

void GEOrealityCheck::priorFilter(const GeoVaLs & gv) const {
  oops::Log::debug() << "GEOrealityCheck priorFilter" << std::endl;
  ufo_georealitycheck_prior_f90(key_, gv.toFortran());
}

// -----------------------------------------------------------------------------

void GEOrealityCheck::postFilter(const ioda::ObsVector & hofxb) const {
  oops::Log::debug() << "GEOrealityCheck postFilter" << std::endl;
  ufo_georealitycheck_post_f90(key_, hofxb.size(), hofxb.toFortran());
}

// -----------------------------------------------------------------------------

void GEOrealityCheck::print(std::ostream & os) const {
  os << "GEOrealityCheck::print not yet implemented " << key_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
