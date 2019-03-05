/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */
#include "ufo/gnssro/QC/ROgeorealityCheck.h"

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, ROgeorealityCheck> >
     makerROgeorealityChk_("ROgeoreality Check");
// -----------------------------------------------------------------------------

ROgeorealityCheck::ROgeorealityCheck(const ioda::ObsSpace & os,
                                 const eckit::Configuration & config) {
  oops::Log::debug() << "ROgeorealityCheck contructor starting" << std::endl;
  const eckit::Configuration * conf = &config;
  ufo_rogeorealitycheck_create_f90(key_, os, conf);
  oops::Log::debug() << "ROgeorealityCheck contructor key = " << key_ << std::endl;
}

// -----------------------------------------------------------------------------

ROgeorealityCheck::~ROgeorealityCheck() {
  oops::Log::debug() << "ROgeorealityCheck destructor key = " << key_ << std::endl;
  ufo_rogeorealitycheck_delete_f90(key_);
}

// -----------------------------------------------------------------------------

void ROgeorealityCheck::priorFilter(const GeoVaLs & gv) const {
  oops::Log::debug() << "ROgeorealityCheck priorFilter" << std::endl;
  ufo_rogeorealitycheck_prior_f90(key_, gv.toFortran());
}

// -----------------------------------------------------------------------------

void ROgeorealityCheck::postFilter(const ioda::ObsVector & hofxb) const {
  oops::Log::debug() << "ROgeorealityCheck postFilter" << std::endl;
  ufo_rogeorealitycheck_post_f90(key_, hofxb.size(), hofxb.toFortran());
}

// -----------------------------------------------------------------------------

void ROgeorealityCheck::print(std::ostream & os) const {
  os << "ROgeorealityCheck::print not yet implemented " << key_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
