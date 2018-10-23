/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/ObsCheck.h"

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"

#include "ufo/Fortran.h"
#include "ufo/FortranObsCheck.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsCheck::ObsCheck(const ioda::ObsSpace & obsdb, const oops::Variables & var,
                   const util::DateTime & t1, const util::DateTime & t2) {
  oops::Log::trace() << "ObsCheck contructor starting " << t1 << " " << t2 << std::endl;
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
//  ufo_obsdb_getobscheck_f90(obsdb, var.toFortran(), &p1, &p2, keyObsCheck_);
  oops::Log::trace() << "ObsCheck contructor key = " << keyObsCheck_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsCheck::ObsCheck(const eckit::Configuration & config) {
  oops::Log::trace() << "ObsCheck contructor config starting" << std::endl;
  const eckit::Configuration * conf = &config;
  ufo_obscheck_setup_f90(keyObsCheck_, &conf);
  oops::Log::trace() << "ObsCheck contructor config key = " << keyObsCheck_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsCheck::ObsCheck(const ioda::ObsSpace & os) {
  oops::Log::trace() << "ObsCheck ObsSpace starting" << std::endl;
  oops::Log::trace() << "ObsCheck ObsSpace end " << std::endl;
}

// -----------------------------------------------------------------------------

ObsCheck::~ObsCheck() {
  ufo_obscheck_delete_f90(keyObsCheck_);
  oops::Log::trace() << "ObsCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsCheck::print(std::ostream & os) const {
  os << "ObsCheck::print not implemented";
}

// -----------------------------------------------------------------------------

void ObsCheck::postFilter(const GeoVaLs & gv, const ioda::ObsVector & ov,
                          const ioda::ObsSpace & os) const {
  oops::Log::trace() << "ObsCheck postFilter starting" << std::endl;
  ufo_postFilter_f90(gv.toFortran(), ov.toFortran(), os);
  oops::Log::trace() << "ObsCheck postFilter end" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsCheck::priorFilter(const ioda::ObsSpace & os) const {
  oops::Log::trace() << "ObsCheck priorFilter starting" << std::endl;
  ufo_priorFilter_f90(os);
  oops::Log::trace() << "ObsCheck priorFilter end" << std::endl;
}

// -----------------------------------------------------------------------------
}  // namespace ufo
