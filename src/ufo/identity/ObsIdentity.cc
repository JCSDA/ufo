/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/identity/ObsIdentity.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsIdentity> makerIdentity_("Identity");
// -----------------------------------------------------------------------------

ObsIdentity::ObsIdentity(const ioda::ObsSpace & odb,
                         const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperObsIdentity_(0), odb_(odb), varin_()
{
  const eckit::Configuration * configc = &config;
  ufo_identity_setup_f90(keyOperObsIdentity_, &configc, odb.obsvariables(), varin_);

  oops::Log::trace() << "ObsIdentity created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsIdentity::~ObsIdentity() {
  ufo_identity_delete_f90(keyOperObsIdentity_);
  oops::Log::trace() << "ObsIdentity destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsIdentity::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                              ObsDiagnostics &) const {
  ufo_identity_simobs_f90(keyOperObsIdentity_, gom.toFortran(), odb_,
                          ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsIdentity: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsIdentity::print(std::ostream & os) const {
  os << "ObsIdentity::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
