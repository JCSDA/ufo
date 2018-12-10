/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/generic/identity/ObsIdentity.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"


namespace ufo {

// -----------------------------------------------------------------------------

ObsIdentity::ObsIdentity(const ioda::ObsSpace & odb,
                                     const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_(), varout_()
{
  // from config
  const std::vector<std::string> vvin{"ocean_upper_level_temperature"};
  varin_.reset(new oops::Variables(vvin));
  // from config
  const std::vector<std::string> vvout{"obs_sst"};
  varout_.reset(new oops::Variables(vvout));
  const eckit::Configuration * configc = &config;
  ufo_identity_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsIdentity created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsIdentity::~ObsIdentity() {
  ufo_identity_delete_f90(keyOper_);
  oops::Log::trace() << "ObsIdentity destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsIdentity::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                    const ObsBias & bias) const {
  ufo_identity_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran(),
                      bias.toFortran());
  oops::Log::trace() << "ObsIdentity: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsIdentity::print(std::ostream & os) const {
  os << "ObsIdentity::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
