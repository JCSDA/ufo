/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/identity/ObsIdentity.h"

#include <ostream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "oops/util/Logger.h"

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"


namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsIdentity> makerSST_("SeaSurfaceTemp");
static ObsOperatorMaker<ObsIdentity> makerSSS_("SeaSurfaceSalinity");
// -----------------------------------------------------------------------------

ObsIdentity::ObsIdentity(const ioda::ObsSpace & odb,
                         const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperObsIdentity_(0), odb_(odb), varin_()
{
  int c_name_size = 200;
  char *buffin = new char[c_name_size];
  const eckit::Configuration * configc = &config;

  const oops::Variables & vars = odb.obsvariables();
  const eckit::Configuration * varconfig = &vars.toFortran();

  ufo_identity_setup_f90(keyOperObsIdentity_, &configc, &varconfig, buffin, c_name_size);

  std::string vstr_in(buffin);
  std::vector<std::string> vvin;
  boost::split(vvin, vstr_in, boost::is_any_of("\t"));
  varin_.reset(new oops::Variables(vvin));

  oops::Log::trace() << "ObsIdentity created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsIdentity::~ObsIdentity() {
  ufo_identity_delete_f90(keyOperObsIdentity_);
  oops::Log::trace() << "ObsIdentity destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsIdentity::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec) const {
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
