/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/atmsfcinterp/ObsAtmSfcInterp.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"


namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAtmSfcInterp> makerAtmSfcInterp_("2mTemp");
// -----------------------------------------------------------------------------

ObsAtmSfcInterp::ObsAtmSfcInterp(const ioda::ObsSpace & odb,
                       const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_(), varout_()
{
  // TODO(anyone): list the variables for GeoVaLs that are needed for the observation
  //       operator below in vv (e.g., vv{"temperature", "humidity"})
  const std::vector<std::string> vvin{""};
  varin_.reset(new oops::Variables(vvin));
  const std::vector<std::string> vvout{""};
  varout_.reset(new oops::Variables(vvout));
  const eckit::Configuration * configc = &config;
  ufo_atmsfcinterp_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsAtmSfcInterp created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmSfcInterp::~ObsAtmSfcInterp() {
  ufo_atmsfcinterp_delete_f90(keyOper_);
  oops::Log::trace() << "ObsAtmSfcInterp destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmSfcInterp::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                              const ObsBias & bias) const {
  ufo_atmsfcinterp_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran(),
                      bias.toFortran());
  oops::Log::trace() << "ObsAtmSfcInterp: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmSfcInterp::print(std::ostream & os) const {
  os << "ObsAtmSfcInterp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
