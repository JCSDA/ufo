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
#include <boost/algorithm/string.hpp>

#include "oops/util/Logger.h"

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"


namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAtmSfcInterp> maker2mTemp_("2mTemp");
// -----------------------------------------------------------------------------

ObsAtmSfcInterp::ObsAtmSfcInterp(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperAtmSfcInterp_(0),
    odb_(odb), varin_(), varout_()
{
  int c_name_size = 200;
  char *buffin = new char[c_name_size];
  char *buffout = new char[c_name_size];
  const eckit::Configuration * configc = &config;

  ufo_atmsfcinterp_setup_f90(keyOperAtmSfcInterp_, &configc, buffin, buffout, c_name_size);

  std::string vstr_in(buffin), vstr_out(buffout);
  std::vector<std::string> vvin;
  std::vector<std::string> vvout;
  boost::split(vvin, vstr_in, boost::is_any_of("\t"));
  boost::split(vvout, vstr_out, boost::is_any_of("\t"));
  varin_.reset(new oops::Variables(vvin));
  varout_.reset(new oops::Variables(vvout));

  oops::Log::trace() << "ObsAtmSfcInterp created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmSfcInterp::~ObsAtmSfcInterp() {
  ufo_atmsfcinterp_delete_f90(keyOperAtmSfcInterp_);
  oops::Log::trace() << "ObsAtmSfcInterp destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmSfcInterp::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                              const ObsBias & bias) const {
  ufo_atmsfcinterp_simobs_f90(keyOperAtmSfcInterp_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran(),
                      bias.toFortran());
  oops::Log::trace() << "ObsAtmSfcInterp: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmSfcInterp::print(std::ostream & os) const {
  os << "ObsAtmSfcInterp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
