/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/atmosphere/radiosonde/ObsRadiosonde.h"

#include <ostream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp> 

#include "oops/util/Logger.h"

#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadiosonde> makerRadiosonde_("Radiosonde");
// -----------------------------------------------------------------------------

ObsRadiosonde::ObsRadiosonde(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperRadiosonde_(0), odb_(odb), varin_(), varout_()
{
  const eckit::Configuration * configc = &config;
  ufo_radiosonde_setup_f90(keyOperRadiosonde_, &configc);
  
  char *buffin;
  char *buffout;
  ufo_radiosonde_getvars_f90(&configc,buffin,buffout);

  std::vector<std::string> vv,vout;
  std::string vstr_in(buffin), vstr_out(buffout);
  boost::split(vv, vstr_in, boost::is_any_of("\t"));
  boost::split(vout, vstr_out, boost::is_any_of("\t"));

  varin_.reset(new oops::Variables(vv));
  varout_.reset(new oops::Variables(vout));

  // Specify variables in C++ - should this be an option?
  //const std::vector<std::string> vv{"virtual_temperature", "atmosphere_ln_pressure_coordinate"};
  //varin_.reset(new oops::Variables(vv));
  //const std::vector<std::string> vout{"air_temperature"};
  //varout_.reset(new oops::Variables(vout));

  oops::Log::trace() << "ObsRadiosonde created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadiosonde::~ObsRadiosonde() {
  ufo_radiosonde_delete_f90(keyOperRadiosonde_);
  oops::Log::trace() << "ObsRadiosonde destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadiosonde::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                const ObsBias & bias) const {
  ufo_radiosonde_simobs_f90(keyOperRadiosonde_, gom.toFortran(), odb_,
                            ovec.size(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsRadiosonde::print(std::ostream & os) const {
  os << "ObsRadiosonde::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
