/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/geos-aero/ObsGeosAod.h"

#include <ostream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"


namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGeosAod> makerGeosAod_("GeosAod");
// -----------------------------------------------------------------------------

ObsGeosAod::ObsGeosAod(const ioda::ObsSpace & odb,
                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_(), varout_()
{
  int c_name_size = 800;
  char *buffin = new char[c_name_size];
  char *buffout = new char[c_name_size];
  const eckit::Configuration * configc = &config;

  ufo_geosaod_setup_f90(keyOper_, &configc, buffin, buffout, c_name_size);

  std::string vstr_in(buffin), vstr_out(buffout);
  std::vector<std::string> vvin;
  std::vector<std::string> vvout;
  boost::split(vvin, vstr_in, boost::is_any_of("\t"));
  boost::split(vvout, vstr_out, boost::is_any_of("\t"));
  varin_.reset(new oops::Variables(vvin));
  varout_.reset(new oops::Variables(vvout));

  oops::Log::trace() << "ObsGeosAod created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGeosAod::~ObsGeosAod() {
  ufo_geosaod_delete_f90(keyOper_);
  oops::Log::trace() << "ObsGeosAod destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGeosAod::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                              const ObsBias & bias) const {
  ufo_geosaod_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran(),
                      bias.toFortran());
  oops::Log::trace() << "ObsGeosAod: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGeosAod::print(std::ostream & os) const {
  os << "ObsGeosAod::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
