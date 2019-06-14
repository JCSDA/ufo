/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/marine/marinevertinterp/ObsMarineVertInterp.h"

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
static ObsOperatorMaker<ObsMarineVertInterp> makerMarineVertInterp_("InsituSalinity");
// -----------------------------------------------------------------------------

ObsMarineVertInterp::ObsMarineVertInterp(const ioda::ObsSpace & odb,
                                         const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_()
{
  int c_name_size = 200;
  char *buffin = new char[c_name_size];
  const eckit::Configuration * configc = &config;

  const oops::Variables & observed = odb.obsvariables();
  const eckit::Configuration * varconfig = &observed.toFortran();
  ufo_marinevertinterp_setup_f90(keyOper_, &configc, &varconfig, buffin, c_name_size);

  std::string vstr_in(buffin);
  std::vector<std::string> vvin;
  boost::split(vvin, vstr_in, boost::is_any_of("\t"));
  varin_.reset(new oops::Variables(vvin));

  oops::Log::trace() << "ObsMarineVertInterp created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsMarineVertInterp::~ObsMarineVertInterp() {
  ufo_marinevertinterp_delete_f90(keyOper_);
  oops::Log::trace() << "ObsMarineVertInterp destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterp::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec) const {
  ufo_marinevertinterp_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsMarineVertInterp: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterp::print(std::ostream & os) const {
  os << "ObsMarineVertInterp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
