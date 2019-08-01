/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/geos_aero/ObsGeosAod.h"

#include <ostream>
#include <set>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "oops/util/Logger.h"
#include "oops/util/IntSetParser.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGeosAod> makerGeosAod_("GeosAod");
// -----------------------------------------------------------------------------

ObsGeosAod::ObsGeosAod(const ioda::ObsSpace & odb,
                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_()
{
  const eckit::Configuration * configc = &config;

  const oops::Variables & observed = odb.obsvariables();
  const eckit::Configuration * varconfig = &observed.toFortran();
//  ufo_geosaod_setup_f90(keyOper_, &configc, &varconfig, buffin, c_name_size);

//  std::string vstr_in(buffin);
//  std::vector<std::string> vvin;
//  boost::split(vvin, vstr_in, boost::is_any_of("\t"));
//  varin_.reset(new oops::Variables(vvin));

  ufo_geosaod_setup_f90(keyOper_, &configc, &varconfig, varin_);

  oops::Log::trace() << "ObsGeosAod created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGeosAod::~ObsGeosAod() {
  ufo_geosaod_delete_f90(keyOper_);
  oops::Log::trace() << "ObsGeosAod destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGeosAod::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec) const {
  ufo_geosaod_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsGeosAod: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGeosAod::print(std::ostream & os) const {
  os << "ObsGeosAod::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
