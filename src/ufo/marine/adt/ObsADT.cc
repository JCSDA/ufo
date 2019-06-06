/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/marine/adt/ObsADT.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsADT> makerADT_("ADT");
// -----------------------------------------------------------------------------

ObsADT::ObsADT(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vvin{"sea_surface_height_above_geoid"};
  varin_.reset(new oops::Variables(vvin));
  const eckit::Configuration * configc = &config;
  ufo_adt_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsADT created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsADT::~ObsADT() {
  ufo_adt_delete_f90(keyOper_);
  oops::Log::trace() << "ObsADT destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADT::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec) const {
  ufo_adt_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsADT: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADT::print(std::ostream & os) const {
  os << "ObsADT::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
