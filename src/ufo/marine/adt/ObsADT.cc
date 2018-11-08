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

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsADT> makerADT_("ADT");
// -----------------------------------------------------------------------------

ObsADT::ObsADT(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperADT_(0), odb_(odb), varin_(), varout_()
{
  const eckit::Configuration * configc = &config;
  ufo_adt_setup_f90(keyOperADT_, &configc);

  const std::vector<std::string> vv{"sea_surface_height_above_geoid"};
  varin_.reset(new oops::Variables(vv));

  const std::vector<std::string> vout{"zz"};
  varout_.reset(new oops::Variables(vout));

  oops::Log::trace() << "ObsADT created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsADT::~ObsADT() {
  ufo_adt_delete_f90(keyOperADT_);
  oops::Log::trace() << "ObsADT destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADT::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec, const ObsBias & bias) const {
  ufo_adt_simobs_f90(keyOperADT_, gom.toFortran(), odb_,
                     ovec.size(), ovec.toFortran(),
                     bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsADT::print(std::ostream & os) const {
  os << "ObsADT::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
