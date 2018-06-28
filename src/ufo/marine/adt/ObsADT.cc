/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsADT.h"

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "ioda/ObsSpace.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ioda/ObsVector.h"
#include "oops/util/ObjectCounter.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsADT> makerADT_("ADT");
// -----------------------------------------------------------------------------

ObsADT::ObsADT(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperADT_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_adt_setup_f90(keyOperADT_, &configc);
  const std::vector<std::string> vv{"sea_surface_height_above_geoid"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsADT created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsADT::~ObsADT() {
  ufo_adt_delete_f90(keyOperADT_);
  oops::Log::trace() << "ObsADT destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADT::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec, const ObsBias & bias) const {
  ufo_adt_eqv_f90(keyOperADT_, gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsADT::print(std::ostream & os) const {
  os << "ObsADT::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
