/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsSeaSurfaceTemp.h"

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
static ObsOperatorMaker<ObsSeaSurfaceTemp> makerSeaSurfaceTemp_("SeaSurfaceTemp");
// -----------------------------------------------------------------------------

ObsSeaSurfaceTemp::ObsSeaSurfaceTemp(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperSeaSurfaceTemp_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_seasurfacetemp_setup_f90(keyOperSeaSurfaceTemp_, &configc);
  const std::vector<std::string> vv{"ocean_upper_level_temperature"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsSeaSurfaceTemp created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaSurfaceTemp::~ObsSeaSurfaceTemp() {
  ufo_seasurfacetemp_delete_f90(keyOperSeaSurfaceTemp_);
  oops::Log::trace() << "ObsSeaSurfaceTemp destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaSurfaceTemp::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec, const ObsBias & bias) const {
  ufo_seasurfacetemp_eqv_f90(keyOperSeaSurfaceTemp_, gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsSeaSurfaceTemp::print(std::ostream & os) const {
  os << "ObsSeaSurfaceTemp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
