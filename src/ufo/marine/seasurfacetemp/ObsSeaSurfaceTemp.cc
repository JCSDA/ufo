/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/marine/seasurfacetemp/ObsSeaSurfaceTemp.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsSeaSurfaceTemp> makerSeaSurfaceTemp_("SeaSurfaceTemp");
// -----------------------------------------------------------------------------

ObsSeaSurfaceTemp::ObsSeaSurfaceTemp(const ioda::ObsSpace & odb,
                                     const eckit::Configuration & config)
  : keyOperSeaSurfaceTemp_(0), odb_(odb), varin_(), varout_()
{
  const eckit::Configuration * configc = &config;
  ufo_seasurfacetemp_setup_f90(keyOperSeaSurfaceTemp_, &configc);

  const std::vector<std::string> vv{"ocean_upper_level_temperature"};
  varin_.reset(new oops::Variables(vv));

  const std::vector<std::string> vout{"zz"};
  varout_.reset(new oops::Variables(vout));

  oops::Log::trace() << "ObsSeaSurfaceTemp created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaSurfaceTemp::~ObsSeaSurfaceTemp() {
  ufo_seasurfacetemp_delete_f90(keyOperSeaSurfaceTemp_);
  oops::Log::trace() << "ObsSeaSurfaceTemp destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaSurfaceTemp::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                    const ObsBias & bias) const {
  ufo_seasurfacetemp_simobs_f90(keyOperSeaSurfaceTemp_, gom.toFortran(), odb_,
                             ovec.size(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

void ObsSeaSurfaceTemp::print(std::ostream & os) const {
  os << "ObsSeaSurfaceTemp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
