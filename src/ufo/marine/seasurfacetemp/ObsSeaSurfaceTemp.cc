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

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"


namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsSeaSurfaceTemp> makerSeaSurfaceTemp_("SeaSurfaceTemp");
// -----------------------------------------------------------------------------

ObsSeaSurfaceTemp::ObsSeaSurfaceTemp(const ioda::ObsSpace & odb,
                                     const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_(), varout_()
{
  const std::vector<std::string> vvin{"ocean_upper_level_temperature"};
  varin_.reset(new oops::Variables(vvin));
  const std::vector<std::string> vvout{"obs_sst"};
  varout_.reset(new oops::Variables(vvout));
  const eckit::Configuration * configc = &config;
  ufo_seasurfacetemp_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsSeaSurfaceTemp created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaSurfaceTemp::~ObsSeaSurfaceTemp() {
  ufo_seasurfacetemp_delete_f90(keyOper_);
  oops::Log::trace() << "ObsSeaSurfaceTemp destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaSurfaceTemp::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                    const ObsBias & bias) const {
  ufo_seasurfacetemp_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran(),
                      bias.toFortran());
  oops::Log::trace() << "ObsSeaSurfaceTemp: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaSurfaceTemp::print(std::ostream & os) const {
  os << "ObsSeaSurfaceTemp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
