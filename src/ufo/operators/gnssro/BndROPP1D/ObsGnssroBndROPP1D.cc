/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/gnssro/BndROPP1D/ObsGnssroBndROPP1D.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGnssroBndROPP1D> makerGnssroBndROPP1D_("GnssroBndROPP1D");
// -----------------------------------------------------------------------------

ObsGnssroBndROPP1D::ObsGnssroBndROPP1D(const ioda::ObsSpace & odb,
                                       const Parameters_ & params)
  : ObsOperatorBase(odb), keyOperGnssroBndROPP1D_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vv{"air_temperature", "specific_humidity", "air_pressure",
                                    "geopotential_height", "surface_altitude"};
  varin_.reset(new oops::Variables(vv));

  ufo_gnssro_bndropp1d_setup_f90(keyOperGnssroBndROPP1D_,
                                 params.options.value().toConfiguration());
  oops::Log::trace() << "ObsGnssroBndROPP1D created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndROPP1D::~ObsGnssroBndROPP1D() {
  ufo_gnssro_bndropp1d_delete_f90(keyOperGnssroBndROPP1D_);
  oops::Log::trace() << "ObsGnssroBndROPP1D destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP1D::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                     ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  ufo_gnssro_bndropp1d_simobs_f90(keyOperGnssroBndROPP1D_, gom.toFortran(), odb_,
                                  ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP1D::print(std::ostream & os) const {
  os << "ObsGnssroBndROPP1D::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
