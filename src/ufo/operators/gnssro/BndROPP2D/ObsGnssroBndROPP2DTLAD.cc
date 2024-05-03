/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/gnssro/BndROPP2D/ObsGnssroBndROPP2DTLAD.h"

#include <ostream>
#include <string>
#include <vector>


#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsGnssroBndROPP2DTLAD> makerGnssroBndROPP2DTL_("GnssroBndROPP2D");
// -----------------------------------------------------------------------------

ObsGnssroBndROPP2DTLAD::ObsGnssroBndROPP2DTLAD(const ioda::ObsSpace & odb,
                                               const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOperGnssroBndROPP2D_(0), varin_()
{
  ufo_gnssro_bndropp2d_tlad_setup_f90(keyOperGnssroBndROPP2D_,
                                      params.options.value().toConfiguration());
  const std::vector<std::string> vv{"air_temperature", "specific_humidity", "air_pressure"};

  varin_.reset(new oops::Variables(vv));
  oops::Log::info() << "ObsGnssroBndROPP2DTLAD vars: " << *varin_ << std::endl;
  oops::Log::trace() << "ObsGnssroBndROPP2DTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndROPP2DTLAD::~ObsGnssroBndROPP2DTLAD() {
  ufo_gnssro_bndropp2d_tlad_delete_f90(keyOperGnssroBndROPP2D_);
  oops::Log::trace() << "ObsGnssroBndROPP2DTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2DTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                           const QCFlags_t & qc_flags) {
  ufo_gnssro_bndropp2d_tlad_settraj_f90(keyOperGnssroBndROPP2D_, geovals.toFortran(), obsspace());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2DTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                           const QCFlags_t & qc_flags) const {
  ufo_gnssro_bndropp2d_simobs_tl_f90(keyOperGnssroBndROPP2D_, geovals.toFortran(), obsspace(),
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2DTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                           const QCFlags_t & qc_flags) const {
  ufo_gnssro_bndropp2d_simobs_ad_f90(keyOperGnssroBndROPP2D_, geovals.toFortran(), obsspace(),
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2DTLAD::print(std::ostream & os) const {
  os << "ObsGnssroBndROPP2DTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
