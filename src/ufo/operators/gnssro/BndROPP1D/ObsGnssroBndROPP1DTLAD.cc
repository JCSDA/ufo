/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/gnssro/BndROPP1D/ObsGnssroBndROPP1DTLAD.h"

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
static LinearObsOperatorMaker<ObsGnssroBndROPP1DTLAD> makerGnssroBndROPP1DTL_("GnssroBndROPP1D");
// -----------------------------------------------------------------------------

ObsGnssroBndROPP1DTLAD::ObsGnssroBndROPP1DTLAD(const ioda::ObsSpace & odb,
                                               const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOperGnssroBndROPP1D_(0), varin_()
{
  ufo_gnssro_bndropp1d_tlad_setup_f90(keyOperGnssroBndROPP1D_,
                                      params.options.value().toConfiguration());
  const std::vector<std::string> vv{"air_temperature", "specific_humidity", "air_pressure"};

  varin_.reset(new oops::Variables(vv));
  oops::Log::info() << "ObsGnssroBndROPP1DTLAD vars: " << *varin_ << std::endl;
  oops::Log::trace() << "ObsGnssroBndROPP1DTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndROPP1DTLAD::~ObsGnssroBndROPP1DTLAD() {
  ufo_gnssro_bndropp1d_tlad_delete_f90(keyOperGnssroBndROPP1D_);
  oops::Log::trace() << "ObsGnssroBndROPP1DTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP1DTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                           const QCFlags_t & qc_flags) {
  ufo_gnssro_bndropp1d_tlad_settraj_f90(keyOperGnssroBndROPP1D_, geovals.toFortran(), obsspace());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP1DTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                           const QCFlags_t & qc_flags) const {
  ufo_gnssro_bndropp1d_simobs_tl_f90(keyOperGnssroBndROPP1D_, geovals.toFortran(), obsspace(),
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP1DTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                           const QCFlags_t & qc_flags) const {
  ufo_gnssro_bndropp1d_simobs_ad_f90(keyOperGnssroBndROPP1D_, geovals.toFortran(), obsspace(),
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP1DTLAD::print(std::ostream & os) const {
  os << "ObsGnssroBndROPP1DTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
