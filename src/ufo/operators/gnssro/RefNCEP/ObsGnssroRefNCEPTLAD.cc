/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/gnssro/RefNCEP/ObsGnssroRefNCEPTLAD.h"

#include <ostream>
#include <string>
#include <vector>


#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// ----------------------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsGnssroRefNCEPTLAD> makerGnssroRefNCEPTL_("GnssroRefNCEP");
// ----------------------------------------------------------------------------------------

ObsGnssroRefNCEPTLAD::ObsGnssroRefNCEPTLAD(const ioda::ObsSpace & odb,
                                           const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOperGnssroRefNCEP_(0), varin_()
{
  ufo_gnssro_refncep_tlad_setup_f90(keyOperGnssroRefNCEP_,
                                    params.options.value().toConfiguration());
  const std::vector<std::string> vv{"air_temperature", "specific_humidity", "air_pressure"};

  varin_.reset(new oops::Variables(vv));
  oops::Log::info() << "ObsGnssroRefNCEPTLAD vars: " << *varin_ << std::endl;
  oops::Log::trace() << "ObsGnssroRefNCEPTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroRefNCEPTLAD::~ObsGnssroRefNCEPTLAD() {
  ufo_gnssro_refncep_tlad_delete_f90(keyOperGnssroRefNCEP_);
  oops::Log::trace() << "ObsGnssroRefNCEPTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroRefNCEPTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                         const QCFlags_t & qc_flags) {
  ufo_gnssro_refncep_tlad_settraj_f90(keyOperGnssroRefNCEP_, geovals.toFortran(), obsspace());
}

// -----------------------------------------------------------------------------

void ObsGnssroRefNCEPTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                         const QCFlags_t & qc_flags) const {
  ufo_gnssro_refncep_simobs_tl_f90(keyOperGnssroRefNCEP_, geovals.toFortran(), obsspace(),
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroRefNCEPTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                         const QCFlags_t & qc_flags) const {
  ufo_gnssro_refncep_simobs_ad_f90(keyOperGnssroRefNCEP_, geovals.toFortran(), obsspace(),
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroRefNCEPTLAD::print(std::ostream & os) const {
  os << "ObsGnssroRefNCEPTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
