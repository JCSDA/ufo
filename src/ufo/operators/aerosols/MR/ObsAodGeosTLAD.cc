/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/aerosols/MR/ObsAodGeosTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"

#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAodGeosTLAD> makerAodGeosTL_("AodGeos");
// -----------------------------------------------------------------------------

ObsAodGeosTLAD::ObsAodGeosTLAD(const ioda::ObsSpace & odb,
                               const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOper_(0), varin_()
{
  ufo_aodgeos_tlad_setup_f90(keyOper_, params.toConfiguration(),
                             odb.obsvariables(), varin_);

  oops::Log::trace() << "ObsAodGeosTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodGeosTLAD::~ObsAodGeosTLAD() {
  ufo_aodgeos_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsAodGeosTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodGeosTLAD::setTrajectory(const GeoVaLs & geovals,
                                   ObsDiagnostics &) {
  oops::Log::trace() << "ObsAodGeosTLAD: trajectory entering" << std::endl;
  ufo_aodgeos_tlad_settraj_f90(keyOper_, geovals.toFortran(), obsspace());
  oops::Log::trace() << "ObsAodGeosTLAD: set trajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodGeosTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                   const QCFlags_t & qc_flags) const {
  ufo_aodgeos_simobs_tl_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsAodGeosTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodGeosTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                   const QCFlags_t & qc_flags) const {
  ufo_aodgeos_simobs_ad_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsAodGeosTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodGeosTLAD::print(std::ostream & os) const {
  os << "ObsAodGeosTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
