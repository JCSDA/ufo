/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/marinevertinterp/ObsMarineVertInterpTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsMarineVertInterpTLAD> makerMarinevertinterpTL_("MarineVertInterp");
// -----------------------------------------------------------------------------

ObsMarineVertInterpTLAD::ObsMarineVertInterpTLAD(const ioda::ObsSpace & odb,
                                                 const ObsMarineVertInterpTLADParameters & params)
  : LinearObsOperatorBase(odb), keyOper_(0), varin_()
{
  ufo_marinevertinterp_tlad_setup_f90(keyOper_, params.toConfiguration(),
                                      odb.assimvariables(), varin_);
  oops::Log::trace() << "ObsMarineVertInterpTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsMarineVertInterpTLAD::~ObsMarineVertInterpTLAD() {
  ufo_marinevertinterp_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsMarineVertInterpTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterpTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &) {
  ufo_marinevertinterp_tlad_settraj_f90(keyOper_, geovals.toFortran(), obsspace());
  oops::Log::trace() << "ObsMarineVertInterpTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterpTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_marinevertinterp_simobs_tl_f90(keyOper_, geovals.toFortran(), obsspace(),
                                      ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsMarineVertInterpTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterpTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_marinevertinterp_simobs_ad_f90(keyOper_, geovals.toFortran(), obsspace(),
                                      ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsMarineVertInterpTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterpTLAD::print(std::ostream & os) const {
  os << "ObsMarinevertinterpTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
