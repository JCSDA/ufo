/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/geos_aero/ObsGeosAodTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"

#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsGeosAodTLAD> makerGeosAodTL_("GeosAod");
// -----------------------------------------------------------------------------

ObsGeosAodTLAD::ObsGeosAodTLAD(const ioda::ObsSpace & odb,
                               const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_()
{
  ufo_geosaod_tlad_setup_f90(keyOper_, config, odb.obsvariables(), varin_);

  oops::Log::trace() << "ObsGeosAodTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGeosAodTLAD::~ObsGeosAodTLAD() {
  ufo_geosaod_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsGeosAodTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGeosAodTLAD::setTrajectory(const GeoVaLs & geovals,
                                   const ObsBias & bias, ObsDiagnostics &) {
  oops::Log::trace() << "ObsGeosAodTLAD: trajectory entering" << std::endl;
  ufo_geosaod_tlad_settraj_f90(keyOper_, geovals.toFortran(), odb_);
  oops::Log::trace() << "ObsGeosAodTLAD: set trajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGeosAodTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_geosaod_simobs_tl_f90(keyOper_, geovals.toFortran(), odb_,
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsGeosAodTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGeosAodTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_geosaod_simobs_ad_f90(keyOper_, geovals.toFortran(), odb_,
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsGeosAodTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGeosAodTLAD::print(std::ostream & os) const {
  os << "ObsGeosAodTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
