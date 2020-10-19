/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "tools/new_obsop/example/ObsExampleTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsExampleTLAD> makerExampleTL_("Example");
// -----------------------------------------------------------------------------

ObsExampleTLAD::ObsExampleTLAD(const ioda::ObsSpace & odb,
                               const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_()
{
  ufo_example_tlad_setup_f90(keyOper_, config, odb.obsvariables(), varin_);
  oops::Log::trace() << "ObsExampleTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsExampleTLAD::~ObsExampleTLAD() {
  ufo_example_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsExampleTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExampleTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                   ObsDiagnostics & ydiags) {
  ufo_example_tlad_settraj_f90(keyOper_, geovals.toFortran(), odb_, ydiags.toFortran());
  oops::Log::trace() << "ObsExampleTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExampleTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_example_simobs_tl_f90(keyOper_, geovals.toFortran(), odb_,
                            ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsExampleTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExampleTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_example_simobs_ad_f90(keyOper_, geovals.toFortran(), odb_,
                            ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsExampleTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExampleTLAD::print(std::ostream & os) const {
  os << "ObsExampleTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
