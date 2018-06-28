/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

//TODO: through the file replace ObsExampleTLAD with <Your_Obs_Operator_Name>TLAD
#include "ObsExampleTLAD.h"

#include <ostream>
#include <string>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsExampleTLAD::ObsExampleTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOper_(0), varin_(), odb_(odb)
{
  // TODO: list the variables for GeoVaLs that are needed for the observation 
  //       operator below in vv (e.g., vv{"temperature", "humidity"})
  const std::vector<std::string> vv{""};
  varin_.reset(new oops::Variables(vv));
  const eckit::Configuration * configc = &config;
  // TODO: replace ufo_example_tlad_setup_f90 with the call to your Fortran routine
  //       to setup tl/ad obs operator (defined in ObsExampleTLAD.interface.F90)
  ufo_example_tlad_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsExampleTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsExampleTLAD::~ObsExampleTLAD() {
  // TODO: replace ufo_example_tlad_delete_f90 with the call to your Fortran routine
  //       to destruct tl/ad observation operator (defined in ObsExampleTLAD.interface.F90)
  ufo_example_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsExampleTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExampleTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  // TODO: replace ufo_example_tlad_settraj_f90 with the call to your Fortran routine
  //       to set trajectory for tl/ad (defined in ObsExampleTLAD.interface.F90)
  ufo_example_tlad_settraj_f90(keyOper_, geovals.toFortran(), odb_.toFortran());
  oops::Log::trace() << "ObsExampleTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExampleTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                             const ObsBiasIncrement & bias) const {
  // TODO: replace ufo_example_tlad_eqv_tl_f90 with the call to your Fortran routine
  //       to apply tl observation operator (defined in ObsExampleTLAD.interface.F90)
  ufo_example_tlad_eqv_tl_f90(keyOper_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
  oops::Log::trace() << "ObsExampleTLAD: tangent linear observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExampleTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                             ObsBiasIncrement & bias) const {
  // TODO: replace ufo_example_tlad_eqv_ad_f90 with the call to your Fortran routine
  //       to apply ad observation operator (defined in ObsExampleTLAD.interface.F90)
  ufo_example_tlad_eqv_ad_f90(keyOper_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
  oops::Log::trace() << "ObsExampleTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExampleTLAD::print(std::ostream & os) const {
  os << "ObsExampleTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSEXAMPLETLAD_H_
