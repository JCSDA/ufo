/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/radarreflectivity/directZDA/ObsDirectZDATLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsDirectZDATLAD> makerDirectZDATL_("DirectZDA");
// -----------------------------------------------------------------------------

ObsDirectZDATLAD::ObsDirectZDATLAD(const ioda::ObsSpace & odb,
                                   const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOperDirectZDA_(0), varin_()
{
  ufo_directZDA_tlad_setup_f90(keyOperDirectZDA_, params.toConfiguration(),
                                                  odb.obsvariables(), varin_);
  oops::Log::trace() << "ObsDirectZDATLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsDirectZDATLAD::~ObsDirectZDATLAD() {
  ufo_directZDA_tlad_delete_f90(keyOperDirectZDA_);
  oops::Log::trace() << "ObsDirectZDATLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsDirectZDATLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                     const QCFlags_t & qc_flags_t) {
  oops::Log::trace() << "ObsDirectZDATLAD::setTrajectory entering" << std::endl;

  ufo_directZDA_tlad_settraj_f90(keyOperDirectZDA_, geovals.toFortran(), obsspace() );

  oops::Log::trace() << "ObsDirectZDATLAD::setTrajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsDirectZDATLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                     const QCFlags_t & qc_flags_t) const {
  ufo_directZDA_simobs_tl_f90(keyOperDirectZDA_, geovals.toFortran(), obsspace() ,
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsDirectZDATLAD: simulateObsTL exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsDirectZDATLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                     const QCFlags_t & qc_flags) const {
  ufo_directZDA_simobs_ad_f90(keyOperDirectZDA_, geovals.toFortran(), obsspace() ,
                              ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsDirectZDATLAD: simulateObsAD exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsDirectZDATLAD::print(std::ostream & os) const {
  os << "ObsDirectZDATLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
