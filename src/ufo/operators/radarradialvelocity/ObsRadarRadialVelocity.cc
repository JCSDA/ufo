/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/radarradialvelocity/ObsRadarRadialVelocity.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadarRadialVelocity> makerRadarRadialVelocity_("RadarRadialVelocity");
// -----------------------------------------------------------------------------

ObsRadarRadialVelocity::ObsRadarRadialVelocity(const ioda::ObsSpace & odb,
                       const ObsRadarRadialVelocityParameters & params)
  : ObsOperatorBase(odb), keyOper_(0), odb_(odb), varin_()
{
  ufo_radarradialvelocity_setup_f90(keyOper_, params.toConfiguration(),
                                    odb.assimvariables(), varin_);

  oops::Log::trace() << "ObsRadarRadialVelocity created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadarRadialVelocity::~ObsRadarRadialVelocity() {
  ufo_radarradialvelocity_delete_f90(keyOper_);
  oops::Log::trace() << "ObsRadarRadialVelocity destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarRadialVelocity::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                         ObsDiagnostics &, const QCFlags_t & qc_flags) const {
  ufo_radarradialvelocity_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsRadarRadialVelocity: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarRadialVelocity::print(std::ostream & os) const {
  os << "ObsRadarRadialVelocity::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
