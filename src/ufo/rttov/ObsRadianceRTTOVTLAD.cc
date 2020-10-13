/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/rttov/ObsRadianceRTTOVTLAD.h"

#include <algorithm>
#include <ostream>
#include <set>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsRadianceRTTOVTLAD> makerRTTOVTL_("RTTOV");
// -----------------------------------------------------------------------------

ObsRadianceRTTOVTLAD::ObsRadianceRTTOVTLAD(const ioda::ObsSpace & odb,
                                         const eckit::Configuration & config)
  : keyOperRadianceRTTOV_(0), odb_(odb), varin_()
{
  // parse channels from the config and create variable names
  const oops::Variables & observed = odb.obsvariables();
  std::vector<int> channels_list = observed.channels();

  // call Fortran setup routine
  ufo_radiancerttov_tlad_setup_f90(keyOperRadianceRTTOV_, config,
                                  channels_list.size(), channels_list[0], varin_);

  oops::Log::trace() << "ObsRadianceRTTOVTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceRTTOVTLAD::~ObsRadianceRTTOVTLAD() {
  ufo_radiancerttov_tlad_delete_f90(keyOperRadianceRTTOV_);
  oops::Log::trace() << "ObsRadianceRTTOVTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                        ObsDiagnostics & ydiags) {
  ufo_radiancerttov_tlad_settraj_f90(keyOperRadianceRTTOV_, geovals.toFortran(), odb_,
                                    ydiags.toFortran());
  oops::Log::trace() << "ObsRadianceRTTOVTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_radiancerttov_simobs_tl_f90(keyOperRadianceRTTOV_, geovals.toFortran(), odb_,
                             ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsRadianceRTTOVTLAD::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_radiancerttov_simobs_ad_f90(keyOperRadianceRTTOV_, geovals.toFortran(), odb_,
                             ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsRadianceRTTOVTLAD::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVTLAD::print(std::ostream & os) const {
  os << "ObsRadianceRTTOVTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
