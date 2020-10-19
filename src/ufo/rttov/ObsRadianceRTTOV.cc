/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/rttov/ObsRadianceRTTOV.h"

#include <algorithm>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadianceRTTOV> makerRTTOV_("RTTOV");

// -----------------------------------------------------------------------------

ObsRadianceRTTOV::ObsRadianceRTTOV(const ioda::ObsSpace & odb,
                                 const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperRadianceRTTOV_(0),
    odb_(odb), varin_()
{
  // parse channels from the config and create variable names
  const oops::Variables & observed = odb.obsvariables();
  std::vector<int> channels_list = observed.channels();

  // call Fortran setup routine
  ufo_radiancerttov_setup_f90(keyOperRadianceRTTOV_, config,
                             channels_list.size(), channels_list[0], varin_);

  oops::Log::info() << "ObsRadianceRTTOV channels: " << channels_list << std::endl;
  oops::Log::trace() << "ObsRadianceRTTOV created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceRTTOV::~ObsRadianceRTTOV() {
  ufo_radiancerttov_delete_f90(keyOperRadianceRTTOV_);
  oops::Log::trace() << "ObsRadianceRTTOV destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOV::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                  ObsDiagnostics & dvec) const {
  ufo_radiancerttov_simobs_f90(keyOperRadianceRTTOV_, gom.toFortran(), odb_,
                          ovec.nvars(), ovec.nlocs(), ovec.toFortran(),
                          dvec.toFortran());
  oops::Log::trace() << "ObsRadianceRTTOV simulateObs done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOV::print(std::ostream & os) const {
  os << "ObsRadianceRTTOV::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
