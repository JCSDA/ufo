/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/crtm/ObsAodCRTM.h"

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
static ObsOperatorMaker<ObsAodCRTM> makerAOD_("AodCRTM");

// -----------------------------------------------------------------------------

ObsAodCRTM::ObsAodCRTM(const ioda::ObsSpace & odb,
                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperAodCRTM_(0), odb_(odb), varin_()
{
  // parse channels from the config and create variable names
  const oops::Variables & observed = odb.obsvariables();
  std::vector<int> channels_list = observed.channels();

  // call Fortran setup routine
  ufo_aodcrtm_setup_f90(keyOperAodCRTM_, config, channels_list.size(), channels_list[0], varin_);
  oops::Log::info() << "ObsAodCRTM variables: " << varin_ << std::endl;
  oops::Log::info() << "ObsAodCRTM channels: " << channels_list << std::endl;
  oops::Log::trace() << "ObsAodCRTM created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodCRTM::~ObsAodCRTM() {
  ufo_aodcrtm_delete_f90(keyOperAodCRTM_);
  oops::Log::trace() << "ObsAodCRTM destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCRTM::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             ObsDiagnostics &) const {
  ufo_aodcrtm_simobs_f90(keyOperAodCRTM_, gom.toFortran(), odb_,
                          ovec.nvars(), ovec.nlocs(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAodCRTM::print(std::ostream & os) const {
  os << "ObsAodCRTM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
