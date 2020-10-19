/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/crtm/ObsAodLUTs.h"

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
static ObsOperatorMaker<ObsAodLUTs> makerAOD_("AodLUTs");

// -----------------------------------------------------------------------------

ObsAodLUTs::ObsAodLUTs(const ioda::ObsSpace & odb,
                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperAodLUTs_(0), odb_(odb), varin_()
{
  // parse channels from the config and create variable names
  const oops::Variables & observed = odb.obsvariables();
  std::vector<int> channels_list = observed.channels();

  // call Fortran setup routine
  ufo_aodluts_setup_f90(keyOperAodLUTs_, config, channels_list.size(), channels_list[0], varin_);
  oops::Log::info() << "ObsAodLUTs variables: " << varin_ << std::endl;
  oops::Log::info() << "ObsAodLUTs channels: " << channels_list << std::endl;
  oops::Log::trace() << "ObsAodLUTs created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodLUTs::~ObsAodLUTs() {
  ufo_aodluts_delete_f90(keyOperAodLUTs_);
  oops::Log::trace() << "ObsAodLUTs destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodLUTs::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             ObsDiagnostics &) const {
  ufo_aodluts_simobs_f90(keyOperAodLUTs_, gom.toFortran(), odb_,
                          ovec.nvars(), ovec.nlocs(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAodLUTs::print(std::ostream & os) const {
  os << "ObsAodLUTs::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
