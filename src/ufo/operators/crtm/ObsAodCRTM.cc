/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/crtm/ObsAodCRTM.h"

#include <vector>

#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/crtm/ObsAodCRTM.interface.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAodCRTM> makerAOD_("AodCRTM");

// -----------------------------------------------------------------------------

ObsAodCRTM::ObsAodCRTM(const ioda::ObsSpace & odb,
                       const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOperAodCRTM_(0), odb_(odb), varin_(),
    parameters_(parameters)
{
  // parse channels from the config and create variable names
  const oops::ObsVariables & observed = odb.assimvariables();
  std::vector<int> channels_list = observed.channels();

  // call Fortran setup routine
  ufo_aodcrtm_setup_f90(keyOperAodCRTM_, parameters_.toConfiguration(),
                        channels_list.size(), channels_list[0], varin_);
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
                             ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  ufo_aodcrtm_simobs_f90(keyOperAodCRTM_, gom.toFortran(), odb_,
                          ovec.nvars(), ovec.nlocs(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAodCRTM::print(std::ostream & os) const {
  os << "ObsAodCRTM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
