/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/crtm/ObsRadianceCRTM.h"

#include <vector>

#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/crtm/ObsRadianceCRTM.interface.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadianceCRTM> makerCRTM_("CRTM");

// -----------------------------------------------------------------------------

ObsRadianceCRTM::ObsRadianceCRTM(const ioda::ObsSpace & odb,
                                 const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOperRadianceCRTM_(0),
    odb_(odb), varin_(), parameters_(parameters)
{
  // parse channels from the config and create variable names
  const oops::Variables & observed = odb.assimvariables();
  std::vector<int> channels_list = observed.channels();

  // call Fortran setup routine
  ufo_radiancecrtm_setup_f90(keyOperRadianceCRTM_, parameters_.toConfiguration(),
                             channels_list.size(), channels_list[0], varin_, odb.comm());

  oops::Log::info() << "ObsRadianceCRTM channels: " << channels_list << std::endl;
  oops::Log::trace() << "ObsRadianceCRTM created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceCRTM::~ObsRadianceCRTM() {
  ufo_radiancecrtm_delete_f90(keyOperRadianceCRTM_);
  oops::Log::trace() << "ObsRadianceCRTM destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTM::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                  ObsDiagnostics & dvec) const {
  ufo_radiancecrtm_simobs_f90(keyOperRadianceCRTM_, gom.toFortran(), odb_,
                          ovec.nvars(), ovec.nlocs(), ovec.toFortran(),
                          dvec.toFortran());
  oops::Log::trace() << "ObsRadianceCRTM simulateObs done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTM::print(std::ostream & os) const {
  os << "ObsRadianceCRTM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
