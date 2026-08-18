/*
 * (C) Copyright 2025-2026 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/crtm/ObsExtCoeffProfCRTM.h"

#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/dateFunctions.h"
#include "oops/util/TimeWindow.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/crtm/ObsExtCoeffProfCRTM.interface.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsExtCoeffProfCRTM>
       makerExtCoeffProfCRTM_("ExtinctionCoefficientProfileCRTM");

// -----------------------------------------------------------------------------

ObsExtCoeffProfCRTM::ObsExtCoeffProfCRTM(const ioda::ObsSpace & odb,
                       const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOperExtCoeffProfCRTM_(0), odb_(odb), varin_(),
    parameters_(parameters)
{
  // parse channels from the config and create variable names
  const oops::ObsVariables & observed = odb.assimvariables();
  std::vector<int> channels_list = observed.channels();

  // get a single central of middle time from observation space
  const util::DateTime midPoint = odb.timeWindow().midpoint();
  std::string year, month, day, hour, minute, second;
  midPoint.toYYYYMMDDhhmmss(year, month, day, hour,  minute,  second);
  // Julian Day Number since noon Universal Time (UT) on January 1, 4713 BCE
  uint64_t midPointJulday = util::datefunctions::dateToJulian(std::stoi(year),
                                                              std::stoi(month),
                                                              std::stoi(day));

  // call Fortran setup routine
  ufo_extcoeffprofcrtm_setup_f90(keyOperExtCoeffProfCRTM_, parameters_.toConfiguration(),
                        channels_list.size(), channels_list[0], midPointJulday,
                        varin_, odb.comm());
  oops::Log::info() << "ObsExtCoeffProfCRTM channels (levels): " << channels_list << std::endl;
  oops::Log::trace() << "ObsExtCoeffProfCRTM constructor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsExtCoeffProfCRTM::~ObsExtCoeffProfCRTM() {
  ufo_extcoeffprofcrtm_delete_f90(keyOperExtCoeffProfCRTM_);
  oops::Log::trace() << "ObsExtCoeffProfCRTM destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExtCoeffProfCRTM::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  int nlevs = parameters_.nProfileLevels.value();
  ufo_extcoeffprofcrtm_simobs_f90(keyOperExtCoeffProfCRTM_, gom.toFortran(), odb_,
                          ovec.nvars(), ovec.nlocs(), nlevs, ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsExtCoeffProfCRTM::print(std::ostream & os) const {
  os << "ObsExtCoeffProfCRTM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
