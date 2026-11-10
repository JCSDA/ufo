/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/crtm/ObsAodCRTM.h"

#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/dateFunctions.h"
#include "oops/util/TimeWindow.h"

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

  // get a single central of middle time from observation space
  const util::DateTime midPoint = odb.timeWindow().midpoint();
  std::string year, month, day, hour, minute, second;
  midPoint.toYYYYMMDDhhmmss(year, month, day, hour,  minute,  second);
  // Julian Day Number since noon Universal Time (UT) on January 1, 4713 BCE
  uint64_t midPointJulday = util::datefunctions::dateToJulian(std::stoi(year),
                                                              std::stoi(month),
                                                              std::stoi(day));

  // call Fortran setup routine
  ufo_aodcrtm_setup_f90(keyOperAodCRTM_, parameters_.toConfiguration(),
                        channels_list.size(), channels_list[0], midPointJulday,
                        varin_, odb.comm());
  oops::Log::info() << "ObsAodCRTM variables: " << varin_ << std::endl;
  oops::Log::info() << "ObsAodCRTM channels: " << channels_list << std::endl;
  oops::Log::trace() << "ObsAodCRTM constructor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodCRTM::~ObsAodCRTM() {
  ufo_aodcrtm_delete_f90(keyOperAodCRTM_);
  oops::Log::trace() << "ObsAodCRTM destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCRTM::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  ufo_aodcrtm_simobs_f90(keyOperAodCRTM_, gom.toFortran(), odb_,
                          ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsAodCRTM::simulateObsAD done" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCRTM::print(std::ostream & os) const {
  os << "ObsAodCRTM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
