/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/crtm/ObsAodCRTMTLAD.h"

#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/util/dateFunctions.h"
#include "oops/util/Logger.h"
#include "oops/util/TimeWindow.h"
#include "ufo/GeoVaLs.h"
#include "ufo/operators/crtm/ObsAodCRTMTLAD.interface.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAodCRTMTLAD> makerAodTL_("AodCRTM");
// -----------------------------------------------------------------------------

ObsAodCRTMTLAD::ObsAodCRTMTLAD(const ioda::ObsSpace & odb,
                               const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOperAodCRTM_(0), varin_()
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

  ufo_aodcrtm_tlad_setup_f90(keyOperAodCRTM_, params.toConfiguration(),
                             channels_list.size(), channels_list[0], midPointJulday,
                             varin_);
  oops::Log::info() << "ObsAodCRTMTLAD variables: " << varin_ << std::endl;
  oops::Log::info() << "ObsAodCRTMTLAD channels: " << channels_list << std::endl;
  oops::Log::trace() << "ObsAodCRTMTLAD constructor done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodCRTMTLAD::~ObsAodCRTMTLAD() {
  ufo_aodcrtm_tlad_delete_f90(keyOperAodCRTM_);
  oops::Log::trace() << "ObsAodCRTMTLAD destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCRTMTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                   const QCFlags_t & qc_flags) {
  ufo_aodcrtm_tlad_settraj_f90(keyOperAodCRTM_, geovals.toFortran(), obsspace());
  oops::Log::trace() << "ObsAodCRTMTLAD::setTrajectory done" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCRTMTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_aodcrtm_simobs_tl_f90(keyOperAodCRTM_, geovals.toFortran(), obsspace(),
                             ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsAodCRTMTLAD::simulateObsTL done" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCRTMTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_aodcrtm_simobs_ad_f90(keyOperAodCRTM_, geovals.toFortran(), obsspace(),
                             ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsAodCRTMTLAD::simulateObsAD done" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCRTMTLAD::print(std::ostream & os) const {
  os << "ObsAodCRTMTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
