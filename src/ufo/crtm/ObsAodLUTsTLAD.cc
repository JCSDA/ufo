/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/crtm/ObsAodLUTsTLAD.h"

#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAodLUTsTLAD> makerAodTL_("AodLUTs");
// -----------------------------------------------------------------------------

ObsAodLUTsTLAD::ObsAodLUTsTLAD(const ioda::ObsSpace & odb,
                               const eckit::Configuration & config)
  : LinearObsOperatorBase(odb), keyOperAodLUTs_(0), varin_()
{
  // parse channels from the config and create variable names
  const oops::Variables & observed = odb.obsvariables();
  std::vector<int> channels_list = observed.channels();

  ufo_aodluts_tlad_setup_f90(keyOperAodLUTs_, config,
                             channels_list.size(), channels_list[0], varin_);
  oops::Log::info() << "ObsAodLUTsTLAD variables: " << varin_ << std::endl;
  oops::Log::info() << "ObsAodLUTsTLAD channels: " << channels_list << std::endl;
  oops::Log::trace() << "ObsAodLUTsTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodLUTsTLAD::~ObsAodLUTsTLAD() {
  ufo_aodluts_tlad_delete_f90(keyOperAodLUTs_);
  oops::Log::trace() << "ObsAodLUTsTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodLUTsTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                   ObsDiagnostics &) {
  ufo_aodluts_tlad_settraj_f90(keyOperAodLUTs_, geovals.toFortran(), obsspace());
}

// -----------------------------------------------------------------------------

void ObsAodLUTsTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_aodluts_simobs_tl_f90(keyOperAodLUTs_, geovals.toFortran(), obsspace(),
                             ovec.nvars(), ovec.nlocs(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAodLUTsTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_aodluts_simobs_ad_f90(keyOperAodLUTs_, geovals.toFortran(), obsspace(),
                             ovec.nvars(), ovec.nlocs(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAodLUTsTLAD::print(std::ostream & os) const {
  os << "ObsAodLUTsTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
