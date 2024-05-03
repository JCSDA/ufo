/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
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
#include "ufo/operators/gnssro/BendMetOffice/ObsGnssroBendMetOfficeTLAD.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsGnssroBendMetOfficeTLAD>
    makerGnssroBendMetOfficeTL_("GnssroBendMetOffice");
// -----------------------------------------------------------------------------

ObsGnssroBendMetOfficeTLAD::ObsGnssroBendMetOfficeTLAD(const ioda::ObsSpace & odb,
                                               const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), keyOperGnssroBendMetOffice_(0), varin_()
{
  const std::vector<std::string> vv{"air_pressure_levels", "specific_humidity"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::info() << "ObsGnssroBendMetOfficeTLAD vars: " << *varin_ << std::endl;

  std::set<int> channelset = oops::parseIntSet(parameters.channelList);
  std::vector<int> channels;
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels));

  ufo_gnssro_bendmetoffice_tlad_setup_f90(keyOperGnssroBendMetOffice_,
                                          parameters.vertInterpOPS,
                                          parameters.pseudoLevels,
                                          parameters.minTempGrad,
                                          channels.size(),
                                          channels[0],
                                          parameters.noSuperCheck);

  oops::Log::trace() << "ObsGnssroBendMetOfficeTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBendMetOfficeTLAD::~ObsGnssroBendMetOfficeTLAD() {
  ufo_gnssro_bendmetoffice_tlad_delete_f90(keyOperGnssroBendMetOffice_);
  oops::Log::trace() << "ObsGnssroBendMetOfficeTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBendMetOfficeTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                               const QCFlags_t & qc_flags) {
  ufo_gnssro_bendmetoffice_tlad_settraj_f90(keyOperGnssroBendMetOffice_, geovals.toFortran(),
                                            obsspace());
}

// -----------------------------------------------------------------------------

void ObsGnssroBendMetOfficeTLAD::simulateObsTL(
        const GeoVaLs & geovals, ioda::ObsVector & ovec, const QCFlags_t & qc_flags) const {
  ufo_gnssro_bendmetoffice_simobs_tl_f90(keyOperGnssroBendMetOffice_, geovals.toFortran(),
                                         obsspace(), ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBendMetOfficeTLAD::simulateObsAD(
        GeoVaLs & geovals, const ioda::ObsVector & ovec, const QCFlags_t & qc_flags) const {
  ufo_gnssro_bendmetoffice_simobs_ad_f90(keyOperGnssroBendMetOffice_, geovals.toFortran(),
                                         obsspace(), ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroBendMetOfficeTLAD::print(std::ostream & os) const {
  os << "ObsGnssroBendMetOfficeTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
