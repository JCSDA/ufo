/*
 * (C) British Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/scatwind/NeutralMetOffice/ObsScatwindNeutralMetOfficeTLAD.h"

#include <ostream>
#include <vector>

#include "oops/util/Logger.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsScatwindNeutralMetOfficeTLAD>
  makerScatwindNeutralMetOfficeTL_("ScatwindNeutralMetOffice");
// -----------------------------------------------------------------------------

ObsScatwindNeutralMetOfficeTLAD::ObsScatwindNeutralMetOfficeTLAD(const ioda::ObsSpace & odb,
                                                                 const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), keyOperScatwindNeutralMetOffice_(0),
    varin_(), parameters_(parameters)
{
  // parse channels from the config and create variable names
  const std::vector<int> channels_list = odb.assimvariables().channels();

  ufo_scatwind_neutralmetoffice_tlad_setup_f90(keyOperScatwindNeutralMetOffice_,
                                               parameters_.surfaceTypeCheck,
                                               parameters_.surfaceTypeSea,
                                               odb.assimvariables(),
                                               varin_,
                                               channels_list.size(),
                                               channels_list[0],
                                               parameters_.toConfiguration());

  oops::Log::trace() << "ObsScatwindNeutralMetOfficeTLAD created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsScatwindNeutralMetOfficeTLAD::~ObsScatwindNeutralMetOfficeTLAD() {
  ufo_scatwind_neutralmetoffice_tlad_delete_f90(keyOperScatwindNeutralMetOffice_);
  oops::Log::trace() << "ObsScatwindNeutralMetOfficeTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsScatwindNeutralMetOfficeTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                                    const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsScatwindNeutralMetOfficeTLAD::setTrajectory entering" << std::endl;

  ufo_scatwind_neutralmetoffice_tlad_settraj_f90(keyOperScatwindNeutralMetOffice_,
                                                 geovals.toFortran(), obsspace());

  oops::Log::trace() << "ObsScatwindNeutralMetOfficeTLAD::setTrajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsScatwindNeutralMetOfficeTLAD::simulateObsTL(const GeoVaLs & geovals,
                                                    ioda::ObsVector & ovec,
                                                    const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsScatwindNeutralMetOfficeTLAD::simulateObsTL entering" << std::endl;

  ufo_scatwind_neutralmetoffice_simobs_tl_f90(keyOperScatwindNeutralMetOffice_,
                                              geovals.toFortran(), obsspace(),
                                              ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsScatwindNeutralMetOfficeTLAD::simulateObsTL exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsScatwindNeutralMetOfficeTLAD::simulateObsAD(GeoVaLs & geovals,
                                                    const ioda::ObsVector & ovec,
                                                    const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsScatwindNeutralMetOfficeTLAD::simulateObsAD entering" << std::endl;

  ufo_scatwind_neutralmetoffice_simobs_ad_f90(keyOperScatwindNeutralMetOffice_,
                                              geovals.toFortran(), obsspace(),
                                              ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsScatwindNeutralMetOfficeTLAD::simulateObsAD exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsScatwindNeutralMetOfficeTLAD::print(std::ostream & os) const {
  os << "ObsScatwindNeutralMetOfficeTLAD::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
