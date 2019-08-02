/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/crtm/ObsAodCRTMTLAD.h"

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
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAodCRTMTLAD> makerAodTL_("Aod");
// -----------------------------------------------------------------------------

ObsAodCRTMTLAD::ObsAodCRTMTLAD(const ioda::ObsSpace & odb,
                               const eckit::Configuration & config)
  : keyOperAodCRTM_(0), varin_(), odb_(odb), channels_(odb.obsvariables().channels())
{
  const std::vector<std::string> vv{
    "sulf", "bc1", "bc2", "oc1", "oc2", "dust1", "dust2", "dust3", "dust4", "dust5",
    "seas1", "seas2", "seas3", "seas4"};
  varin_.reset(new oops::Variables(vv));

  const eckit::LocalConfiguration obsOpts(config, "ObsOptions");
  const eckit::Configuration * configOpts = &obsOpts;
  const eckit::Configuration * configOper = &config;
  ufo_aodcrtm_tlad_setup_f90(keyOperAodCRTM_, &configOpts, &configOper);
  oops::Log::trace() << "ObsAodCRTMTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodCRTMTLAD::~ObsAodCRTMTLAD() {
  ufo_aodcrtm_tlad_delete_f90(keyOperAodCRTM_);
  oops::Log::trace() << "ObsAodCRTMTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCRTMTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_aodcrtm_tlad_settraj_f90(keyOperAodCRTM_, geovals.toFortran(), odb_,
                                channels_.size(), channels_[0]);
}

// -----------------------------------------------------------------------------

void ObsAodCRTMTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                    const ObsBiasIncrement & bias) const {
  ufo_aodcrtm_simobs_tl_f90(keyOperAodCRTM_, geovals.toFortran(), odb_,
                             ovec.size(), ovec.toFortran(),
                             channels_.size(), channels_[0]);
}

// -----------------------------------------------------------------------------

void ObsAodCRTMTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                    ObsBiasIncrement & bias) const {
  ufo_aodcrtm_simobs_ad_f90(keyOperAodCRTM_, geovals.toFortran(), odb_,
                             ovec.size(), ovec.toFortran(),
                             channels_.size(), channels_[0]);
}

// -----------------------------------------------------------------------------

void ObsAodCRTMTLAD::print(std::ostream & os) const {
  os << "ObsAodCRTMTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
