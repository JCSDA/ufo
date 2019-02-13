/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/atmosphere/crtm/ObsRadianceCRTMTLAD.h"

#include <ostream>
#include <set>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/utils/IntSetParser.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsRadianceCRTMTLAD> makerCRTMTL_("CRTM");
// -----------------------------------------------------------------------------

ObsRadianceCRTMTLAD::ObsRadianceCRTMTLAD(const ioda::ObsSpace & odb,
                                         const eckit::Configuration & config)
  : keyOperRadianceCRTM_(0), varin_(), odb_(odb)
{
  const std::vector<std::string> vv{"air_temperature"};
  varin_.reset(new oops::Variables(vv));

  // parse channels from the config and create variable names
  std::string chlist = config.getString("channels");
  std::set<int> channels = parseIntSet(chlist);
  channels_.reserve(channels.size());
  for (const int jj : channels) {
    channels_.push_back(jj);
  }

  const eckit::LocalConfiguration obsOptions(config, "ObsOptions");
  const eckit::Configuration * configc = &obsOptions;
  ufo_radiancecrtm_tlad_setup_f90(keyOperRadianceCRTM_, &configc);
  oops::Log::trace() << "ObsRadianceCRTMTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceCRTMTLAD::~ObsRadianceCRTMTLAD() {
  ufo_radiancecrtm_tlad_delete_f90(keyOperRadianceCRTM_);
  oops::Log::trace() << "ObsRadianceCRTMTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_radiancecrtm_tlad_settraj_f90(keyOperRadianceCRTM_, geovals.toFortran(), odb_,
                                channels_.size(), channels_[0]);
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                    const ObsBiasIncrement & bias) const {
  ufo_radiancecrtm_simobs_tl_f90(keyOperRadianceCRTM_, geovals.toFortran(), odb_,
                             ovec.size(), ovec.toFortran(),
                             channels_.size(), channels_[0]);
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                    ObsBiasIncrement & bias) const {
  ufo_radiancecrtm_simobs_ad_f90(keyOperRadianceCRTM_, geovals.toFortran(), odb_,
                             ovec.size(), ovec.toFortran(),
                             channels_.size(), channels_[0]);
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::print(std::ostream & os) const {
  os << "ObsRadianceCRTMTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
