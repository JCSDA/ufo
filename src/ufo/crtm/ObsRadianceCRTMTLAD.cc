/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/crtm/ObsRadianceCRTMTLAD.h"

#include <algorithm>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

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
static LinearObsOperatorMaker<ObsRadianceCRTMTLAD> makerCRTMTL_("CRTM");
// -----------------------------------------------------------------------------

ObsRadianceCRTMTLAD::ObsRadianceCRTMTLAD(const ioda::ObsSpace & odb,
                                         const eckit::Configuration & config)
  : keyOperRadianceCRTM_(0), varin_(), odb_(odb)
{
  int c_name_size = 2000;
  char *buffin = new char[c_name_size];

  // parse channels from the config and create variable names
  const oops::Variables & observed = odb.obsvariables();
  std::vector<int> channels_list = observed.channels();

  // call Fortran setup routine
  const eckit::Configuration * configc = &config;
  ufo_radiancecrtm_tlad_setup_f90(keyOperRadianceCRTM_, &configc, channels_list.size(),
                                  channels_list[0], buffin, c_name_size);

  std::string vstr_in(buffin);
  std::vector<std::string> vvin;
  boost::split(vvin, vstr_in, boost::is_any_of("\t"));
  varin_.reset(new oops::Variables(vvin));

  oops::Log::trace() << "ObsRadianceCRTMTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceCRTMTLAD::~ObsRadianceCRTMTLAD() {
  ufo_radiancecrtm_tlad_delete_f90(keyOperRadianceCRTM_);
  oops::Log::trace() << "ObsRadianceCRTMTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_radiancecrtm_tlad_settraj_f90(keyOperRadianceCRTM_, geovals.toFortran(), odb_);
  oops::Log::trace() << "ObsRadianceCRTMTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                    const ObsBiasIncrement & bias) const {
  ufo_radiancecrtm_simobs_tl_f90(keyOperRadianceCRTM_, geovals.toFortran(), odb_,
                             ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsRadianceCRTMTLAD::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                    ObsBiasIncrement & bias) const {
  ufo_radiancecrtm_simobs_ad_f90(keyOperRadianceCRTM_, geovals.toFortran(), odb_,
                             ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsRadianceCRTMTLAD::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::print(std::ostream & os) const {
  os << "ObsRadianceCRTMTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
