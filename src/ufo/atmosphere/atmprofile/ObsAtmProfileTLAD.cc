/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/atmosphere/atmprofile/ObsAtmProfileTLAD.h"

#include <ostream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

#include <boost/scoped_ptr.hpp>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAtmProfileTLAD> makerRadiosondeTL_("Radiosonde");
static LinearObsOperatorMaker<ObsAtmProfileTLAD> makerAircraftTL_("Aircraft");
static LinearObsOperatorMaker<ObsAtmProfileTLAD> makerSatwindTL_("Satwind");
// -----------------------------------------------------------------------------

ObsAtmProfileTLAD::ObsAtmProfileTLAD(const ioda::ObsSpace & odb,
                                     const eckit::Configuration & config)
  : keyOperAtmProfile_(0), varin_(), odb_(odb)
{
  int c_name_size = 200;
  char *buffin = new char[c_name_size];
  char *buffout = new char[c_name_size];
  const eckit::Configuration * configc = &config;

  ufo_atmprofile_tlad_setup_f90(keyOperAtmProfile_, &configc, buffin, buffout, c_name_size);

  std::string vstr_in(buffin), vstr_out(buffout);
  std::vector<std::string> vvin;
  std::vector<std::string> vvout;
  boost::split(vvin, vstr_in, boost::is_any_of("\t"));
  boost::split(vvout, vstr_out, boost::is_any_of("\t"));
  varin_.reset(new oops::Variables(vvin));
  varout_.reset(new oops::Variables(vvout));

  oops::Log::trace() << "ObsAtmProfileTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmProfileTLAD::~ObsAtmProfileTLAD() {
  ufo_atmprofile_tlad_delete_f90(keyOperAtmProfile_);
  oops::Log::trace() << "ObsAtmProfileTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmProfileTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_atmprofile_tlad_settraj_f90(keyOperAtmProfile_, geovals.toFortran(), odb_);
}

// -----------------------------------------------------------------------------

void ObsAtmProfileTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                      const ObsBiasIncrement & bias) const {
  ufo_atmprofile_simobs_tl_f90(keyOperAtmProfile_, geovals.toFortran(), odb_,
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAtmProfileTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                      ObsBiasIncrement & bias) const {
  ufo_atmprofile_simobs_ad_f90(keyOperAtmProfile_, geovals.toFortran(), odb_,
                               ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsAtmProfileTLAD::print(std::ostream & os) const {
  os << "ObsAtmProfileTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
