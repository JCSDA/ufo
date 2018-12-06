/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/atmosphere/atmprofile/ObsAtmProfile.h"

#include <ostream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAtmProfile> makerRadiosonde_("Radiosonde");
static ObsOperatorMaker<ObsAtmProfile> makerAircraft_("Aircraft");
static ObsOperatorMaker<ObsAtmProfile> makerSatwnd_("Satwind");
// -----------------------------------------------------------------------------

ObsAtmProfile::ObsAtmProfile(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperAtmProfile_(0), odb_(odb), varin_(), varout_()
{
  int c_name_size = 200;
  char *buffin = new char[c_name_size];
  char *buffout = new char[c_name_size];
  const eckit::Configuration * configc = &config;

  ufo_atmprofile_setup_f90(keyOperAtmProfile_, &configc, buffin, buffout, c_name_size);

  std::string vstr_in(buffin), vstr_out(buffout);
  std::vector<std::string> vvin;
  std::vector<std::string> vvout;
  boost::split(vvin, vstr_in, boost::is_any_of("\t"));
  boost::split(vvout, vstr_out, boost::is_any_of("\t"));
  varin_.reset(new oops::Variables(vvin));
  varout_.reset(new oops::Variables(vvout));

  oops::Log::trace() << "ObsAtmProfile created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAtmProfile::~ObsAtmProfile() {
  ufo_atmprofile_delete_f90(keyOperAtmProfile_);
  oops::Log::trace() << "ObsAtmProfile destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAtmProfile::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                const ObsBias & bias) const {
  ufo_atmprofile_simobs_f90(keyOperAtmProfile_, gom.toFortran(), odb_,
                            ovec.size(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

Locations * ObsAtmProfile::locateObs(const util::DateTime & t1,
                                     const util::DateTime & t2) const {
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  int keylocs;
  ufo_atmprofile_locateobs_f90(keyOperAtmProfile_, odb_, &p1, &p2, keylocs);

  return new Locations(keylocs);
}

// -----------------------------------------------------------------------------

void ObsAtmProfile::print(std::ostream & os) const {
  os << "ObsAtmProfile::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
