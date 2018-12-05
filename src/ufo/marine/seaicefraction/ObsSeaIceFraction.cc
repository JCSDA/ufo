/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/marine/seaicefraction/ObsSeaIceFraction.h"

#include <ostream>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsSeaIceFraction> makerSeaIceFraction_("SeaIceFraction");
// -----------------------------------------------------------------------------

ObsSeaIceFraction::ObsSeaIceFraction(const ioda::ObsSpace & odb,
                                     const eckit::Configuration & config)
  : keyOperSeaIceFraction_(0), odb_(odb), varin_(), varout_()
{
  const eckit::Configuration * configc = &config;
  ufo_seaicefraction_setup_f90(keyOperSeaIceFraction_, &configc);

  const std::vector<std::string> vv{"ice_concentration"};
  varin_.reset(new oops::Variables(vv));

  const std::vector<std::string> vout{"zz"};
  varout_.reset(new oops::Variables(vout));

  oops::Log::trace() << "ObsSeaIceFraction created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceFraction::~ObsSeaIceFraction() {
  ufo_seaicefraction_delete_f90(keyOperSeaIceFraction_);
  oops::Log::trace() << "ObsSeaIceFraction destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFraction::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             const ObsBias & bias) const {
  ufo_seaicefraction_simobs_f90(keyOperSeaIceFraction_, gom.toFortran(), odb_,
                         ovec.size(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------

Locations * ObsSeaIceFraction::locateObs(const util::DateTime & t1,
                                   const util::DateTime & t2) const {
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  int keylocs;
  ufo_seaicefraction_locateobs_f90(keyOper_, odb_, &p1, &p2, keylocs);

  return new Locations(keylocs);
}

// -----------------------------------------------------------------------------

void ObsSeaIceFraction::print(std::ostream & os) const {
  os << "ObsSeaIceFraction::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
