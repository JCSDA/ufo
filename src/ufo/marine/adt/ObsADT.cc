/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/marine/adt/ObsADT.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"


namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsADT> makerADT_("ADT");
// -----------------------------------------------------------------------------

ObsADT::ObsADT(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_(), varout_()
{
  const std::vector<std::string> vvin{"sea_surface_height_above_geoid"};
  varin_.reset(new oops::Variables(vvin));
  const std::vector<std::string> vvout{"obs_absolute_dynamic_topography"};
  varout_.reset(new oops::Variables(vvout));
  const eckit::Configuration * configc = &config;
  ufo_adt_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsADT created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsADT::~ObsADT() {
  ufo_adt_delete_f90(keyOper_);
  oops::Log::trace() << "ObsADT destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADT::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                         const ObsBias & bias) const {
  ufo_adt_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran(),
                     bias.toFortran());
  oops::Log::trace() << "ObsADT: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

Locations * ObsADT::locateObs(const util::DateTime & t1,
                              const util::DateTime & t2) const {
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  int keylocs;
  ufo_adt_locateobs_f90(keyOper_, odb_, &p1, &p2, keylocs);

  return new Locations(keylocs);
}

// -----------------------------------------------------------------------------

void ObsADT::print(std::ostream & os) const {
  os << "ObsADT::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
