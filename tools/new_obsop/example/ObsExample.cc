/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/example/ObsExample.h"

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
static ObsOperatorMaker<ObsExample> makerExample_("Example");
// -----------------------------------------------------------------------------

ObsExample::ObsExample(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_(), varout_()
{
  // TODO(anyone): list the variables for GeoVaLs that are needed for the observation
  //       operator below in vv (e.g., vv{"temperature", "humidity"})
  const std::vector<std::string> vvin{""};
  varin_.reset(new oops::Variables(vvin));
  const std::vector<std::string> vvout{""};
  varout_.reset(new oops::Variables(vvout));
  const eckit::Configuration * configc = &config;
  ufo_example_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsExample created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsExample::~ObsExample() {
  ufo_example_delete_f90(keyOper_);
  oops::Log::trace() << "ObsExample destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExample::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                              const ObsBias & bias) const {
  ufo_example_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran(),
                      bias.toFortran());
  oops::Log::trace() << "ObsExample: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

Locations * ObsExample::locateObs(const util::DateTime & t1,
                                   const util::DateTime & t2) const {
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  int keylocs;
  ufo_example_locateobs_f90(keyOper_, odb_, &p1, &p2, keylocs);

  return new Locations(keylocs);
}

// -----------------------------------------------------------------------------

void ObsExample::print(std::ostream & os) const {
  os << "ObsExample::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
