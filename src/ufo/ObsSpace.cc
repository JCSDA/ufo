/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ObsSpace.h"

#include <map>
#include <string>

#include "eckit/config/Configuration.h"

#include "util/abor1_cpp.h"
#include "util/Logger.h"

//#include "ObsVector.h"

namespace ufo {
// -----------------------------------------------------------------------------

ObsSpace::ObsSpace(const eckit::Configuration & config,
                   const util::DateTime & bgn, const util::DateTime & end)
  : oops::ObsSpaceBase(config, bgn, end), winbgn_(bgn), winend_(end)
{
  oops::Log::trace() << "ufo::ObsSpace config  = " << config << std::endl;

  const eckit::Configuration * configc = &config;
  ufo_obsdb_setup_f90(keyOspace_, &configc);

  obsname_ = config.getString("ObsType");
  oops::Log::trace() << "ufo::ObsSpace contructed name = " << obsname_ << std::endl;
}

// -----------------------------------------------------------------------------

//void ObsSpace::printJo(const ObsVector &, const ObsVector &) {
//  oops::Log::info() << "ObsSpace::printJo not implemented" << std::endl;
//}

// -----------------------------------------------------------------------------

ObsSpace::~ObsSpace() {
}

// -----------------------------------------------------------------------------

void ObsSpace::print(std::ostream & os) const {
  os << "ObsSpace::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
