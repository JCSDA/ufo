/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/groundgnss/ZenithTotalDelayROPP/ObsGroundgnssROPP.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGroundgnssROPP> makerGroundgnssROPP_("GroundgnssROPP");
// -----------------------------------------------------------------------------

ObsGroundgnssROPP::ObsGroundgnssROPP(const ioda::ObsSpace & odb,
                                     const Parameters_ & params)
  : ObsOperatorBase(odb), keyOperGroundgnssROPP_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vv{"air_temperature", "specific_humidity", "air_pressure",
                                    "geopotential_height", "surface_altitude"};
  varin_.reset(new oops::Variables(vv));

  ufo_groundgnss_ropp_setup_f90(keyOperGroundgnssROPP_, params.toConfiguration());
  oops::Log::trace() << "ObsGroundgnssROPP created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGroundgnssROPP::~ObsGroundgnssROPP() {
  ufo_groundgnss_ropp_delete_f90(keyOperGroundgnssROPP_);
  oops::Log::trace() << "ObsGroundgnssROPP destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGroundgnssROPP::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                     ObsDiagnostics & d,
                                     const QCFlags_t & qc_flags) const {
  ufo_groundgnss_ropp_simobs_f90(keyOperGroundgnssROPP_, gom.toFortran(), odb_,
                                  ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGroundgnssROPP::print(std::ostream & os) const {
  os << "ObsGroundgnssROPP::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
