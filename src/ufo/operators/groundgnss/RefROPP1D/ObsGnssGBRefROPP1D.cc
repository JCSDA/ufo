/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/groundgnss/RefROPP1D/ObsGnssGBRefROPP1D.h"

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
static ObsOperatorMaker<ObsGnssGBRefROPP1D> makerGnssGBRefROPP1D_("GnssGBRefROPP1D");
// -----------------------------------------------------------------------------

ObsGnssGBRefROPP1D::ObsGnssGBRefROPP1D(const ioda::ObsSpace & odb,
                                       const Parameters_ & params)
  : ObsOperatorBase(odb), keyOperGnssGBRefROPP1D_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vv{"air_temperature", "specific_humidity", "air_pressure",
                                    "geopotential_height", "surface_altitude"};
  varin_.reset(new oops::Variables(vv));

  ufo_gnssgb_refropp1d_setup_f90(keyOperGnssGBRefROPP1D_, params.toConfiguration());
  oops::Log::trace() << "ObsGnssGBRefROPP1D created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssGBRefROPP1D::~ObsGnssGBRefROPP1D() {
  ufo_gnssgb_refropp1d_delete_f90(keyOperGnssGBRefROPP1D_);
  oops::Log::trace() << "ObsGnssGBRefROPP1D destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssGBRefROPP1D::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                     ObsDiagnostics &, const QCFlags_t & qc_flags) const
{
  ufo_gnssgb_refropp1d_simobs_f90(keyOperGnssGBRefROPP1D_, gom.toFortran(), odb_,
                                  ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssGBRefROPP1D::print(std::ostream & os) const {
  os << "ObsGnssGBRefROPP1D::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
