/*
 * (C) Copyright 2024 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/radardopplerwind/ObsRadarDopplerWind.h"

#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/radardopplerwind/ObsRadarDopplerWindParameters.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadarDopplerWind> obsRadarDopplerWindMaker_("RadarDopplerWind");
// -----------------------------------------------------------------------------

ObsRadarDopplerWind::ObsRadarDopplerWind(const ioda::ObsSpace & odb,
                                         const Parameters_ & params)
  : ObsOperatorBase(odb), odb_(odb),
    data_(odb, params)
{
  oops::Log::trace() << "ObsRadarDopplerWind constructed" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadarDopplerWind::~ObsRadarDopplerWind() {
  oops::Log::trace() << "ObsRadarDopplerWind destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWind::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                      ObsDiagnostics &,
                                      const QCFlags_t &) const {
  oops::Log::trace() << "ObsRadarDopplerWind: simulateObs entered" << std::endl;

  data_.cacheVertCoordGeoVaLs(gv);
  data_.fillHofX(odb_, gv, ovec);

  oops::Log::trace() << "ObsRadarDopplerWind: simulateObs finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWind::print(std::ostream & os) const {
  data_.print(os);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
