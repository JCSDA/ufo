/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/AnalyticInit.h"

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// \brief Constructor for tests
AnalyticInit::AnalyticInit(const eckit::Configuration & config) : config_(config)
{ }
// -----------------------------------------------------------------------------
/*! \brief Analytic initialization for GeoVaLs
 *
 * \details This ufo::AnalyticInit constructor was introduced in May, 2018 for use with
 * the interpolation test.   If "analytic_init" is not specified in the
 * configuration then this does nothing.  If "analytic_init" **is** specified, then
 * the values are replaced by values computed directly from one of several idealized
 * analytic states.
 *
 * \date May, 2018: Created (M. Miesch, JCSDA)
 * \date June, 2018: Split off from constructor into independent method
 *                   (M. Miesch, JCSDA)
 */
void AnalyticInit::fillGeoVaLs(const Locations & locs, GeoVaLs & geovals) const
{
  oops::Log::trace() << "AnalyticInit::analytic_init starting" << std::endl;
  if (config_.has("analytic_init")) {
      ufo_geovals_analytic_init_f90(geovals.toFortran(), locs.toFortran(), config_);
  }
  oops::Log::trace() << "AnalyticInit::analytic_init done" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace ufo
