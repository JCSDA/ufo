/*
 * (C) Copyright 2024 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/radarreflectivity/genericreflectivity/ObsRadarReflectivity.h"
#include "ufo/operators/radarreflectivity/genericreflectivity/ObsRadarReflectivityParameters.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadarReflectivity> obsRadarReflectivityMaker_("RadarReflectivity");
// -----------------------------------------------------------------------------

ObsRadarReflectivity::ObsRadarReflectivity(const ioda::ObsSpace & odb,
                                           const Parameters_ & params)
  : ObsOperatorBase(odb), odb_(odb)
{
  oops::Log::trace() << "ObsRadarReflectivity constructor entered" << std::endl;

  // Enable this operator to be used as part of the Composite operator.
  // Indices of operator variables.
  std::vector<int> operatorVarIndices;
  // Get operator variables and their indices in the full list of simulated variables.
  getOperatorVariables(params.variables.value(), odb.assimvariables(),
                       operatorVars_, operatorVarIndices);
  // Only reflectivity may be used. Its index in the full list of simulated variables
  // is passed to the reflectivity algorithm constructor.
  if (operatorVars_ != oops::ObsVariables({"reflectivity"})) {
      throw eckit::UserError("This operator can only be used "
                             "to simulate reflectivity.",
                             Here());
    }

  // Instantiate reflectivity algorithm.
  const ReflectivityAlgorithmParametersWrapper & algparams =
    params.reflectivityAlgorithmParameters.value();
  reflectivityAlgorithm_ =
    ReflectivityAlgorithmFactory::create(algparams.reflectivityAlgorithmName,
                                         odb,
                                         operatorVarIndices[0],
                                         requiredVars_,
                                         requiredVarsTLPlaceholder_);

  oops::Log::trace() << "ObsRadarReflectivity constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadarReflectivity::~ObsRadarReflectivity() {
  oops::Log::trace() << "ObsRadarReflectivity destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarReflectivity::simulateObs(const GeoVaLs & gv,
                                       ioda::ObsVector & ovec,
                                       ObsDiagnostics & diag,
                                       const QCFlags_t & flag) const {
  oops::Log::trace() << "ObsRadarReflectivity: simulateObs entered" << std::endl;
  reflectivityAlgorithm_->simulateObs(gv, ovec, diag, flag);
  oops::Log::trace() << "ObsRadarReflectivity: simulateObs finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarReflectivity::print(std::ostream & os) const {
  reflectivityAlgorithm_->print(os);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
