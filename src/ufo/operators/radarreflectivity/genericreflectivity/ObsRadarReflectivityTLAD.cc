/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/radarreflectivity/genericreflectivity/ObsRadarReflectivityParameters.h"
#include "ufo/operators/radarreflectivity/genericreflectivity/ObsRadarReflectivityTLAD.h"

#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsRadarReflectivityTLAD>
makerRadarReflectivityTL_("RadarReflectivity");
// -----------------------------------------------------------------------------

ObsRadarReflectivityTLAD::ObsRadarReflectivityTLAD(const ioda::ObsSpace & odb,
                                                   const Parameters_ & params)
  : LinearObsOperatorBase(odb, VariableNameMap(params.AliasFile.value())),
    odb_(odb)
{
  oops::Log::trace() << "ObsRadarReflectivityTLAD constructor entered" << std::endl;

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
                                         requiredVarsPlaceholder_,
                                         requiredVarsTL_);

  oops::Log::trace() << "ObsRadarReflectivityTLAD constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadarReflectivityTLAD::~ObsRadarReflectivityTLAD() {
  oops::Log::trace() << "ObsRadarReflectivityTLAD destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarReflectivityTLAD::setTrajectory(const GeoVaLs & gv,
                                             ObsDiagnostics & diag,
                                             const QCFlags_t & flag) {
  oops::Log::trace() << "ObsRadarReflectivityTLAD: setTrajectory entered" << std::endl;
  reflectivityAlgorithm_->setTrajectory(gv, diag, flag);
  oops::Log::trace() << "ObsRadarReflectivityTLAD: setTrajectory finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarReflectivityTLAD::simulateObsTL(const GeoVaLs & dx,
                                             ioda::ObsVector & dy,
                                             const QCFlags_t & flag) const {
  oops::Log::trace() << "ObsRadarReflectivityTLAD: simulateObsTL entered" << std::endl;
  reflectivityAlgorithm_->simulateObsTL(dx, dy, flag);
  oops::Log::trace() << "ObsRadarReflectivityTLAD: simulateObsTL finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarReflectivityTLAD::simulateObsAD(GeoVaLs & dx,
                                             const ioda::ObsVector & dy,
                                             const QCFlags_t & flag) const {
  oops::Log::trace() << "ObsRadarReflectivityTLAD: simulateObsAD entered" << std::endl;
  reflectivityAlgorithm_->simulateObsAD(dx, dy, flag);
  oops::Log::trace() << "ObsRadarReflectivityTLAD: simulateObsAD finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarReflectivityTLAD::print(std::ostream & os) const {
  reflectivityAlgorithm_->print(os);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
