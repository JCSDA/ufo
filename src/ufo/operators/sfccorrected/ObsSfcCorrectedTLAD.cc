/*
 * (C) Crown Copyright 2024, Met Office
 * (C) Copyright 2024 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <memory>
#include <ostream>
#include <utility>
#include <vector>

#include "ufo/operators/sfccorrected/ObsSfcCorrectedTLAD.h"

#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/filters/Variables.h"
#include "ufo/GeoVaLs.h"
#include "ufo/operators/sfccorrected/EvalSurfacePressure.h"
#include "ufo/operators/sfccorrected/EvalSurfaceTemperature.h"
#include "ufo/utils/OperatorUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsSfcCorrectedTLAD> makerSfcCorrected_("SfcCorrected");
// -----------------------------------------------------------------------------

ObsSfcCorrectedTLAD::ObsSfcCorrectedTLAD(const ioda::ObsSpace & odb,
                                         const Parameters_ & params)
  : LinearObsOperatorBase(odb), requiredVars_(), odb_(odb), params_(params)
{
  oops::Log::trace() << "ObsSfcCorrectedTLAD constructor start." << std::endl;

  // Get the variables to simulate hofx for will be a subset of the assimilated variables
  getOperatorVariables(params_.variables.value(), odb.assimvariables(),
                       operatorVars_, operatorVarIndices_);

  // Create one surface operator per variable for the configured correction scheme
  operators_ = SurfaceOperatorFactory::createForVariables(operatorVars_, params_, requiredVars_);

  oops::Log::trace() << "ObsSfcCorrectedTLAD constructor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsSfcCorrectedTLAD::~ObsSfcCorrectedTLAD() {
  oops::Log::trace() << "ObsSfcCorrectedTLAD destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSfcCorrectedTLAD::setTrajectory(const GeoVaLs & gv,
                                    ObsDiagnostics &,
                                    const QCFlags_t &) {
  oops::Log::trace() << "ObsSfcCorrectedTLAD::setTrajectory start." << std::endl;
  for (size_t i = 0; i < operatorVars_.size(); ++i) {
    // Set Trajectory
    operators_[i]->settraj(gv, odb_);
  }
  oops::Log::trace() << "ObsSfcCorrectedTLAD::setTrajectory done." << std::endl;
}
// -----------------------------------------------------------------------------

void ObsSfcCorrectedTLAD::simulateObsTL(const GeoVaLs & gv, ioda::ObsVector & ovec) const {
  oops::Log::trace() << "ObsSfcCorrectedTLAD::simulateObsTL start." << std::endl;
  // Look over variables to calculate hofx
  std::vector<double> hofxTL(ovec.nlocs());
  for (size_t i = 0; i < operatorVars_.size(); ++i) {
     int jvar = operatorVarIndices_[i];
    // Calculate hofxTL
    operators_[i]->TL(gv, hofxTL);

    // Populate the observation vector
    for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
      const size_t idx = jloc * ovec.nvars() + jvar;
      ovec[idx] = hofxTL[jloc];
    }
  }
  oops::Log::trace() << "ObsSfcCorrectedTLAD::simulateObsTL done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSfcCorrectedTLAD::simulateObsAD(GeoVaLs & gv, const ioda::ObsVector & ovec) const {
  oops::Log::trace() << "ObsSfcCorrectedTLAD::simulateObsAD start." << std::endl;
  // Look over variables to calculate hofx
  std::vector<double> hofxAD(ovec.nlocs());
  for (size_t i = 0; i < operatorVars_.size(); ++i) {
     int jvar = operatorVarIndices_[i];

    // Populate the observation vector
    for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
      const size_t idx = jloc * ovec.nvars() + jvar;
      hofxAD[jloc] = ovec[idx];
    }
    operators_[i]->AD(gv, hofxAD);
  }
  oops::Log::trace() << "ObsSfcCorrectedTLAD::simulateObsAD done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSfcCorrectedTLAD::print(std::ostream & os) const {
  os << "ObsSfcCorrectedTLAD::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
