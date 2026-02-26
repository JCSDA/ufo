/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/densityreduction/ObsDensityReductionTLAD.h"

#include <ostream>
#include <utility>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/operators/densityreduction/ObsDensityReductionParameters.h"
#include "ufo/ScopedDefaultGeoVaLFormatChange.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsDensityReductionTLAD>
 makerDensityReductionTL_("Density Reduction");
// -----------------------------------------------------------------------------

ObsDensityReductionTLAD::ObsDensityReductionTLAD
(const ioda::ObsSpace & odb, const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), odb_(odb)
{
  oops::Log::trace() << "ObsDensityReductionTLAD constructor start" << std::endl;

  const eckit::LocalConfiguration &operatorConfig = parameters.oper.value();
  ObsOperatorParametersWrapper operatorParams;
  operatorParams.deserialize(operatorConfig);
  operator_ = std::unique_ptr<LinearObsOperatorBase>
    (LinearObsOperatorFactory::create(odb, operatorParams.operatorParameters));

  requiredVars_ += operator_->requiredVars();

  oops::Log::trace() << "ObsDensityReductionTLAD constructor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsDensityReductionTLAD::~ObsDensityReductionTLAD() {
  oops::Log::trace() << "ObsDensityReductionTLAD destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsDensityReductionTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics & ydiags,
                                     const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsDensityReductionTLAD::setTrajectory start" << std::endl;

  ScopedDefaultGeoVaLFormatChange change(geovals, GeoVaLFormat::REDUCED);
  operator_->setTrajectory(geovals, ydiags, qc_flags);
  ScopedDefaultGeoVaLFormatChange changeback(geovals, GeoVaLFormat::SAMPLED);

  oops::Log::trace() << "ObsDensityReductionTLAD::setTrajectory done" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsDensityReductionTLAD::simulateObsTL(const GeoVaLs & geovals,
                                            ioda::ObsVector & ovec) const {
  oops::Log::trace() << "ObsDensityReductionTLAD::simulateObsTL start" << std::endl;

  ScopedDefaultGeoVaLFormatChange change(geovals, GeoVaLFormat::REDUCED);
  operator_->simulateObsTL(geovals, ovec);
  ScopedDefaultGeoVaLFormatChange changeback(geovals, GeoVaLFormat::SAMPLED);

  oops::Log::trace() << "ObsDensityReductionTLAD::simulateObsTL done " <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsDensityReductionTLAD::simulateObsAD(GeoVaLs & geovals,
                                            const ioda::ObsVector & ovec) const {
  oops::Log::trace() << "ObsDensityReductionTLAD::simulateObsAD start" << std::endl;

  ScopedDefaultGeoVaLFormatChange change(geovals, GeoVaLFormat::REDUCED);
  operator_->simulateObsAD(geovals, ovec);
  ScopedDefaultGeoVaLFormatChange changeback(geovals, GeoVaLFormat::SAMPLED);

  oops::Log::trace() << "ObsDensityReductionTLAD::simulateObsAD done" <<  std::endl;
}

// -----------------------------------------------------------------------------

oops::ObsVariables ObsDensityReductionTLAD::simulatedVars() const {
  return operator_->simulatedVars();
}

// -----------------------------------------------------------------------------

void ObsDensityReductionTLAD::print(std::ostream & os) const {
  os << "ObsDensityReduction with operator:" << std::endl;
  os << *operator_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
