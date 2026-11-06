/*
 * (C) Copyright 2019 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/timeoper/ObsTimeOperTLAD.h"

#include <algorithm>
#include <ostream>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/operators/timeoper/ObsTimeOperUtil.h"
#include "ufo/ScopedDefaultGeoVaLFormatChange.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsTimeOperTLAD> makerTimeOperTL_("TimeOperLinInterp");
// -----------------------------------------------------------------------------

ObsTimeOperTLAD::ObsTimeOperTLAD(const ioda::ObsSpace & odb,
                                 const Parameters_ & parameters)
  : LinearObsOperatorBase(odb),
    actualoperator_(LinearObsOperatorFactory::create(
                      odb,
                      oops::validateAndDeserialize<LinearObsOperatorParametersWrapper>(
                        parameters.obsOperator.value()).operatorParameters)),
    timeWeights_(timeWeightCreate(odb, parameters))
{
  oops::Log::trace() << "ObsTimeOperTLAD constructor done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsTimeOperTLAD::~ObsTimeOperTLAD() {
  oops::Log::trace() << "ObsTimeOperTLAD destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTimeOperTLAD::setTrajectory(const GeoVaLs & geovals,
                                    ObsDiagnostics & ydiags,
                                    const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsTimeOperTLAD::setTrajectory start" << std::endl;

  // oops::Log::debug() << "ObsTimeOperTLAD::setTrajectory input geovals "
  //                    << geovals << std::endl;

  GeoVaLs gvt(geovals);
  timeInterpolate(gvt, timeWeights_);
  ScopedDefaultGeoVaLFormatChange change(gvt, GeoVaLFormat::REDUCED);
  actualoperator_->setTrajectory(gvt, ydiags, qc_flags);

  oops::Log::trace() << "ObsTimeOperTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTimeOperTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  oops::Log::trace() << "ObsTimeOperTLAD::simulateObsTL start" << std::endl;

  // oops::Log::debug() << "ObsTimeOperTLAD::setTrajectory input geovals "
  //                    << geovals << std::endl;

  GeoVaLs gvt(geovals);
  timeInterpolate(gvt, timeWeights_);
  ScopedDefaultGeoVaLFormatChange change(gvt, GeoVaLFormat::REDUCED);
  actualoperator_->simulateObsTL(gvt, ovec);

  oops::Log::trace() << "ObsTimeOperTLAD::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTimeOperTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  oops::Log::trace() << "ObsTimeOperTLAD::simulateObsAD start" << std::endl;

  // oops::Log::debug() << "ObsTimeOperTLAD::simulateObsAD input geovals "
  //                    << geovals << std::endl;

  GeoVaLs gvad(geovals);
  gvad.zero();
  timeInterpolate(gvad, timeWeights_);  // allocates the reduced format
  {
    ScopedDefaultGeoVaLFormatChange change(gvad, GeoVaLFormat::REDUCED);
    actualoperator_->simulateObsAD(gvad, ovec);
  }
  timeInterpolateAD(geovals, gvad, timeWeights_);

  // oops::Log::debug() << "ObsTimeOperTLAD::simulateObsAD final geovals "
  //                    << geovals << std::endl;

  oops::Log::trace() << "ObsTimeOperTLAD::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTimeOperTLAD::print(std::ostream & os) const {
  os << "ObsTimeOperTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
