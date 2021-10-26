/*
 * (C) Copyright 2019 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/timeoper/ObsTimeOperTLAD.h"

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
#include "ufo/timeoper/ObsTimeOperUtil.h"

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
  oops::Log::trace() << "ObsTimeOperTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsTimeOperTLAD::~ObsTimeOperTLAD() {
  oops::Log::trace() << "ObsTimeOperTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTimeOperTLAD::setTrajectory(const GeoVaLs & geovals,
                                    ObsDiagnostics & ydiags) {
  oops::Log::trace() << "ObsTimeOperTLAD::setTrajectory entering" << std::endl;

  oops::Log::debug() << "ObsTimeOperTLAD::setTrajectory input geovals "
                     << geovals << std::endl;

  GeoVaLs gv1(obsspace().distribution(), geovals.getVars());
  GeoVaLs gv2(obsspace().distribution(), geovals.getVars());
  geovals.split(gv1, gv2);

  oops::Log::debug() << "ObsTimeOperTLAD::setTrajectory split geovals gv1 "
                     << gv1 << std::endl;

  oops::Log::debug() << "ObsTimeOperTLAD::setTrajectory split geovals gv2 "
                     << gv2 << std::endl;

  gv1 *= timeWeights_[0];
  gv2 *= timeWeights_[1];
  gv1 += gv2;

  oops::Log::debug() << "ObsTimeOperTLAD::setTrajectory final geovals gv1 "
                     << gv1 << std::endl;

  actualoperator_->setTrajectory(gv1, ydiags);

  oops::Log::debug() << gv1;

  oops::Log::trace() << "ObsTimeOperTLAD::setTrajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTimeOperTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  oops::Log::trace() << "ObsTimeOperTLAD::simulateObsTL entering" << std::endl;

  oops::Log::debug() << "ObsTimeOperTLAD::setTrajectory input geovals "
                     << geovals << std::endl;

  GeoVaLs gv1(obsspace().distribution(), geovals.getVars());
  GeoVaLs gv2(obsspace().distribution(), geovals.getVars());
  geovals.split(gv1, gv2);

  oops::Log::debug() << "ObsTimeOperTLAD::simulateObsTL split geovals gv1 "
                     << gv1 << std::endl;

  oops::Log::debug() << "ObsTimeOperTLAD::simulateObsTL split geovals gv2 "
                     << gv2 << std::endl;

  gv1 *= timeWeights_[0];
  gv2 *= timeWeights_[1];
  gv1 += gv2;

  oops::Log::debug() << "ObsTimeOperTLAD::simulateObsTL final geovals gv1 "
                     << gv1 << std::endl;

  actualoperator_->simulateObsTL(gv1, ovec);

  oops::Log::trace() << "ObsTimeOperTLAD::simulateObsTL exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTimeOperTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  oops::Log::trace() << "ObsTimeOperTLAD::simulateObsAD entering" << std::endl;

  oops::Log::debug() << "ObsTimeOperTLAD::simulateObsAD input geovals "
                     << geovals << std::endl;

  GeoVaLs gv1(obsspace().distribution(), geovals.getVars());
  GeoVaLs gv2(obsspace().distribution(), geovals.getVars());
  geovals.split(gv1, gv2);

  oops::Log::debug() << "ObsTimeOperTLAD::simulateObsAD split geovals gv1 "
                     << gv1 << std::endl;

  oops::Log::debug() << "ObsTimeOperTLAD::simulateObsAD split geovals gv2 "
                     << gv2 << std::endl;

  actualoperator_->simulateObsAD(gv1, ovec);

  gv2 = gv1;
  gv1 *= timeWeights_[0];
  gv2 *= timeWeights_[1];

  geovals.merge(gv1, gv2);

  oops::Log::debug() << "ObsTimeOperTLAD::simulateObsAD final geovals "
                     << geovals << std::endl;

  oops::Log::trace() << "ObsTimeOperTLAD::simulateObsAD exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTimeOperTLAD::print(std::ostream & os) const {
  os << "ObsTimeOperTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
