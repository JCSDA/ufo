/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/adt/ObsADTTLAD.h"

#include <cmath>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/distribution/Accumulator.h"
#include "ioda/distribution/Distribution.h"
#include "ioda/ObsVector.h"

#include "ufo/filters/Variable.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/OperatorUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsADTTLAD> makerADTTL_("ADT");
// -----------------------------------------------------------------------------

ObsADTTLAD::ObsADTTLAD(const ioda::ObsSpace & odb, const Parameters_ & params)
  : LinearObsOperatorBase(odb, VariableNameMap(params.AliasFile.value(),
                                               params.variableMaps.value())),
    odb_(odb)
{
  oops::Log::trace() << "ObsADTTLAD constructor start" << std::endl;

  std::vector<int> operatorVarIndices;
  getOperatorVariables(params.variables.value(), odb.assimvariables(),
    operatorVars_, operatorVarIndices);
  requiredVars_.push_back("sea_surface_height_above_geoid");  // aka var_abs_topo

  // sanity check to make sure adt is the only variable
  ASSERT_MSG(
    operatorVars_.size() == 1 && operatorVars_[0] == "absoluteDynamicTopography",
    "ADT can only work on variable \"absoluteDynamicTopography\"");
  ASSERT(operatorVarIndices.size() == 1);
  operatorVarIndex_ = operatorVarIndices[0];

  oops::Log::trace() << "ObsADTTLAD constructor done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsADTTLAD::~ObsADTTLAD() {
  oops::Log::trace() << "ObsADTTLAD destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                               const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsADTTLAD::setTrajectory start" << std::endl;

  const double missing = util::missingValue<double>();
  const size_t nlocs = odb_.nlocs();

  // trajectory geoval (sea surface height above geoid)
  std::vector<double> vec(nlocs);
  geovals.getAtLevel(vec, oops::Variable{"sea_surface_height_above_geoid"}, 0);

  // observations
  std::vector<double> obs;
  odb_.get_db("ObsValue", odb_.assimvariables().variables()[operatorVarIndex_], obs);

  // Cache the locations that enter the global-mean offset: those that pass QC
  // and have finite obs and geoval. This must match the set used by the nonlinear
  // ObsADT::simulateObs so that the TL/AD are its exact linearization.
  qcmask_.assign(nlocs, 0);
  auto accumCnt = odb_.distribution()->createAccumulator<int>();
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    // Note: only using QC flags for passed observations, passive
    // observations are excluded from the global-mean offset.
    if (qc_flags[operatorVarIndex_][jloc] == 0 &&
        obs[jloc] != missing &&
        vec[jloc] != missing &&
        std::isfinite(obs[jloc]) &&
        std::isfinite(vec[jloc])) {
      qcmask_[jloc] = 1;
      accumCnt->addTerm(jloc, 1);
    }
  }
  qccount_ = accumCnt->computeResult();

  oops::Log::trace() << "ObsADTTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  oops::Log::trace() << "ObsADTTLAD::simulateObsTL start" << std::endl;

  const double missing = util::missingValue<double>();

  // get geovals (SSH increment)
  std::vector<double> vec(ovec.nlocs());
  geovals.getAtLevel(vec, oops::Variable{"sea_surface_height_above_geoid"}, 0);

  // Global offset = mean of the SSH increment over the QC-filtered set that
  // was cached at setTrajectory.
  double offset = 0;
  auto accumVal = odb_.distribution()->createAccumulator<double>();
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    if (qcmask_[jloc]) {
      accumVal->addTerm(jloc, vec[jloc]);
    }
  }
  if (qccount_ > 0) {
    offset = accumVal->computeResult() / qccount_;
  }

  // subtract offset from geoval (applied at every non-missing location)
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    const size_t idx = jloc * ovec.nvars() + operatorVarIndex_;
    ovec[idx] = vec[jloc];
    if (ovec[idx] != missing) ovec[idx] -= offset;
  }

  oops::Log::trace() << "ObsADTTLAD::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  oops::Log::trace() << "ObsADTTLAD::simulateObsAD start" << std::endl;

  const double missing = util::missingValue<double>();

  // get geovals (adjoint accumulator)
  std::vector<double> vec(ovec.nlocs());
  geovals.getAtLevel(vec, oops::Variable{"sea_surface_height_above_geoid"}, 0);

  auto accumVal = odb_.distribution()->createAccumulator<double>();
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    const size_t idx = jloc * ovec.nvars() + operatorVarIndex_;
    if (ovec[idx] != missing) {
      accumVal->addTerm(jloc, ovec[idx]);
    }
  }
  const double A = accumVal->computeResult();
  const double meanAdj = (qccount_ > 0) ? A / qccount_ : 0.0;

  // accumulate into geoval
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    const size_t idx = jloc * ovec.nvars() + operatorVarIndex_;
    if (ovec[idx] != missing) vec[jloc] += ovec[idx];
    if (qcmask_[jloc]) vec[jloc] -= meanAdj;
  }
  geovals.putAtLevel(vec, oops::Variable{"sea_surface_height_above_geoid"}, 0);

  oops::Log::trace() << "ObsADTTLAD::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::print(std::ostream & os) const {
  os << "ObsADTTLAD operator" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
