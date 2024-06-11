/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/adt/ObsADTTLAD.h"

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
  : LinearObsOperatorBase(odb, VariableNameMap(params.AliasFile.value())),
    odb_(odb)
{
  oops::Log::trace() << "ObsADTTLAD constructor starting" << std::endl;

  std::vector<int> operatorVarIndices;
  getOperatorVariables(params.variables.value(), odb.assimvariables(),
    operatorVars_, operatorVarIndices);
  requiredVars_.push_back("sea_surface_height_above_geoid");

  // sanity check to make sure adt is the only variable
  ASSERT_MSG(
    operatorVars_.size() == 1 && operatorVars_[0] == "absoluteDynamicTopography",
    "ADT can only work on variable \"absoluteDynamicTopography\"");
  ASSERT(operatorVarIndices.size() == 1);
  operatorVarIndex_ = operatorVarIndices[0];

  oops::Log::trace() << "ObsADTTLAD constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsADTTLAD::~ObsADTTLAD() {
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                               const QCFlags_t & qc_flags) {
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                               const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsADTTLAD: simulateObsTL starting" << std::endl;

  const double missing = util::missingValue<double>();

  // get geovals
  std::vector<double> vec(ovec.nlocs());
  geovals.getAtLevel(vec, oops::Variable{"sea_surface_height_above_geoid"}, 0);

  // calculate global offsets
  double offset = 0;
  auto accumVal = odb_.distribution()->createAccumulator<double>();
  auto accumCnt = odb_.distribution()->createAccumulator<int>();
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    if (vec[jloc] != missing) {
      accumVal->addTerm(jloc, vec[jloc]);
      accumCnt->addTerm(jloc, 1);
    }
  }
  auto count = accumCnt->computeResult();
  if (count > 0) {
    offset = accumVal->computeResult() / count;
  }
  oops::Log::debug() << "ObsADT simulateObsTL offset: " << offset << std::endl;

  // subtract offset from geoval
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    const size_t idx = jloc * ovec.nvars() + operatorVarIndex_;
    ovec[idx] = vec[jloc];
    if (ovec[idx] != missing) ovec[idx] -= offset;
  }

  oops::Log::trace() << "ObsADTTLAD: simulateObsTL finished" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                               const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsADTTLAD: simulateObsAD starting" << std::endl;

  const double missing = util::missingValue<double>();

  // get geovals
  std::vector<double> vec(ovec.nlocs());
  geovals.getAtLevel(vec, oops::Variable{"sea_surface_height_above_geoid"}, 0);

  // calculate global offsets
  double offset = 0;
  auto accumVal = odb_.distribution()->createAccumulator<double>();
  auto accumCnt = odb_.distribution()->createAccumulator<int>();
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    if (ovec[jloc] != missing) {
      accumVal->addTerm(jloc, ovec[jloc]);
      accumCnt->addTerm(jloc, 1);
    }
  }
  auto count = accumCnt->computeResult();
  if (count > 0) {
    offset = accumVal->computeResult() / count;
  }
  oops::Log::debug() << "ObsADT simulateObsAD offset: " << offset << std::endl;

  // subtract offset from geoval
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    const size_t idx = jloc * ovec.nvars() + operatorVarIndex_;
    if (ovec[idx] != missing)
      vec[idx] += ovec[jloc] - offset;
  }
  geovals.putAtLevel(vec, oops::Variable{"sea_surface_height_above_geoid"}, 0);

  oops::Log::trace() << "ObsADTTLAD: simulateObsAD finished" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::print(std::ostream & os) const {
  os << "ObsADTTLAD operator" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
