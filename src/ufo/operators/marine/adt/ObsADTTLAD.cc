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
  std::vector<int> operatorVarIndices;
  getOperatorVariables(params.variables.value(), odb.assimvariables(),
    operatorVars_, operatorVarIndices)

  // sanity check to make sure adt is the only variable
  ASSERT_MSG(
    operatorVars_.size() == 1 && operatorVars_[0] == "absoluteDynamicTopography",
    "ADT can only work on variable \"absoluteDynamicTopography\"");
  ASSERT(operatorVarIndices.size() == 1);
  operatorVarIndex_ = operatorVarIndices[0];

  requiredVars_.push_back("sea_surface_height_above_geoid");
}

// -----------------------------------------------------------------------------

ObsADTTLAD::~ObsADTTLAD() {
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &) {
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  std::vector<double> vec(ovec.nlocs());
  const double missing = util::missingValue(missing);

  // get geoval
  const std::string varname = ovec.varnames().variables()[operatorVarIndex_];
  geovals.getAtLevel(vec, "sea_surface_height_above_geoid", 0);

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
  oops::Log::debug() << "ObsADT simulateObsTL offset: " << offset<<std::endl;

  // subtract offset from geoval
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    const size_t idx = jloc * ovec.nvars() + operatorVarIndex_;
    ovec[idx] = vec[jloc];
    if (ovec[idx] != missing) ovec[idx] -= offset;
  }
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  std::vector<double> vec(ovec.nlocs());
  const double missing = util::missingValue(missing);

  // get geoval
  const std::string varname = ovec.varnames().variables()[operatorVarIndex_];
  geovals.getAtLevel(vec, "sea_surface_height_above_geoid", 0);

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
  oops::Log::debug() << "ObsADT simulateObsAD offset: " << offset<<std::endl;

  // subtract offset from geoval
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    const size_t idx = jloc * ovec.nvars() + operatorVarIndex_;
    if (ovec[idx] != missing)
      vec[idx] += ovec[jloc] - offset;
  }
  geovals.putAtLevel(vec, "sea_surface_height_above_geoid", 0);
}

// -----------------------------------------------------------------------------

void ObsADTTLAD::print(std::ostream & os) const {
  os << "ObsADTTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
