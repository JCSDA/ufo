/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/adt/ObsADT.h"

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"

#include "ioda/distribution/Accumulator.h"
#include "ioda/distribution/Distribution.h"
#include "ioda/ObsVector.h"

#include "ufo/filters/Variable.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/OperatorUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsADT> makerADT_("ADT");
// -----------------------------------------------------------------------------

ObsADT::ObsADT(const ioda::ObsSpace & odb, const ObsADTParameters & params)
  : ObsOperatorBase(odb, VariableNameMap(params.AliasFile.value())),
    odb_(odb)
{
  oops::Log::trace() << "ObsADT constructor starting" << std::endl;

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

  oops::Log::trace() << "ObsADT constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsADT::~ObsADT() {
}

// -----------------------------------------------------------------------------

void ObsADT::simulateObs(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                         ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsADT: simulateObs starting" << std::endl;

  const double missing = util::missingValue<double>();

  // get geovals
  std::vector<double> vec(ovec.nlocs());
  geovals.getAtLevel(vec, oops::Variable{"sea_surface_height_above_geoid"}, 0);

  // get obs
  std::vector<double> obs;
  odb_.get_db("ObsValue", ovec.varnames().variables()[operatorVarIndex_], obs);

  // calculate global offsets
  double offset = 0;
  auto accumVal = odb_.distribution()->createAccumulator<double>();
  auto accumCnt = odb_.distribution()->createAccumulator<int>();
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    // TODO(travis) also remove obs that have failed QC before this point?
    // Not doing it now to keep answers from changing with this PR.
    if (obs[jloc] != missing && vec[jloc] != missing) {
      accumVal->addTerm(jloc, vec[jloc] - obs[jloc]);
      accumCnt->addTerm(jloc, 1);
    }
  }
  auto count = accumCnt->computeResult();
  if (count > 0) {
    offset = accumVal->computeResult() / count;
  }
  oops::Log::debug() << "ObsADT: simulateObs offset: " << offset << std::endl;

  // subtract offset from geoval
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    const size_t idx = jloc * ovec.nvars() + operatorVarIndex_;
    ovec[idx] = vec[jloc];
    if (ovec[idx] != missing) ovec[idx] -= offset;
  }

  oops::Log::trace() << "ObsADT: simulateObs finished" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADT::print(std::ostream & os) const {
  os << "ObsADT operator" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
