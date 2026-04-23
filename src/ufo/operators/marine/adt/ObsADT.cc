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
  oops::Log::trace() << "ObsADT constructor start" << std::endl;

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

  oops::Log::trace() << "ObsADT constructor start" << std::endl;
}

// -----------------------------------------------------------------------------

ObsADT::~ObsADT() {
  oops::Log::trace() << "ObsADT destructor start" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADT::simulateObs(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                         ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsADT::simulateObs start" << std::endl;

  const double missing = util::missingValue<double>();
  const size_t nlocs = ovec.nlocs();

  // -----------------------
  // 1. GeoVaLs
  // -----------------------
  std::vector<double> vec(nlocs);
  geovals.getAtLevel(vec,
                     oops::Variable{"sea_surface_height_above_geoid"}, 0);

  // -----------------------
  // 2. Observations
  // -----------------------
  std::vector<double> obs;
  odb_.get_db("ObsValue",
              ovec.varnames().variables()[operatorVarIndex_],
              obs);

  // -----------------------
  // 3. Compute global offset (QC-filtered)
  // -----------------------
  double offset = 0.0;

  auto accumVal = odb_.distribution()->createAccumulator<double>();
  auto accumCnt = odb_.distribution()->createAccumulator<int>();

  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    // ✅ CORRECT QC ACCESS FOR ObsDataVector<int>
    const int qc = qc_flags[operatorVarIndex_][jloc];

    if (qc == 0 &&
        obs[jloc] != missing &&
        vec[jloc] != missing &&
        std::isfinite(obs[jloc]) &&
        std::isfinite(vec[jloc])) {

      accumVal->addTerm(jloc, vec[jloc] - obs[jloc]);
      accumCnt->addTerm(jloc, 1);
    }
  }

  const int count = accumCnt->computeResult();

  if (count > 0) {
    offset = accumVal->computeResult() / count;
  }

  // -----------------------
  // 4. Apply offset (no QC filtering here!)
  // -----------------------
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    const size_t idx = jloc * ovec.nvars() + operatorVarIndex_;

    ovec[idx] = vec[jloc];

    if (ovec[idx] != missing && std::isfinite(ovec[idx])) {
      ovec[idx] -= offset;
    }
  }

  oops::Log::trace() << "ObsADT::simulateObs done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsADT::print(std::ostream & os) const {
  os << "ObsADT operator" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
