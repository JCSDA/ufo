/*
 * (C) Copyright 2017-2018 UCAR
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
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsADT> makerADT_("ADT");
// -----------------------------------------------------------------------------

ObsADT::ObsADT(const ioda::ObsSpace & odb, const ObsADTParameters & params)
  : ObsOperatorBase(odb, VariableNameMap(params.AliasFile.value())),
    odb_(odb)
{
  getOperatorVariables(params.variables.value(), odb.assimvariables(),
    operatorVars_, operatorVarIndices_);

  // sanity check to make sure adt is the only variable
  ASSERT_MSG(
    operatorVars_.size() == 1 && operatorVars_[0] == "absoluteDynamicTopography",
    "ADT can only work on variable \"absoluteDynamicTopography\"");

  requiredVars_.push_back("sea_surface_height_above_geoid");
}

// -----------------------------------------------------------------------------

ObsADT::~ObsADT() {
}

// -----------------------------------------------------------------------------

void ObsADT::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                         ObsDiagnostics &) const {
  oops::Log::trace() << "ObsADT:simulateObs starting" << std::endl;

  std::vector<double> vec(ovec.nlocs());
  const double missing = util::missingValue(missing);

  for (int jvar : operatorVarIndices_) {

    const std::string varname = ovec.varnames().variables()[jvar];
    gv.getAtLevel(vec, "sea_surface_height_above_geoid", 0);
    std::vector<double> obs;
    odb_.get_db("ObsValue", varname, obs);

    // calculate global offsets
    double offset = 0;
    auto accumVal = odb_.distribution()->createAccumulator<double>();
    auto accumCnt = odb_.distribution()->createAccumulator<int>();
    for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
      // TODO(travis) also remove obs that have failed QC before this point?
      // Not doing it now to keep answers from changing with this PR.
      if (obs[jloc] != missing && vec[jloc] != missing) {
        accumVal->addTerm(jloc, obs[jloc] - vec[jloc]);
        accumCnt->addTerm(jloc, 1);
      }
    }
    auto count = accumCnt->computeResult();
    if (count > 0) {
      offset = accumVal->computeResult() / count;
    }

    // subtract offset from geoval
    for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
      const size_t idx = jloc * ovec.nvars() + jvar;
      ovec[idx] = vec[jloc] + offset;
    }
  }
}

// -----------------------------------------------------------------------------

void ObsADT::print(std::ostream & os) const {
  os << "ObsADT::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
