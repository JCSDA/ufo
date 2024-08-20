/*
 * (C) Copyright 2023- UCAR
 # (C) Crown Copyright 2024 - UK Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/logarithm/ObsLogarithm.h"

#include <algorithm>
#include <cassert>
#include <ostream>
#include <sstream>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsLogarithm> obsLogarithmMaker_("Logarithm");
// -----------------------------------------------------------------------------

ObsLogarithm::ObsLogarithm(const ioda::ObsSpace& odb,
                           const Parameters_& parameters)
    : ObsOperatorBase(odb, VariableNameMap(parameters.AliasFile.value())),
      odb_(odb) {
  oops::Log::trace() << "ObsLogarithm constructor starting" << std::endl;

  getOperatorVariables(parameters.variables.value(), odb.assimvariables(),
                       operatorVars_, operatorVarIndices_);
  requiredVars_ += nameMap_.convertName(operatorVars_);

  if (parameters.logBase.value() != boost::none) {
    // If parameters.logBase.value() is boost::none, logBase_ will be 0.0 (set
    // in the header) and the log will be taken to the base e.
    logBase_ = parameters.logBase.value().value();
    if (logBase_ <= 0.0 || logBase_ == 1.0) {
      std::ostringstream msg;
      msg << "Invalid log base: " << logBase_;
      throw eckit::BadValue(msg.str(), Here());
    }
  }
  oops::Log::trace() << "ObsLogarithm constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsLogarithm::~ObsLogarithm() {
  oops::Log::trace() << "ObsLogarithm destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsLogarithm::simulateObs(const GeoVaLs& gv, ioda::ObsVector& ovec,
                               ObsDiagnostics& odiags,
                               const QCFlags_t& qc_flags) const {
  oops::Log::trace() << "ObsLogarithm: simulateObs starting" << std::endl;

  assert(logBase_ >= 0.0);
  assert(logBase_ != 1.0);
  std::vector<double> vec(ovec.nlocs());
  const double missing = util::missingValue<double>();
  for (int jvar : operatorVarIndices_) {
    const oops::Variable var =
        nameMap_.convertName(ovec.varnames().variables()[jvar]);
    // Fill vec with the GeoVaL at the level closest to the Earth's surface.
    gv.getAtLevel(vec, var, gv.nlevs(var) - 1);
    for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
      const size_t idx = jloc * ovec.nvars() + jvar;
      if (vec[jloc] == missing || vec[jloc] <= 0.0) {
        ovec[idx] = missing;
      } else if (logBase_ == 0.0) {
        // natural log
        ovec[idx] = std::log(vec[jloc]);
      } else if (logBase_ == 2.0) {
        ovec[idx] = std::log2(vec[jloc]);
      } else if (logBase_ == 10.0) {
        ovec[idx] = std::log10(vec[jloc]);
      } else {
        ovec[idx] = std::log(vec[jloc]) / std::log(logBase_);
      }
    }
  }

  oops::Log::trace() << "ObsLogarithm: simulateObs finished" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsLogarithm::print(std::ostream& os) const {
  os << "ObsLogarithm operator" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
