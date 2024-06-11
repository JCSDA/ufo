/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/product/ObsProduct.h"

#include <algorithm>
#include <ostream>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsProduct> obsProductMaker_("Product");
// -----------------------------------------------------------------------------

ObsProduct::ObsProduct(const ioda::ObsSpace & odb,
                         const Parameters_ & parameters)
  : ObsOperatorBase(odb, VariableNameMap(parameters.AliasFile.value())), odb_(odb)
{
  oops::Log::trace() << "ObsProduct constructor starting" << std::endl;

  getOperatorVariables(parameters.variables.value(), odb.assimvariables(),
                       operatorVars_, operatorVarIndices_);

  if (parameters.geovalVariable.value() != boost::none) {
      // Check that there is only one simulated variable. When specifying geovalVariable,
      // only one simulated variable is currently supported, this could be extended in future.
      const bool initialized = parameters.variables.value().is_initialized();
      if ((initialized && (parameters.variables.value().value().size() > 1))
         || (!initialized && odb.assimvariables().size() > 1)) {
          throw eckit::UserError(
            "More than one simulated variable is not currently supported",
            Here());
      }
      geovalName_ = parameters.geovalVariable.value().value();
      requiredVars_.push_back(parameters.geovalVariable.value().value());
      operatorVarIndices_.resize(1);
     } else {
      requiredVars_ += nameMap_.convertName(operatorVars_);
  }

  // Add scaling variable to the list
  variableNameToScaleHofxBy_ = parameters.variableNameToScaleHofxBy.value();
  variableGroupToScaleHofxBy_ = parameters.variableGroupToScaleHofxBy.value();
  if (variableGroupToScaleHofxBy_ == "GeoVaLs") {
    // If the group is GeoVaLs (default) then push back to requiredVars
    requiredVars_.push_back(variableNameToScaleHofxBy_);
  }

  // Set scaling variable exponent
  if (parameters.scalingVariableExponent.value() != boost::none) {
      scalingVariableExponent_ = parameters.scalingVariableExponent.value().value();
  }

  oops::Log::trace() << "ObsProduct constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsProduct::~ObsProduct() {
  oops::Log::trace() << "ObsProduct destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProduct::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                              ObsDiagnostics & odiags,
                              const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsProduct: simulateObs starting" << std::endl;

  // Get variable that will scale h(x)
  std::vector<double> scalingVariable(ovec.nlocs());
  if (variableGroupToScaleHofxBy_ == "GeoVaLs") {
    gv.get(scalingVariable, oops::Variable{variableNameToScaleHofxBy_});
  } else {
    // Get from the observation space
    odb_.get_db(variableGroupToScaleHofxBy_, variableNameToScaleHofxBy_, scalingVariable);
  }

  // Set missing values to 1.0
  const double missing = util::missingValue<double>();
  for (double& element : scalingVariable) {
    element = (element == missing) ? 1.0 : element;
  }

  // If exponent is set, do (scalingVariable)^a
  if (scalingVariableExponent_ != 0) {
      const bool exponentIsFraction = static_cast<int>(scalingVariableExponent_) !=
                                      scalingVariableExponent_;
      for (double& element : scalingVariable) {
          if (element < 0 && exponentIsFraction) {
              oops::Log::warning() << "Trying to raise a negative number to non-integer exponent,"
                                      " '" << element << "' in scaling geovals set to missing"
                                   << std::endl;
              element = missing;
          } else if (scalingVariableExponent_ < 0 &&
                     oops::is_close_absolute(element, 0.0, 1e-10, 0, oops::TestVerbosity::SILENT)) {
              oops::Log::warning() << "Trying to divide by zero, '"
                                   << element << "' in scaling geovals set to missing"
                                   << std::endl;
              element = missing;
          } else {
              element = std::pow(element, scalingVariableExponent_);
          }
      }
  }

  std::vector<double> vec(ovec.nlocs());
  for (int jvar : operatorVarIndices_) {
    const oops::Variable var = (geovalName_ == "")?
                      nameMap_.convertName(ovec.varnames().variables()[jvar])
                      : oops::Variable{geovalName_};
    // Get GeoVaL at the level closest to the Earth's surface.
    gv.getAtLevel(vec, var, gv.nlevs(var) - 1);
    for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
      const size_t idx = jloc * ovec.nvars() + jvar;
      if (scalingVariable[jloc] != missing) {
          ovec[idx] = vec[jloc] * scalingVariable[jloc];
      } else {
          ovec[idx] = missing;
      }
    }
  }

  oops::Log::trace() << "ObsProduct: simulateObs finished" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProduct::print(std::ostream & os) const {
  os << "ObsProduct operator" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
