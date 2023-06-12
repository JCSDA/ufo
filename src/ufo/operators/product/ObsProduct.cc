/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/product/ObsProduct.h"

#include <ostream>
#include <vector>

#include "ioda/ObsVector.h"

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
  : ObsOperatorBase(odb, VariableNameMap(parameters.AliasFile.value()))
{
  oops::Log::trace() << "ObsProduct constructor starting" << std::endl;

  getOperatorVariables(parameters.variables.value(), odb.assimvariables(),
                       operatorVars_, operatorVarIndices_);
  requiredVars_ += nameMap_.convertName(operatorVars_);

  // Add scaling variable to the list
  requiredVars_.push_back(parameters.geovalsToScaleHofxBy.value());

  oops::Log::trace() << "ObsProduct constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsProduct::~ObsProduct() {
  oops::Log::trace() << "ObsProduct destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProduct::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                              ObsDiagnostics &) const {
  oops::Log::trace() << "ObsProduct: simulateObs starting" << std::endl;

  // Scale h(x) by a certain GeoVaLs
  std::vector<double> scalingGeoVaLs(ovec.nlocs());
  gv.get(scalingGeoVaLs, requiredVars_.variables().back());

  // Set missing values to 1.0
  auto missing = util::missingValue(scalingGeoVaLs[0]);
  for (double& element : scalingGeoVaLs) {
    element = (element == missing) ? 1.0 : element;
  }

  std::vector<double> vec(ovec.nlocs());
  for (int jvar : operatorVarIndices_) {
    const std::string varname = nameMap_.convertName(ovec.varnames().variables()[jvar]);
    // Get GeoVaL at the level closest to the Earth's surface.
    gv.getAtLevel(vec, varname, gv.nlevs(varname) - 1);
    for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
      const size_t idx = jloc * ovec.nvars() + jvar;
      ovec[idx] = vec[jloc] * scalingGeoVaLs[jloc];
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
