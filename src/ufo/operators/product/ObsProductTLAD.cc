/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/product/ObsProductTLAD.h"

#include <ostream>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsProductTLAD> makerProductTL_("Product");
// -----------------------------------------------------------------------------

ObsProductTLAD::ObsProductTLAD(const ioda::ObsSpace & odb,
                                 const Parameters_ & parameters)
  : LinearObsOperatorBase(odb, VariableNameMap(parameters.AliasFile.value()))
{
  oops::Log::trace() << "ObsProductTLAD constructor starting" << std::endl;

  getOperatorVariables(parameters.variables.value(), odb.assimvariables(),
                       operatorVars_, operatorVarIndices_);

  if (parameters.geovalVariable.value() != boost::none) {
      geovalName_ = parameters.geovalVariable.value().value();
      requiredVars_.push_back(parameters.geovalVariable.value().value());
      operatorVarIndices_.resize(1);
     } else {
      requiredVars_ += nameMap_.convertName(operatorVars_);
  }

  // Save scaling variable name
  scalingGeoVar_ = parameters.geovalsToScaleHofxBy.value();

  // Initialize the vector that will hold the scaling GeoVaLs with the correct length.
  scalingGeoVaLs_.resize(odb.nlocs());

  // Set geovals Exponent
  if (parameters.geovalsExponent.value() != boost::none) {
      geovalsExponent_ = parameters.geovalsExponent.value().value();
  }

  oops::Log::trace() << "ObsProductTLAD constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsProductTLAD::~ObsProductTLAD() {
  oops::Log::trace() << "ObsProductTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProductTLAD::setTrajectory(const GeoVaLs & gv, ObsDiagnostics &) {
  // Save the scaling geovals
  gv.get(scalingGeoVaLs_, scalingGeoVar_);

  // Set missing values to 1.0
  auto missing = util::missingValue(scalingGeoVaLs_[0]);
  for (double& element : scalingGeoVaLs_) {
    element = (element == missing) ? 1.0 : element;
  }

  // If exponent is set, do (geovals)^a
  if (geovalsExponent_ != 0) {
      const bool exponentIsFraction = static_cast<int>(geovalsExponent_) != geovalsExponent_;
      for (double& element : scalingGeoVaLs_) {
          if (element < 0 && exponentIsFraction) {
              oops::Log::warning() << "Trying to raise a negative number to non-integer exponent,"
                                      " '" << element << "' in scaling geovals set to missing"
                                   << std::endl;
              element = missing;
          } else if (geovalsExponent_ < 0 &&
                     oops::is_close_absolute(element, 0.0, 1e-10, 0, oops::TestVerbosity::SILENT)) {
              oops::Log::warning() << "Trying to divide by zero, '"
                                   << element << "' in scaling geovals set to missing"
                                   << std::endl;
              element = missing;
          } else {
              element = std::pow(element, geovalsExponent_);
          }
      }
  }
  oops::Log::trace() << "ObsProductTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProductTLAD::simulateObsTL(const GeoVaLs & dx, ioda::ObsVector & dy) const {
  oops::Log::trace() << "ObsProductTLAD: TL observation operator starting" << std::endl;

  std::vector<double> vec(dy.nlocs());
  const auto missing = util::missingValue(scalingGeoVaLs_[0]);
  for (int jvar : operatorVarIndices_) {
      const std::string varname = (geovalName_ == "")?
                        nameMap_.convertName(dy.varnames().variables()[jvar])
                        : geovalName_;
    // Fill dy with dx at the level closest to the Earth's surface.
    dx.getAtLevel(vec, varname, dx.nlevs(varname) - 1);
    for (size_t jloc = 0; jloc < dy.nlocs(); ++jloc) {
      const size_t idx = jloc * dy.nvars() + jvar;
      if (scalingGeoVaLs_[jloc] != missing) {
          dy[idx] = vec[jloc] * scalingGeoVaLs_[jloc];
      }
    }
  }

  oops::Log::trace() << "ObsProductTLAD: TL observation operator finished" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProductTLAD::simulateObsAD(GeoVaLs & dx, const ioda::ObsVector & dy) const {
  oops::Log::trace() << "ObsProductTLAD: adjoint observation operator starting" << std::endl;

  const double missing = util::missingValue(missing);

  std::vector<double> vec(dy.nlocs());
  for (int jvar : operatorVarIndices_) {
      const std::string varname = (geovalName_ == "")?
                        nameMap_.convertName(dy.varnames().variables()[jvar])
                        : geovalName_;
    // Get current value of dx at the level closest to the Earth's surface.
    dx.getAtLevel(vec, varname, dx.nlevs(varname) - 1);
    // Increment dx with non-missing values of dy.
    for (size_t jloc = 0; jloc < dy.nlocs(); ++jloc) {
      const size_t idx = jloc * dy.nvars() + jvar;
      if (dy[idx] != missing && scalingGeoVaLs_[jloc] != missing)
        vec[jloc] += dy[idx] * scalingGeoVaLs_[jloc];
    }
    // Store new value of dx.
    dx.putAtLevel(vec, varname, dx.nlevs(varname) - 1);
  }

  oops::Log::trace() << "ObsProductTLAD: adjoint observation operator finished" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProductTLAD::print(std::ostream & os) const {
  os << "ObsProductTLAD operator" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
