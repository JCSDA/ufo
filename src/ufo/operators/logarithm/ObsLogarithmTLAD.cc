/*
 * (C) Copyright 2023- UCAR
 * (C) Crown Copyright 2024 - UK Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/logarithm/ObsLogarithmTLAD.h"

#include <cassert>
#include <ostream>
#include <sstream>
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
static LinearObsOperatorMaker<ObsLogarithmTLAD> makerLogarithmTL_("Logarithm");
// -----------------------------------------------------------------------------

ObsLogarithmTLAD::ObsLogarithmTLAD(const ioda::ObsSpace& odb,
                                   const Parameters_& parameters)
    : LinearObsOperatorBase(odb, VariableNameMap(parameters.AliasFile.value())),
      odb_(odb) {
  oops::Log::trace() << "ObsLogarithmTLAD constructor starting" << std::endl;

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

  oops::Log::trace() << "ObsLogarithmTLAD constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsLogarithmTLAD::~ObsLogarithmTLAD() {
  oops::Log::trace() << "ObsLogarithmTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsLogarithmTLAD::setTrajectory(const GeoVaLs& gv, ObsDiagnostics&,
                                     const QCFlags_t&) {
  oops::Log::trace() << "ObsLogarithmTLAD::setTrajectory starting" << std::endl;

  assert(logBase_ >= 0.0);
  assert(logBase_ != 1.0);
  const size_t nlocs = odb_.nlocs();
  const size_t nreqvars = requiredVars_.size();
  assert(nreqvars == operatorVars_.size());
  dHdx_.resize(nlocs * nreqvars);
  const double missing = util::missingValue<double>();

  std::vector<double> vec(nlocs);
  for (size_t jreqvar = 0; jreqvar < nreqvars; ++jreqvar) {
    // Fill vec with the GeoVaL at the level closest to the Earth's surface -
    // this is the value x about which we linearise.
    gv.getAtLevel(vec, requiredVars_[jreqvar],
                  gv.nlevs(requiredVars_[jreqvar]) - 1);
    // Now save dHdx_ for each x
    for (size_t jloc = 0; jloc < nlocs; ++jloc) {
      const size_t idx = jloc * nreqvars + jreqvar;  // unique index for each x
      if (vec[jloc] == missing || vec[jloc] <= 0.0) {
        dHdx_[idx] = missing;
      } else if (logBase_ == 0.0) {
        // Special case for natural logarithm (logBase_ == 0.0)
        dHdx_[idx] = 1.0 / vec[jloc];
      } else {
        dHdx_[idx] = 1.0 / (vec[jloc] * std::log(logBase_));
      }
    }
  }
  oops::Log::trace() << "ObsLogarithmTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsLogarithmTLAD::simulateObsTL(const GeoVaLs& dx, ioda::ObsVector& dy,
                                     const QCFlags_t& qc_flags) const {
  oops::Log::trace() << "ObsLogarithmTLAD: TL observation operator starting"
                     << std::endl;

  assert(dy.nlocs() == odb_.nlocs());
  assert(logBase_ >= 0.0);
  assert(logBase_ != 1.0);
  const size_t nlocs = odb_.nlocs();
  const size_t nreqvars = requiredVars_.size();
  assert(nreqvars == operatorVars_.size());
  assert(nreqvars == operatorVarIndices_.size());
  const double missing = util::missingValue<double>();

  // Set sections of dy which have the requiredVars_ to an associated value of
  // dHdx_ * dx, where dx is at the level closest to the Earth's surface.
  std::vector<double> vec(nlocs);
  for (size_t jreqvar = 0; jreqvar < nreqvars; ++jreqvar) {
    // Fill vec with the GeoVaL dx at the level closest to the Earth's surface
    dx.getAtLevel(vec, requiredVars_[jreqvar],
                  dx.nlevs(requiredVars_[jreqvar]) - 1);
    const size_t jvar = operatorVarIndices_[jreqvar];  // for indexing dy
    for (size_t jloc = 0; jloc < nlocs; ++jloc) {
      const size_t idx_dy =
          jloc * dy.nvars() + jvar;  // dy.nvars may not equal nreqvars
      const size_t idx_dHdx = jloc * nreqvars + jreqvar;
      if (vec[jloc] == missing || dHdx_[idx_dHdx] == missing) {
        dy[idx_dy] = missing;
      } else {
        dy[idx_dy] = dHdx_[idx_dHdx] * vec[jloc];
      }
    }
  }

  oops::Log::trace() << "ObsLogarithmTLAD: TL observation operator finished"
                     << std::endl;
}

// -----------------------------------------------------------------------------

void ObsLogarithmTLAD::simulateObsAD(GeoVaLs& dx, const ioda::ObsVector& dy,
                                     const QCFlags_t& qc_flags) const {
  oops::Log::trace()
      << "ObsLogarithmTLAD: adjoint observation operator starting" << std::endl;

  assert(dy.nlocs() == odb_.nlocs());
  assert(logBase_ >= 0.0);
  assert(logBase_ != 1.0);
  const size_t nlocs = odb_.nlocs();
  const size_t nreqvars = requiredVars_.size();
  assert(nreqvars == operatorVars_.size());
  assert(nreqvars == operatorVarIndices_.size());
  const double missing = util::missingValue<double>();

  // Increment dx at the level closest to the Earths surface with dHdx_ * dy,
  // skipping over missing values. Note again the dy has a separate index which
  // coverts just the requiredVars_.
  std::vector<double> vec(nlocs);
  for (size_t jreqvar = 0; jreqvar < nreqvars; ++jreqvar) {
    // Fill vec with the GeoVaL dx at the level closest to the Earth's surface
    dx.getAtLevel(vec, requiredVars_[jreqvar],
                  dx.nlevs(requiredVars_[jreqvar]) - 1);
    const size_t jvar = operatorVarIndices_[jreqvar];  // for indexing dy
    for (size_t jloc = 0; jloc < nlocs; ++jloc) {
      const size_t idx_dy =
          jloc * dy.nvars() + jvar;  // dy.nvars may not equal nreqvars
      const size_t idx_dHdx = jloc * nreqvars + jreqvar;
      if (dy[idx_dy] != missing && dHdx_[idx_dHdx] != missing) {
        vec[jloc] += dHdx_[idx_dHdx] * dy[idx_dy];
      }
    }
    // Store the new dx
    dx.putAtLevel(vec, requiredVars_[jreqvar],
                  dx.nlevs(requiredVars_[jreqvar]) - 1);
  }
  oops::Log::trace()
      << "ObsLogarithmTLAD: adjoint observation operator finished" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsLogarithmTLAD::print(std::ostream& os) const {
  os << "ObsLogarithmTLAD operator" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
