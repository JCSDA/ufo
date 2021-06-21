/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ERRORS_OBSERRORDIAGONAL_H_
#define UFO_ERRORS_OBSERRORDIAGONAL_H_

#include <sstream>
#include <string>

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"


namespace ufo {

/// \brief Parameters for diagonal obs errors
class ObsErrorDiagonalParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorDiagonalParameters, Parameters)
 public:
  /// perturbation amplitude multiplier
  oops::Parameter<double> pert{"random amplitude", 1.0, this};
};

// -----------------------------------------------------------------------------
/// \brief Diagonal observation error covariance matrix.
class ObsErrorDiagonal : public util::Printable {
 public:
  static const std::string classname() {return "ufo::ObsErrorDiagonal";}

  ObsErrorDiagonal(const eckit::Configuration &, ioda::ObsSpace &);

/// Update after obs errors potentially changed
  void update(const ioda::ObsVector &);

/// Multiply a Departure by \f$R\f$
  void multiply(ioda::ObsVector &) const;

/// Multiply a Departure by \f$R^{-1}\f$
  void inverseMultiply(ioda::ObsVector &) const;

/// Generate random perturbation
  void randomize(ioda::ObsVector &) const;

/// Save obs errors
  void save(const std::string &) const;

/// Get mean error for Jo table
  double getRMSE() const {return stddev_.rms();}

/// Get obs errors std deviation
  ioda::ObsVector obserrors() const;

/// Return inverseVariance
  ioda::ObsVector inverseVariance() const;

 private:
  void print(std::ostream &) const;
  ioda::ObsVector stddev_;
  ioda::ObsVector inverseVariance_;
  ObsErrorDiagonalParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERRORDIAGONAL_H_
