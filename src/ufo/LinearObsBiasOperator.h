/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_LINEAROBSBIASOPERATOR_H_
#define UFO_LINEAROBSBIASOPERATOR_H_

#include <vector>

#include "oops/util/Printable.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsDiagnostics;

/// Tangent linear and adjoint of the linear combination bias correction operator
class LinearObsBiasOperator : public util::Printable {
 public:
  explicit LinearObsBiasOperator(ioda::ObsSpace &);

  /// Set trajectory (save predictors)
  void setTrajectory(const GeoVaLs &, const ObsBias &, ObsDiagnostics &);
  /// Compute TL of bias correction
  void computeObsBiasTL(const ObsBiasIncrement &, ioda::ObsVector &) const;
  /// Compute adjoint of bias correction
  void computeObsBiasAD(ObsBiasIncrement &, const ioda::ObsVector &) const;

 private:
  /// Print used for logging
  void print(std::ostream &) const override;

  /// ObsSpace used for this bias correction
  ioda::ObsSpace & odb_;

  /// do different records have separate bias-correction coefficients?
  bool byRecord_;

  /// predictors values; set in setTrajectory
  std::vector<ioda::ObsVector> predData_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_LINEAROBSBIASOPERATOR_H_
