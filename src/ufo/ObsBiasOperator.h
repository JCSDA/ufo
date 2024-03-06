/*
 * (C) Copyright 2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIASOPERATOR_H_
#define UFO_OBSBIASOPERATOR_H_
#include "ioda/ObsDataVector.h"
#include "oops/util/Printable.h"

// forward declarations

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsDiagnostics;

/// \brief Application of bias correction.
/// \details Bias correction is computed as a linear combination of bias
/// predictors using bias correction coefficients, both specified in ObsBias class.
class ObsBiasOperator : public util::Printable {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  explicit ObsBiasOperator(ioda::ObsSpace &);

  /// Compute bias correction
  void computeObsBias(const GeoVaLs &, ioda::ObsVector &, const ObsBias &,
                      ObsDiagnostics &, const QCFlags_t &) const;

 private:
  /// Print details (used for logging)
  void print(std::ostream &) const override;

  /// ObsSpace used for computing predictors
  ioda::ObsSpace & odb_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASOPERATOR_H_
