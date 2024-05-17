/*
 * (C) Copyright 2021 Met Office
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_LEGENDRE_H_
#define UFO_PREDICTORS_LEGENDRE_H_

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/predictors/PredictorBase.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// -----------------------------------------------------------------------------

/// Configuration parameters of the Legendre predictor.
class LegendreParameters : public PredictorParametersBase {
  OOPS_CONCRETE_PARAMETERS(LegendreParameters, PredictorParametersBase);

 public:
  /// Number of scan positions.
  oops::RequiredParameter<int> numScanPositions{"number of scan positions", this};

  /// Order of the Legendre polynomial. By default, 1.
  ///
  /// \note If this option is set, a suffix containing its value (even if it's equal to 1) will be
  /// appended to the predictor name.
  oops::OptionalParameter<size_t> order{"order", this};
};

// -----------------------------------------------------------------------------
/**
 *This Legendre predictor is for fitting residual errors between the ends of the
 *scan extremes; thus, the scan position is rescaled between -1 and 1 (xscan variable).
 *Legendre polynomials are calculated with recurrence relations. Polynomials can be 
 *produced for an order (as addressed by input parameter of that name). The data must 
 *contain scan positions in the variable "MetaData/sensorScanPosition" and the number of scan 
 *positons in"MetaData/sensorNumbScanPosition". Both of these are vectors with length
 *nlocs in the current implementation.
 */


class Legendre : public PredictorBase {
 public:
  /// The type of parameters accepted by the constructor of this predictor.
  /// This typedef is used by the PredictorFactory.
  typedef LegendreParameters Parameters_;

  Legendre(const Parameters_ &, const oops::ObsVariables &);
  ~Legendre() {}

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               const ObsBias &,
               ioda::ObsVector &) const override;

 private:
  int order_;
  int nscan_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_LEGENDRE_H_
