/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_LEGENDRE_H_
#define UFO_PREDICTORS_LEGENDRE_H_

#include <vector>
#include "ufo/predictors/PredictorBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// -----------------------------------------------------------------------------
/**
 *This Legendre predictor is for fitting residual errors between the ends of the
 *scan extremes; thus, the scan position is rescaled between -1 and 1 (xscan variable).
 *Legendre polynomials are calculated with recurrence relations. Polynomials can be 
 *produced for an order (as addressed by input parameter of that name). The data must 
 *contain scan positions in the variable "scan_position@MetaData" and the number of scan 
 *positons in"sensor_numb_scan_position@MetaData". Both of these are vectors with length
 *nlocs in the current implementation.
 */


class Legendre : public PredictorBase {
 public:
  Legendre(const eckit::Configuration &, const oops::Variables &);
  ~Legendre() {}

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               ioda::ObsVector &) const override;

 private:
  int order_;
  int nscan_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_LEGENDRE_H_
