/*
 * (C) Copyright 2021- UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_LIGHTNING_OBSLIGHTNING_H_
#define UFO_OPERATORS_LIGHTNING_OBSLIGHTNING_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "ufo/ObsOperatorParametersBase.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/lightning/ObsLightning.interface.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// Configuration options recognized by the lightning operator intended for flash extent density.

class ObsLightningParameters: public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsLightningParameters, ObsOperatorParametersBase)
 public:
  //  oops::Parameter<bool> use_nonlinear
  //    {"use nonlinear operator",
  //     "use nonlinear operator (default is false)",
  //     false,
  //     this};
  oops::Parameter<size_t> num_gridpoints
    {"number of gridpoints",
     "number of gridpoints in one direction representing 15 km",
     5,
     this};
};

/*
 *  This operator was provided by Rong Kong (OU-CAPS) to simulate flash extent density matching
 *  observations from the GOES-R Geostationary Lightning Mapper (GLM). While foundational
 *  research and previous developments are well documented in the citations below, the current
 *  implementation in this code has been newly fitted using data collected from June 11 to June 15, 
 *  2023. This fitting process utilized GSL's RFFS 1-3 hour forecasts, which employ the Thompson
 *  microphysics scheme, provided hourly from 00 to 23 Z. This ensures that the operator is 
 *  fine-tuned and optimized based on the latest model data, incorporating adjustments for 
 *  potential systematic biases.
 *
 *  Relevant publications include:
 *     Kong, R., M. Xue, A. O. Fierro, Y. Jung, C. Liu, E. R. Mansell, and D. R. MacGorman, 2020:
 *          Assimilation of GOES-R Geostationary Lightning Mapper Flash Extent Density Data in
 *          GSI EnKF for the Analysis and Short-Term Forecast of a Mesoscale Convective System.
 *          Mon. Wea. Rev., 148, 2111-2133.
 *     Kong, R., M. Xue, C. Liu, A. O. Fierro, and E. R. Mansell, 2022: Development of New
 *          Observation Operators for Assimilating GOES-R Geostationary Lightning Mapper Flash
 *          Extent Density Data Using GSI EnKF: Tests with Two Convective Events over the United
 *          States. Mon. Wea. Rev., 150, 2091-2110.
 *
 *  The option named 'number of gridpoints' can be used to specify the number of gridpoints, N,
 *  representing a 15 km horizontal distance, and a square of N x N points is used to perform the
 *  integration of graupel mass column needed by the operator.
 */

// -----------------------------------------------------------------------------
/// Lightning observation operator class
class ObsLightning : public ObsOperatorBase,
                     private util::ObjectCounter<ObsLightning> {
 public:
  typedef ObsLightningParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() {return "ufo::ObsLightning";}

  ObsLightning(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsLightning();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  void computeReducedVars(const oops::Variables & reducedVars, GeoVaLs & geovals) const override;

  Locations_ locations() const override;

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  size_t nhoriz_;
  // bool l_fed_nonlinear_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_LIGHTNING_OBSLIGHTNING_H_
