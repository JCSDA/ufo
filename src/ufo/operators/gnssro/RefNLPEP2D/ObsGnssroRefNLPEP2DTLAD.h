/*
 * (C) Copyright 2025 Space Sciences and Engineering, LLC (dba PlanetiQ).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * author Steve Marshall (smarshall@planetiq.com)
 */
#ifndef UFO_OPERATORS_GNSSRO_REFNLPEP2D_OBSGNSSROREFNLPEP2DTLAD_H_
#define UFO_OPERATORS_GNSSRO_REFNLPEP2D_OBSGNSSROREFNLPEP2DTLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/gnssro/RefNLPEP2D/ObsGnssroRefNLPEP2DParameters.h"
#include "ufo/operators/gnssro/utils/GnssroProfileRayPaths.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathOrchestrator.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathParameters.h"
#include "ufo/operators/gnssro/utils/RefractivityCalculator.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

// -----------------------------------------------------------------------------
/// GnssroRefNLPEP2D TL/AD observation operator class
class ObsGnssroRefNLPEP2DTLAD : public LinearObsOperatorBase {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsGnssroRefNLPEP2DParameters Parameters_;
  typedef ObsGnssroRefNLPEP2DOptions ObsOptions_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGnssroRefNLPEP2DTLAD";}


  ObsGnssroRefNLPEP2DTLAD(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsGnssroRefNLPEP2DTLAD() override;

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return *varin_;}

 private:
  void print(std::ostream &) const override;

  // Member data
  std::unique_ptr<oops::Variables> varin_;
  ObsOptions_ opts_;
  GnssroRayPathParameters rpParams_;
  GnssroRayPathOrchestrator orchestrator_;
  std::unique_ptr<RefractivityCalculator> refrCalc_;
  bool trajectorySet_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_GNSSRO_REFNLPEP2D_OBSGNSSROREFNLPEP2DTLAD_H_
