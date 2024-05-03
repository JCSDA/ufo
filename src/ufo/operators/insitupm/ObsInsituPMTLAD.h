/*
 * (C) Copyright 2021.
 *
 * This software is developed by NOAA/NWS/EMC under the Apache 2.0 license
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_INSITUPM_OBSINSITUPMTLAD_H_
#define UFO_OPERATORS_INSITUPM_OBSINSITUPMTLAD_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/insitupm/ObsInsituPMParameters.h"
#include "ufo/operators/insitupm/ObsInsituPMTLAD.interface.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

// -----------------------------------------------------------------------------
/// InsituPM TL/AD observation operator class
class ObsInsituPMTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsInsituPMTLAD> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;

  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsInsituPMParameters Parameters_;

  static const std::string classname() {return "ufo::ObsInsituPMTLAD";}

  ObsInsituPMTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsInsituPMTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_INSITUPM_OBSINSITUPMTLAD_H_
