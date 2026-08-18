/*
 * (C) Copyright 2025-2026 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_CRTM_OBSEXTCOEFFPROFCRTM_H_
#define UFO_OPERATORS_CRTM_OBSEXTCOEFFPROFCRTM_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/Fortran.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/crtm/crtmParameters/ObsExtCoeffProfCRTMParameters.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// ExtCoeffProfCRTM observation for UFO.
class ObsExtCoeffProfCRTM : public ObsOperatorBase,
                    private util::ObjectCounter<ObsExtCoeffProfCRTM> {
 public:
  static const std::string classname() {return "ufo::ObsExtCoeffProfCRTM";}
  typedef ObsExtCoeffProfCRTMParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  ObsExtCoeffProfCRTM(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsExtCoeffProfCRTM();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperExtCoeffProfCRTM_;}
  const int & toFortran() const {return keyOperExtCoeffProfCRTM_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperExtCoeffProfCRTM_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  ObsExtCoeffProfCRTMParameters parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_CRTM_OBSEXTCOEFFPROFCRTM_H_
