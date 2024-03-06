/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_CRTM_OBSAODCRTM_H_
#define UFO_OPERATORS_CRTM_OBSAODCRTM_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/Fortran.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/crtm/crtmParameters/ObsAodCRTMParameters.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// AodCRTM observation for UFO.
class ObsAodCRTM : public ObsOperatorBase,
                    private util::ObjectCounter<ObsAodCRTM> {
 public:
  static const std::string classname() {return "ufo::ObsAodCRTM";}
  typedef ObsAodCRTMParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  ObsAodCRTM(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAodCRTM();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperAodCRTM_;}
  const int & toFortran() const {return keyOperAodCRTM_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperAodCRTM_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  ObsAodCRTMParameters parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_CRTM_OBSAODCRTM_H_
