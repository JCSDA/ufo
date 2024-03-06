/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_CRTM_OBSRADIANCECRTM_H_
#define UFO_OPERATORS_CRTM_OBSRADIANCECRTM_H_

#include <iostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/Fortran.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/crtm/crtmParameters/ObsRadianceCRTMParameters.h"
#include "ufo/operators/crtm/ObsRadianceCRTM.interface.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// RadianceCRTM observation for UFO.
class ObsRadianceCRTM : public ObsOperatorBase,
                    private util::ObjectCounter<ObsRadianceCRTM> {
 public:
  static const std::string classname() {return "ufo::ObsRadianceCRTM";}
  typedef ObsRadianceCRTMParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  ObsRadianceCRTM(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsRadianceCRTM();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperRadianceCRTM_;}
  const int & toFortran() const {return keyOperRadianceCRTM_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperRadianceCRTM_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  ObsRadianceCRTMParameters parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_CRTM_OBSRADIANCECRTM_H_
