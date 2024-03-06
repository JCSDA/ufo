/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GROUNDGNSS_REFROPP1D_OBSGNSSGBREFROPP1D_H_
#define UFO_OPERATORS_GROUNDGNSS_REFROPP1D_OBSGNSSGBREFROPP1D_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/groundgnss/RefROPP1D/ObsGnssGBRefROPP1D.interface.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

class GnssGBRefROPP1DParameters: public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(GnssGBRefROPP1DParameters, ObsOperatorParametersBase)
 public:
  // no options
};

// -----------------------------------------------------------------------------

/// GnssGBRefROPP1D observation operator
class ObsGnssGBRefROPP1D : public ObsOperatorBase,
                        private util::ObjectCounter<ObsGnssGBRefROPP1D> {
 public:
  typedef GnssGBRefROPP1DParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGnssGBRefROPP1D";}

  ObsGnssGBRefROPP1D(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsGnssGBRefROPP1D();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssGBRefROPP1D_;}
  const int & toFortran() const {return keyOperGnssGBRefROPP1D_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssGBRefROPP1D_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_GROUNDGNSS_REFROPP1D_OBSGNSSGBREFROPP1D_H_
