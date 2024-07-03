/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GNSSRO_BNDROPP1D_OBSGNSSROBNDROPP1D_H_
#define UFO_OPERATORS_GNSSRO_BNDROPP1D_OBSGNSSROBNDROPP1D_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/gnssro/BndROPP1D/ObsGnssroBndROPP1D.interface.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

class GnssroBndROPP1DOptionsParameters: public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GnssroBndROPP1DOptionsParameters, Parameters)
 public:
  oops::Parameter<int> useCompress{"use_compress", 0, this};
};

class GnssroBndROPP1DParameters: public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(GnssroBndROPP1DParameters, ObsOperatorParametersBase)
 public:
  oops::Parameter<GnssroBndROPP1DOptionsParameters> options{"obs options", {}, this};
};

// -----------------------------------------------------------------------------

/// GnssroBndROPP1D observation operator
class ObsGnssroBndROPP1D : public ObsOperatorBase,
                        private util::ObjectCounter<ObsGnssroBndROPP1D> {
 public:
  typedef GnssroBndROPP1DParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGnssroBndROPP1D";}

  ObsGnssroBndROPP1D(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsGnssroBndROPP1D();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroBndROPP1D_;}
  const int & toFortran() const {return keyOperGnssroBndROPP1D_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroBndROPP1D_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_BNDROPP1D_OBSGNSSROBNDROPP1D_H_
