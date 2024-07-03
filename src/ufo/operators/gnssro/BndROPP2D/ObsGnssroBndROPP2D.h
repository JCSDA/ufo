/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GNSSRO_BNDROPP2D_OBSGNSSROBNDROPP2D_H_
#define UFO_OPERATORS_GNSSRO_BNDROPP2D_OBSGNSSROBNDROPP2D_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/gnssro/BndROPP2D/ObsGnssroBndROPP2D.interface.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

class GnssroBndROPP2DOptionsParameters: public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GnssroBndROPP2DOptionsParameters, Parameters)
 public:
  oops::Parameter<int> useCompress{"use_compress", 0, this};
  oops::Parameter<size_t> nHoriz{"n_horiz", 31, this};
  oops::Parameter<double> res{"res", 40.0, this};
  oops::Parameter<double> top2D{"top_2d", 20.0, this};
  oops::Parameter<std::string> roType{"ro_type", "spaceborne", this};
};

class GnssroBndROPP2DParameters: public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(GnssroBndROPP2DParameters, ObsOperatorParametersBase)
 public:
  oops::Parameter<GnssroBndROPP2DOptionsParameters> options{"obs options", {}, this};
};

// -----------------------------------------------------------------------------

/// GnssroBndROPP2D observation operator
class ObsGnssroBndROPP2D : public ObsOperatorBase,
                        private util::ObjectCounter<ObsGnssroBndROPP2D> {
 public:
  typedef GnssroBndROPP2DParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGnssroBndROPP2D";}

  ObsGnssroBndROPP2D(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsGnssroBndROPP2D();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  void computeReducedVars(const oops::Variables & reducedVars, GeoVaLs & geovals) const override;

  Locations_ locations() const override;

  int & toFortran() {return keyOperGnssroBndROPP2D_;}
  const int & toFortran() const {return keyOperGnssroBndROPP2D_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroBndROPP2D_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
  const size_t nhoriz_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_BNDROPP2D_OBSGNSSROBNDROPP2D_H_
