/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GNSSRO_REFNCEP_OBSGNSSROREFNCEP_H_
#define UFO_OPERATORS_GNSSRO_REFNCEP_OBSGNSSROREFNCEP_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/gnssro/RefNCEP/ObsGnssroRefNCEP.interface.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

class GnssroRefNCEPOptionsParameters: public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GnssroRefNCEPOptionsParameters, Parameters)
 public:
  oops::Parameter<int> useCompress{"use_compress", 1, this};
};

class GnssroRefNCEPParameters: public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(GnssroRefNCEPParameters, ObsOperatorParametersBase)
 public:
  oops::Parameter<GnssroRefNCEPOptionsParameters> options{"obs options", {}, this};
};

// -----------------------------------------------------------------------------

/// GnssroRefNCEP observation operator
class ObsGnssroRefNCEP : public ObsOperatorBase,
                         private util::ObjectCounter<ObsGnssroRefNCEP> {
 public:
  typedef GnssroRefNCEPParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGnssroRefNCEP";}

  ObsGnssroRefNCEP(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsGnssroRefNCEP();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroRefNCEP_;}
  const int & toFortran() const {return keyOperGnssroRefNCEP_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroRefNCEP_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_REFNCEP_OBSGNSSROREFNCEP_H_
