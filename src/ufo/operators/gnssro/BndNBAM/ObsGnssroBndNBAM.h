/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_GNSSRO_BNDNBAM_OBSGNSSROBNDNBAM_H_
#define UFO_OPERATORS_GNSSRO_BNDNBAM_OBSGNSSROBNDNBAM_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/gnssro/BndNBAM/ObsGnssroBndNBAM.interface.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

class GnssroBndNBAMOptionsParameters: public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GnssroBndNBAMOptionsParameters, Parameters)
 public:
  oops::Parameter<int> useCompress{"use_compress", 1, this};
  oops::Parameter<std::string> vertLayer{"vertlayer", "full", this};
  oops::Parameter<int> srSteps{"sr_steps", 2, this};
// Set "super_ref_qc" as "NBAM_GEOS" for the NBAM used in GEOS GSI (July 2023).
// It should be "NBAM" in GEOS's configuration once it is adapted in GEOS.
  oops::Parameter<std::string> superRefQC{"super_ref_qc", "NBAM", this};
  oops::Parameter<std::string> outputDiags{"output_diags", "false", this};
  oops::Parameter<std::string> gsiVersion{"GSI_version", "EMC", this};
  oops::Parameter<int> nlevAdd{"nlevadd", 13, this};
  oops::Parameter<int> modelTop{"modeltop", 0, this};
  oops::Parameter<int> ngrd{"ngrd", 80, this};
};

class GnssroBndNBAMParameters: public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(GnssroBndNBAMParameters, ObsOperatorParametersBase)
 public:
  oops::Parameter<GnssroBndNBAMOptionsParameters> options{"obs options", {}, this};
};

// -----------------------------------------------------------------------------
// Gnssro BndNBAM observation operator
//   -- to reproduce exactly the operational (2019) NBAM
class ObsGnssroBndNBAM : public ObsOperatorBase,
                        private util::ObjectCounter<ObsGnssroBndNBAM> {
 public:
  typedef GnssroBndNBAMParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGnssroBndNBAM";}

  ObsGnssroBndNBAM(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsGnssroBndNBAM();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroBndNBAM_;}
  const int & toFortran() const {return keyOperGnssroBndNBAM_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroBndNBAM_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_BNDNBAM_OBSGNSSROBNDNBAM_H_
