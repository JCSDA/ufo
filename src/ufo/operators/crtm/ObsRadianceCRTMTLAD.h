/*
 * (C) Copyright 2017-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_CRTM_OBSRADIANCECRTMTLAD_H_
#define UFO_OPERATORS_CRTM_OBSRADIANCECRTMTLAD_H_

#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/Fortran.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/crtm/crtmParameters/ObsRadianceCRTMParameters.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// RadianceCRTM (currently only temperature) observation for UFO.
class ObsRadianceCRTMTLAD : public LinearObsOperatorBase,
                        private util::ObjectCounter<ObsRadianceCRTMTLAD> {
 public:
  typedef ObsRadianceCRTMParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsRadianceCRTMTLAD";}

  ObsRadianceCRTMTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsRadianceCRTMTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperRadianceCRTM_;}
  const int & toFortran() const {return keyOperRadianceCRTM_;}

 private:
  /// \brief Expand the nullifyJacobianParameters' channel-range strings to integer list and
  /// save in the "linear obs operator" sub-config of \p conf.
  ///
  /// \param linearParams validated "linear obs operator" parameters, read only
  /// \param conf the config passed to the Fortran setup routine with expanded channel lists
  static void nullifyJacParamSetup(const CRTMLinearObsOperatorParameters & linearParams,
                                eckit::LocalConfiguration & conf);

  void print(std::ostream &) const override;
  F90hop keyOperRadianceCRTM_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_CRTM_OBSRADIANCECRTMTLAD_H_
