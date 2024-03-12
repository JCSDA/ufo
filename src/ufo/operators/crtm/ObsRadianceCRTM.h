/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_CRTM_OBSRADIANCECRTM_H_
#define UFO_OPERATORS_CRTM_OBSRADIANCECRTM_H_

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/Fortran.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/crtm/crtmParameters/ObsRadianceCRTMParameters.h"

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
  Locations_ locations() const override;
  void computeReducedVars(const oops::Variables &, GeoVaLs &) const override;

  int & toFortran() {return keyOperRadianceCRTM_;}
  const int & toFortran() const {return keyOperRadianceCRTM_;}

 private:
  void fillReducedVarsByMaskedAveraging(GeoVaLs &) const;
  void fillReducedVarsByMaskedCopy(GeoVaLs &) const;

  void print(std::ostream &) const override;

  F90hop keyOperRadianceCRTM_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  ObsRadianceCRTMParameters parameters_;

  bool do_fov_average_;
  // Number of sampling points per semi-axis; number of total samples is proportional to (2n+1)^2
  size_t fov_sample_resol_;
  // Buffer to store FOV sample weights; these are computed once in the method locations, then are
  // used in simulateObs. The weights could be computed in simulateObs, but this would lead to
  // redundant recomputations for applications with several outer minimization loops.
  mutable std::vector<double> sample_weights_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_CRTM_OBSRADIANCECRTM_H_
