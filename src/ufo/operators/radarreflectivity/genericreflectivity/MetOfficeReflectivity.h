/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_METOFFICEREFLECTIVITY_H_
#define UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_METOFFICEREFLECTIVITY_H_

#include <memory>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/operators/radarreflectivity/genericreflectivity/ReflectivityAlgorithmBase.h"

namespace ufo {

/// \brief Parameters associated with the MetOfficeReflectivity class.
class MetOfficeReflectivityParameters : public ReflectivityAlgorithmParametersBase {
  OOPS_CONCRETE_PARAMETERS(MetOfficeReflectivityParameters,
                           ReflectivityAlgorithmParametersBase)

 public:
    oops::RequiredParameter<float>
      rainMultiplier{"rain multiplier",
        "Multiplier for rain component of squared reflectivity.",
        this};
    oops::RequiredParameter<float>
      rainExponent{"rain exponent",
        "Exponent for rain component of squared reflectivity.",
        this};
    oops::RequiredParameter<float>
      iceMultiplier{"ice multiplier",
        "Multiplier for ice component of squared reflectivity.",
        this};
    oops::RequiredParameter<float>
      iceAdditive{"ice additive constant",
        "Additive constant for ice component of squared reflectivity.",
        this};
    oops::RequiredParameter<float>
      iceExponent{"ice exponent",
        "Exponent for ice component of squared reflectivity.",
        this};
    oops::RequiredParameter<float>
      lowerBound{"lower bound",
        "Lower bound of squared reflectivity.",
        this};
};

/// \brief Met Office reflectivity algorithm.
///
/// \details This algorithm computes the H(x) of reflectivity accounting for contributions
/// from both rain and ice.
/// Internally, the model mixing ratios for rain and ice are referred to as
/// `qrain` and `qice` respectively.
/// There are several semi-empirical parameters that control the relationship between
/// the model mixing ratios and the reflectivity.
/// These parameters are tuned for specific combinations of radar sites and NWP models.
///
/// It is possible to use this algorithm in TL/AD operator but that is discouraged
/// due to the high nonlinearity present.
///
/// An example usage of the algorithm is as follows:
///
/// - obs operator:
///     name: RadarReflectivity
///     algorithm:
///       name: Met Office reflectivity
///       rain multiplier: 1000.0
///       rain exponent: 1.1
///       ice multiplier: 0.1
///       ice additive constant: 2.0
///       ice exponent: 2.0
///       lower bound: 1.0
///
class MetOfficeReflectivity : public ReflectivityAlgorithmBase {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  typedef MetOfficeReflectivityParameters Parameters_;

  explicit MetOfficeReflectivity(const Parameters_ &,
                                 const ioda::ObsSpace &,
                                 const int,
                                 oops::Variables &,
                                 oops::Variables &);
  ~MetOfficeReflectivity() {}

 private:  // functions
  /// Observation operator.
  void simulateObsImpl(const GeoVaLs &, ioda::ObsVector &,
                       ObsDiagnostics &, const QCFlags_t &) const override;

  /// Set trajectory for TL/AD operator.
  void setTrajectoryImpl(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;

  /// TL operator.
  void simulateObsTLImpl(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;

  /// AD operator.
  void simulateObsADImpl(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  /// Print operator information.
  void printImpl(std::ostream & os) const override;

 private:  // variables
  /// Parameters. Note: this must not be a const reference because the variable
  /// with which it is initialised goes out of scope.
  const Parameters_ params_;

  /// Stored height GeoVaLs.
  mutable std::vector<std::vector<double>> vec_gv_z_;

  /// Stored interpolated T.
  mutable std::vector<double> vec_interp_T_;

  /// Stored interpolated p.
  mutable std::vector<double> vec_interp_p_;

  /// Stored interpolated rain mixing ratio.
  mutable std::vector<double> vec_interp_qrain_;

  /// Stored interpolated ice mixing ratio.
  mutable std::vector<double> vec_interp_qice_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_METOFFICEREFLECTIVITY_H_
