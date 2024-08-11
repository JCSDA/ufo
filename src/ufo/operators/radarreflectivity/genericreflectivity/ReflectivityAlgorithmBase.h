/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_REFLECTIVITYALGORITHMBASE_H_
#define UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_REFLECTIVITYALGORITHMBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"
#include "oops/util/parameters/HasParameters_.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

/// \brief ReflectivityAlgorithm parameters base class.
class ReflectivityAlgorithmParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(ReflectivityAlgorithmParametersBase, oops::Parameters)
 public:
  /// Name of the reflectivity algorithm
  /// Valid names are specified using a `ReflectivityAlgorithmMaker` in subclasses of
  /// ReflectivityAlgorithmBase.
  oops::RequiredParameter<std::string> reflectivityAlgorithmName{"name", this};
};

/// \brief Concrete class containing the generic parameters used by a subclass of
/// ReflectivityAlgorithmBase *if* a dedicated set of parameters have not been specified for
/// that subclass.
class GenericReflectivityAlgorithmParameters : public ReflectivityAlgorithmParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericReflectivityAlgorithmParameters,
                           ReflectivityAlgorithmParametersBase)
    // There are no generic parameters.
};

/// \brief ReflectivityAlgorithm base class.
///
/// Subclasses of this class must implement the following methods:
/// * `simulateObsImpl`: implementation of the simulateObs method for this algorithm.
/// * `setTrajectoryImpl`: implementation of the setTrajectory method for this algorithm.
/// * `simulateObsTLImpl`: implementation of the simulateObsTL method for this algorithm.
/// * `simulateObsADImpl`: implementation of the simulateObsAD method for this algorithm.
///
/// The `setTrajectoryImpl` routine can be used to store GeoVaLs that represent the
/// full nonlinear trajectory, or any derived quantities such as interpolated values,
/// if they are required in the TL/AD code.
///
/// The `printImpl` method can be overridden to provide information on the
/// algorithm in use.
///
/// Any required variables for the nonlinear and linear operators should be specified
/// in the subclass constructor.
///
/// If the ObsRadarReflectivity operator is used as one component of the Composite operator,
/// the output H(x), TL, and AD vectors are automatically filled at the correct positions.
///
/// For a concrete example of the above see MetOfficeReflectivity.h.
class ReflectivityAlgorithmBase {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;

  explicit ReflectivityAlgorithmBase(const ReflectivityAlgorithmParametersBase &,
                                     const ioda::ObsSpace &,
                                     const int,
                                     oops::Variables &,
                                     oops::Variables &);
  virtual ~ReflectivityAlgorithmBase() {}

  /// Public interface to the `simulateObsImpl` method.
  void simulateObs(const GeoVaLs &,
                   ioda::ObsVector &,
                   ObsDiagnostics &,
                   const QCFlags_t &) const;

  /// Public interface to the `setTrajectoryImpl` method.
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &);

  /// Public interface to the `simulateObsTLImpl` method.
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const;

  /// Public interface to the `simulateObsADImpl` method.
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const;

  /// Public interface to the `printImpl` method.
  void print(std::ostream & os) const;

 protected:
  /// Observation operator.
  /// Abstract function to be overwritten in subclasses.
  virtual void simulateObsImpl(const GeoVaLs &,
                               ioda::ObsVector &,
                               ObsDiagnostics &,
                               const QCFlags_t &) const = 0;

  /// Set trajectory for TL/AD operator.
  /// Abstract function to be overwritten in subclasses.
  virtual void setTrajectoryImpl(const GeoVaLs &,
                                 ObsDiagnostics &,
                                 const QCFlags_t &) = 0;

  /// TL operator.
  /// Abstract function to be overwritten in subclasses.
  virtual void simulateObsTLImpl(const GeoVaLs &,
                                 ioda::ObsVector &,
                                 const QCFlags_t &) const = 0;

  /// AD operator.
  /// Abstract function to be overwritten in subclasses.
  virtual void simulateObsADImpl(GeoVaLs &,
                                 const ioda::ObsVector &,
                                 const QCFlags_t &) const = 0;

  /// Print operator information.
  virtual void printImpl(std::ostream & os) const {
    os << "printing not implemented for this subclass of ReflectivityAlgorithmBase." << std::endl;
  }

  /// ObsSpace.
  const ioda::ObsSpace & obsdb_;

  /// Index of reflectivity in the operator variables. This ensures that H(x), TL and AD values
  /// are filled at the correct positions if the operator is one component of the Composite
  /// operator.
  const int idxReflectivity_;
};

/// Reflectivity algorithm factory.
class ReflectivityAlgorithmFactory {
 public:
  static std::unique_ptr<ReflectivityAlgorithmBase>
    create(const ReflectivityAlgorithmParametersBase &,
           const ioda::ObsSpace &,
           const int,
           oops::Variables &,
           oops::Variables &);

  static std::unique_ptr<ReflectivityAlgorithmParametersBase>
    createParameters(const std::string &name);

  /// \brief Return the names of all algorithms that can be created by one of the
  /// registered makers.
  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~ReflectivityAlgorithmFactory() = default;

 protected:
  explicit ReflectivityAlgorithmFactory(const std::string &);

 private:
  virtual std::unique_ptr<ReflectivityAlgorithmBase>
    make(const ReflectivityAlgorithmParametersBase &,
         const ioda::ObsSpace &,
         const int,
         oops::Variables &,
         oops::Variables &) = 0;

  virtual std::unique_ptr<ReflectivityAlgorithmParametersBase> makeParameters() const = 0;

  static std::map <std::string, ReflectivityAlgorithmFactory*> & getMakers() {
    static std::map <std::string, ReflectivityAlgorithmFactory*> makers_;
    return makers_;
  }
};

template<class T>
class ReflectivityAlgorithmMaker : public ReflectivityAlgorithmFactory {
  typedef oops::TParameters_IfAvailableElseFallbackType_t<T, GenericReflectivityAlgorithmParameters>
    Parameters_;

  std::unique_ptr<ReflectivityAlgorithmBase>
    make(const ReflectivityAlgorithmParametersBase & params,
         const ioda::ObsSpace & obsdb,
         const int idxRefl,
         oops::Variables & reqvars,
         oops::Variables & reqvarsTL) override {
    const auto & stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return std::unique_ptr<ReflectivityAlgorithmBase>
      (new T(stronglyTypedParams, obsdb, idxRefl, reqvars, reqvarsTL));
  }

  std::unique_ptr<ReflectivityAlgorithmParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit ReflectivityAlgorithmMaker(const std::string & name)
    : ReflectivityAlgorithmFactory(name) {}
};

}  // namespace ufo

#endif  // UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_REFLECTIVITYALGORITHMBASE_H_
