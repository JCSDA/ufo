/*
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_SFCCORRECTED_SURFACEOPERATORBASE_H_
#define UFO_OPERATORS_SFCCORRECTED_SURFACEOPERATORBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "ufo/operators/sfccorrected/ObsSfcCorrectedParameters.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class GeoVaLs;

/// \brief SurfaceOperator base class.
///
/// Subclasses of this class must implement the `simobs` method.
/// Subclasses with a linearized operator override `settraj`, `TL` and `AD`
/// defaults for the linearized subclass throw "not yet implemented".
class SurfaceOperatorBase {
 public:
  typedef ObsSfcCorrectedParameters Parameters_;

  explicit SurfaceOperatorBase(const std::string &, const Parameters_ &);
  virtual ~SurfaceOperatorBase() {}

  /// Create the model equivalent of the observation
  virtual void simobs(const ufo::GeoVaLs &,
                      const ioda::ObsSpace &,
                      std::vector<float> &) const = 0;

  /// Model variables required for this operator
  const oops::Variables & requiredVars() const {return requiredVars_;}

  /// Calculate the trajectory
  virtual void settraj(const GeoVaLs &, const ioda::ObsSpace &);

  /// Calculate the tangent linear
  virtual void TL(const GeoVaLs &,
                  std::vector<double> &) const;

  /// Calculate the adjoint
  virtual void AD(GeoVaLs &,
                  std::vector<double> &) const;

 protected:
  // Required variables from the GeoVaLs
  oops::Variables requiredVars_;

  // Save a copy of the parameters
  Parameters_ params_;

  // Operator name, e.g. airTemperatureAt2M_GSL
  std::string name_;
};

/// SurfaceOperator factory.
class SurfaceOperatorFactory {
 public:
  typedef ObsSfcCorrectedParameters Parameters_;

  /// The surface operator type is determined by the combination of the observation variable
  /// and `correctionType`.
  static std::unique_ptr<SurfaceOperatorBase> create(const std::string &,
                                                     const Parameters_ &);

  /// Create one surface operator per variable for the configured `correctionType`
  /// and add each operator's required GeoVaLs to `requiredVars`.
  static std::vector<std::unique_ptr<SurfaceOperatorBase>>
  createForVariables(const oops::ObsVariables &, const Parameters_ &, oops::Variables &);

  /// \brief Return the names of all predictors that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~SurfaceOperatorFactory() = default;

 protected:
  explicit SurfaceOperatorFactory(const std::string &);

 private:
  virtual std::unique_ptr<SurfaceOperatorBase> make(const std::string &,
                                                    const Parameters_ &) = 0;

  static std::map <std::string, SurfaceOperatorFactory*> & getMakers() {
    static std::map <std::string, SurfaceOperatorFactory*> makers_;
    return makers_;
  }
};

/// \brief A subclass of SurfaceOperatorFactory able to create instances of T
/// (a concrete subclass of SurfaceOperatorBase)
template<class T>
class SurfaceOperatorMaker : public SurfaceOperatorFactory {
  typedef ObsSfcCorrectedParameters Parameters_;

  std::unique_ptr<SurfaceOperatorBase> make(const std::string & method,
                                            const Parameters_ & params) override {
    return std::make_unique<T>(method, params);
  }

 public:
  explicit SurfaceOperatorMaker(const std::string & method) :
    SurfaceOperatorFactory(method)
  {}
};

}  // namespace ufo

#endif  // UFO_OPERATORS_SFCCORRECTED_SURFACEOPERATORBASE_H_

