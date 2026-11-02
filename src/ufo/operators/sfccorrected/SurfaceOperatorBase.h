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
#include "oops/util/missingValues.h"
#include "ufo/operators/sfccorrected/ObsSfcCorrectedParameters.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class GeoVaLs;

/// Converts a single profile of geopotential heights at one location to geometric height.
/// Missing levels are left unchanged; nothing is done if the latitude is missing.
template <typename T>
void convertToGeometricHeight(const float lat, std::vector<T> & heights) {
  if (lat == util::missingValue<float>()) return;
  const T missing = util::missingValue<T>();
  for (T & h : heights) {
    if (h != missing) {
      h = static_cast<T>(formulas::Geopotential_to_Geometric_Height(lat, static_cast<float>(h)));
    }
  }
}

/// Converts one geopotential height per location to geometric height, using the latitude
/// at each location. Locations with a missing latitude or height are left unchanged.
template <typename T>
void convertToGeometricHeight(const std::vector<float> & lats, std::vector<T> & heights) {
  const float missingFloat = util::missingValue<float>();
  const T missing = util::missingValue<T>();
  for (size_t iloc = 0; iloc < heights.size(); ++iloc) {
    if (lats[iloc] != missingFloat && heights[iloc] != missing) {
      heights[iloc] = static_cast<T>(formulas::Geopotential_to_Geometric_Height(
          lats[iloc], static_cast<float>(heights[iloc])));
    }
  }
}

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

