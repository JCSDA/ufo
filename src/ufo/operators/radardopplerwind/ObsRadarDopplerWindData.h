/*
 * (C) Copyright 2024 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_RADARDOPPLERWIND_OBSRADARDOPPLERWINDDATA_H_
#define UFO_OPERATORS_RADARDOPPLERWIND_OBSRADARDOPPLERWINDDATA_H_

#include <string>
#include <vector>

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/operators/radardopplerwind/ObsRadarDopplerWindParameters.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {

template <typename comparatorType>
  bool anyEqualTo(const comparatorType comparator, const comparatorType x) {
  return x == comparator;
}

template <typename comparatorType, typename... Args>
  bool anyEqualTo(const comparatorType comparator, const comparatorType x, const Args... y) {
  return x == comparator || anyEqualTo(comparator, y...);
}

/// \brief Data handler for the RadarDopplerWind observation operator and TL/AD code.
class ObsRadarDopplerWindData  {
 public:
  typedef ObsRadarDopplerWindParameters Parameters_;
  typedef std::vector<std::vector<double>> VertCoordGeoVaLs_;

  explicit ObsRadarDopplerWindData(const ioda::ObsSpace &,
                                   const Parameters_ &);

  /// Return required variables for the operator.
  const oops::Variables & requiredVars() const {return requiredVars_;}

  /// Cache vertical coordinate GeoVaLs.
  void cacheVertCoordGeoVaLs(const GeoVaLs &) const;

  /// Return GeoVaLs for vertical coordinate associated with horizontal wind.
  const VertCoordGeoVaLs_ & heightGeoVaLs_uv() const {
    return vec_z_uv_;
  }

  /// Return GeoVaLs for vertical coordinate associated with vertical wind.
  const VertCoordGeoVaLs_ & heightGeoVaLs_w() const {
    return vec_z_w_;
  }

  /// Fill H(x) and tangent linear equivalent.
  void fillHofX(const ioda::ObsSpace &, const GeoVaLs &, ioda::ObsVector &) const;

  void print(std::ostream & os) const {
    os << "ObsRadarDopplerWind operator" << std::endl;
  }

 private:
  /// ObsSpace.
  const ioda::ObsSpace& odb_;

  /// Required variables.
  oops::Variables requiredVars_;

  /// Cache of vertical coordinate GeoVaL associated with horizontal wind.
  mutable VertCoordGeoVaLs_ vec_z_uv_;

  /// Cache of vertical coordinate GeoVaL associated with vertical wind.
  mutable VertCoordGeoVaLs_ vec_z_w_;

  /// Name of vertical coordinate GeoVaL associated with horizontal wind.
  const oops::Variable verticalCoordinate_uv_;

  /// Name of vertical coordinate GeoVaL associated with vertical wind.
  const oops::Variable verticalCoordinate_w_;
};

}  // namespace ufo
#endif  // UFO_OPERATORS_RADARDOPPLERWIND_OBSRADARDOPPLERWINDDATA_H_
