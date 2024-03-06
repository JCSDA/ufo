/*
 * (C) Copyright 2024 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARDOPPLERWIND_OBSRADARDOPPLERWIND_H_
#define UFO_OPERATORS_RADARDOPPLERWIND_OBSRADARDOPPLERWIND_H_

#include <string>

#include "ioda/ObsDataVector.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/radardopplerwind/ObsRadarDopplerWindData.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// \brief Radar Doppler wind observation operator.
///
/// This operator computes the model equivalent of radial velocity [m/s].
///
/// The following variables must have been produced prior to running this operator:
/// * MetaData/sinTilt
/// * MetaData/cosAzimuthCosTilt
/// * MetaData/sinAzimuthCosTilt
/// * MetaData/height
///
/// These variables can be produced with the `RadarBeamGeometry` Variable Transform filter.
/// Tilt is the angle [deg] of the radar beam relative to the horizontal, and
/// and azimuth [deg] is the angle of the radar beam measured clockwise relative to true North.
///
/// If the model has a staggered grid, it is possible to define different associated coordinates
/// for the horizontal (u, v) and vertical (w) wind speeds. That will ensure the correct
/// vertical interpolation is performed.
///
/// An example yaml configuration is as follows:
///  obs operator:
///    name: RadarDopplerWind
///    vertical coordinate for horizontal wind: height_levels
///    vertical coordinate for vertical wind: height
///
class ObsRadarDopplerWind : public ObsOperatorBase,
  private util::ObjectCounter<ObsRadarDopplerWind> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsRadarDopplerWindParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() {return "ufo::ObsRadarDopplerWind";}

  ObsRadarDopplerWind(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsRadarDopplerWind() override;

  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t&) const override;

  const oops::Variables & requiredVars() const override { return data_.requiredVars(); }

 private:
  void print(std::ostream &) const override;

 private:
  /// ObsSpace.
  const ioda::ObsSpace& odb_;

  /// Data handler for the RadarDopplerWind operator and TL/AD code.
  ObsRadarDopplerWindData data_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RADARDOPPLERWIND_OBSRADARDOPPLERWIND_H_
