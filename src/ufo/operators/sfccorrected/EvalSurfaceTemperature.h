/*
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_SFCCORRECTED_EVALSURFACETEMPERATURE_H_
#define UFO_OPERATORS_SFCCORRECTED_EVALSURFACETEMPERATURE_H_

#include <string>
#include <vector>

#include "ufo/operators/sfccorrected/SurfaceOperatorBase.h"

/// Forward declarations
namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

class GeoVaLs;

class airTemperatureAt2M_WRFDA : public SurfaceOperatorBase {
 public:
  explicit airTemperatureAt2M_WRFDA(const std::string &,
                                    const Parameters_ &);
  virtual ~airTemperatureAt2M_WRFDA() {}

  void simobs(const ufo::GeoVaLs &,
              const ioda::ObsSpace &,
              std::vector<float> &) const override;
};

class airTemperatureAt2M_UKMO : public SurfaceOperatorBase {
 public:
  explicit airTemperatureAt2M_UKMO(const std::string &,
                                   const Parameters_ &);
  virtual ~airTemperatureAt2M_UKMO() {}

  void simobs(const ufo::GeoVaLs &,
              const ioda::ObsSpace &,
              std::vector<float> &) const override;
};

class airTemperatureAt2M_GSL : public SurfaceOperatorBase {
 public:
  explicit airTemperatureAt2M_GSL(const std::string &,
                                  const Parameters_ &);
  virtual ~airTemperatureAt2M_GSL() {}

  void simobs(const ufo::GeoVaLs &,
              const ioda::ObsSpace &,
              std::vector<float> &) const override;

/// TL/AD for the GSL surface-temperature correction.
///
/// Forward formula:
///   hofx = T_2m + lapse_rate * (z_sfc - z_obs)
/// where lapse_rate is taken from the YAML option:
///   - Constant      : a fixed scalar from the config.
///   - NoAdjustment  : zero (the second term drops out).
///   - Local         : (T_lev1 - T_top) / (z_top - z_lev1), with
///                     optional min/max clamping.
///
/// Heights and obs_height are treated as fixed (not active increment
/// variables), so they only enter through the saved trajectory.
///
/// dhofx = dT_2m + dlapse_rate * (z_sfc - z_obs).
/// lapse_rate has zero derivative w.r.t. the model state for Constant
/// and NoAdjustment (both are fixed scalars), so dhofx reduces to dT_2m
/// Only "Local" lapse rate depends on model state, via T_lev1 and T_top.
///
/// The "tlad lapse rate option" parameter controls how the lapse rate
/// is differentiated:
///   - MatchSimobs (default): differentiate the configured calculation;
///     in Local mode this adds dT_lev1, dT_top contributions, gated off
///     wherever the trajectory-time clamp was active.
///   - NoAdjustment: freeze the per-location lapse rate at trajectory
///     time and drop any dT_lev1, dT_top contributions (dhofx = dT_2m,
///     same as Constant/NoAdjustment).
///
/// The model's TL/AD may not carry an independent active increment for
/// air_temperature_at_2m (when only exists as a diagnostic),
/// in this case dT_2m is zero and the cases above (Constant,
/// NoAdjustment, or frozen Local) would give dhofx = 0 and no
/// sensitivity to the control vector at all.
/// Setting "tlad spread to lowest level" to true adds dT_lev1 in
/// exactly those cases, so the increment still reaches the 3D prognostic
/// field. Defaults to false (exact linearization, dhofx = dT_2m alone).
  void settraj(const GeoVaLs &,
               const ioda::ObsSpace &) override;
  void TL(const GeoVaLs &,
          std::vector<double> &) const override;
  void AD(GeoVaLs &,
          std::vector<double> &) const override;

 private:
  /// True if the linearization needs dT_lev1/dT_top contributions, i.e.
  /// simobs mode is Local AND tlad mode is MatchSimobs.
  bool needsLevelTemperature_;

  /// TL/AD mode, set in settraj:
  /// Local + MatchSimobs differentiates the lapse rate
  /// optionally spread to the lowest level
  bool tlad_local_path_;
  bool tlad_spread_to_lowest_level_;

  /// Fetch dT_lev1 (and dT_top on the Local path) from the GeoVaLs for TL/AD.
  void getLevelIncrements(const GeoVaLs &,
                          std::vector<double> &,
                          std::vector<double> &) const;

  /// Trajectory state (sized at setTrajectory).
  /// Stored in double for precision in the linear operator.
  std::vector<double> dz_sfc_minus_obs_;   // (z_sfc - z_obs) in m
  std::vector<bool> valid_;                // per-location validity mask

  /// Local-mode trajectory state (only populated when needsLevelTemperature_).
  std::vector<double> inv_lapse_layer_depth_;  // 1 / (z_top - z_lev1), Local lapse-rate layer
  std::vector<bool> lapse_rate_clamped_;      // true where min/max threshold changed the lapse rate
  int surface_level_index_;
  int lapse_top_level_index_;                // top level of the Local lapse-rate layer
};

}  // namespace ufo

#endif  // UFO_OPERATORS_SFCCORRECTED_EVALSURFACETEMPERATURE_H_
