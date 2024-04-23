/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_SUPEROB_SUPEROBRADAR_H_
#define UFO_SUPEROB_SUPEROBRADAR_H_

#include <string>
#include <vector>

#include "ufo/filters/DiagnosticFlag.h"

#include "ufo/superob/SuperObBase.h"

namespace ufo {

/// Parameters associated with the SuperObRadar class.
class SuperObRadarParameters : public SuperObParametersBase {
  OOPS_CONCRETE_PARAMETERS(SuperObRadarParameters, SuperObParametersBase)
 public:
  oops::Parameter<double> numBeamsInSuperObRegion{"number of beams in superob region", 5, this};
  oops::Parameter<double> superObRegionRadialExtent
    {"superob region radial extent [m]", 5000.0, this};
  oops::Parameter<int> superObMinObs{"minimum number of observations in superob region", 5, this};
};

/// \brief Radar superobbing algorithm.
///
/// \details This algorithm computes superobs for radar scans.
/// Each scan is divided into superob regions according to the parameters
/// `number of beams in superob region` and `superob region radial extent [m]`.
/// Observation (O) and background (B) values inside each region are summed,
/// discarding those observations that are masked by a superob template class.
/// (See the `SuperObRadarTemplate` class documentation for futher details.)
/// Note there are typically multiple regions in a scan, so there can be multiple
/// superobs computed.
///
/// If there are insufficient observations inside the region (governed by the parameter
/// `minimum number of observations in superob region`, a superob is not computed.
///
/// In order to compute the superob, the mean innovation (O - B) is calculated and added onto
/// the background value that lies closest to the centre of the superob template.
///
/// Two superob uncertainties are also computed for use in subsequent error assignment:
/// * Total innovation uncertainty,
/// * Background uncertainty.
/// These values are written to the `TotalUncertainty` and `BackgroundUncertainty` groups
/// in the ObsSpace.
///
/// Lastly, a diagnostic flag called `UsedInSuperOb` is used to record the locations
/// whose values of O and B were used to compute the superob in each case.
class SuperObRadar : public SuperObBase {
 public:
  typedef SuperObRadarParameters Parameters_;

  explicit SuperObRadar(const Parameters_ &,
                        const ObsFilterData &,
                        const std::vector<bool> &,
                        const Variables &,
                        const ioda::ObsDataVector<int> &,
                        std::vector<std::vector<bool>> &);
  ~SuperObRadar() {}

 private:
  /// Compute superob values and errors for each record.
  void computeSuperOb(const std::vector<std::size_t> &,
                      const std::vector<float> &,
                      const std::vector<float> &,
                      const ioda::ObsDataRow<int> &,
                      std::vector<float> & superobs,
                      std::vector<bool> &) const override;

  /// Save auxiliary variables the ObsSpace when the algorithm has finished.
  void saveAuxiliaryVariables(const std::string & variableName) const override;

  const Parameters_ & params_;

  std::vector<int> numBeams_;
  std::vector<int> numGates_;
  std::vector<float> beamWidth_;
  std::vector<float> gateWidth_;
  std::vector<float> latitude_;
  std::vector<float> longitude_;
  std::vector<float> minGateRange_;
  std::vector<float> beamTiltAngle_;
  std::vector<float> stationElevation_;
  std::vector<float> gateRange_;
  std::vector<float> beamAzimuthAngle_;

  /// Diagnostic flag indicating whether values of O and B at a location
  /// were used to compute a superob.
  mutable std::vector<DiagnosticFlag> usedInSuperOb_;
  /// Total uncertainty for a superob.
  mutable std::vector<float> totalUncertainty_;
  /// Background uncertainty for a superob.
  mutable std::vector<float> backgroundUncertainty_;
};

}  // namespace ufo

#endif  // UFO_SUPEROB_SUPEROBRADAR_H_
