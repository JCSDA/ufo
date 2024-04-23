/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <Eigen/Dense>

#include <limits>

#include "ufo/superob/SuperObRadar.h"
#include "ufo/superob/SuperObRadarTemplate.h"

#include "ufo/utils/IntegralImage.h"

namespace ufo {

static SuperObMaker<SuperObRadar> makerSuperObRadar_("radar");

SuperObRadar::SuperObRadar(const Parameters_ & params,
                           const ObsFilterData & data,
                           const std::vector<bool> & apply,
                           const Variables & filtervars,
                           const ioda::ObsDataVector<int> & flags,
                           std::vector<std::vector<bool>> & flagged)
  : SuperObBase(params, data, apply, filtervars, flags, flagged),
    params_(params)
{
  oops::Log::trace() << "SuperObRadar constructor" << std::endl;

  const std::size_t nlocs = obsdb_.nlocs();
  const float missing = util::missingValue<float>();

  // Retrieve existing vectors from the ObsSpace.
  numBeams_.resize(nlocs);
  numGates_.resize(nlocs);
  beamWidth_.resize(nlocs);
  gateWidth_.resize(nlocs);
  latitude_.resize(nlocs);
  longitude_.resize(nlocs);
  minGateRange_.resize(nlocs);
  beamTiltAngle_.resize(nlocs);
  stationElevation_.resize(nlocs);
  gateRange_.resize(nlocs);
  beamAzimuthAngle_.resize(nlocs);

  obsdb_.get_db("MetaData", "numberOfBeams", numBeams_);
  obsdb_.get_db("MetaData", "numberOfGates", numGates_);
  obsdb_.get_db("MetaData", "beamWidth", beamWidth_);
  obsdb_.get_db("MetaData", "gateWidth", gateWidth_);
  obsdb_.get_db("MetaData", "stationLatitude", latitude_);
  obsdb_.get_db("MetaData", "stationLongitude", longitude_);
  obsdb_.get_db("MetaData", "minGateRange", minGateRange_);
  obsdb_.get_db("MetaData", "beamTiltAngle", beamTiltAngle_);
  obsdb_.get_db("MetaData", "stationElevation", stationElevation_);
  obsdb_.get_db("MetaData", "gateRange", gateRange_);
  obsdb_.get_db("MetaData", "beamAzimuthAngle", beamAzimuthAngle_);

  // Initialise vectors to be written out to the ObsSpace in the saveAuxiliaryVariables function.
  usedInSuperOb_.resize(nlocs, false);
  totalUncertainty_.resize(nlocs, missing);
  backgroundUncertainty_.resize(nlocs, missing);

  oops::Log::trace() << "SuperObRadar constructor finished" << std::endl;
}

void SuperObRadar::computeSuperOb(const std::vector<std::size_t> & locs,
                                  const std::vector<float> & obs,
                                  const std::vector<float> & hofx,
                                  const ioda::ObsDataRow<int> & flags,
                                  std::vector<float> & superobs,
                                  std::vector<bool> & flagged) const {
  oops::Log::trace() << "SuperObRadar::computeSuperOb" << std::endl;

  const float missing = util::missingValue<float>();

  // Obtain metadata associated with this scan (ObsSpace record).
  // All metadata values are the same, so can just take the first location in the scan
  // regardless of QC decisions.
  const std::size_t firstLocInScan = locs.front();
  SuperObRadarTemplateData scanData;
  scanData.numBeams = numBeams_[firstLocInScan];
  scanData.numGates = numGates_[firstLocInScan];
  scanData.beamWidth = beamWidth_[firstLocInScan];
  scanData.gateWidth = gateWidth_[firstLocInScan];
  scanData.stationLatitude = latitude_[firstLocInScan];
  scanData.stationLongitude = longitude_[firstLocInScan];
  scanData.stationElevation = stationElevation_[firstLocInScan];
  scanData.minGateRange = minGateRange_[firstLocInScan];
  scanData.beamTiltAngle = beamTiltAngle_[firstLocInScan];

  // Fill arrays to be processed with the integral image technique.
  // The arrays are extended at the edges.
  // todo(ctgh): Check the values used to extend the arrays here
  // and in the SuperObRadarTemplate class.
  const int imax = scanData.numGates + 2;
  const int jmax = scanData.numBeams + 3;
  Eigen::ArrayXXd obsArray = Eigen::ArrayXXd::Zero(imax, jmax);
  Eigen::ArrayXXd hofxArray = Eigen::ArrayXXd::Zero(imax, jmax);
  Eigen::ArrayXXd innovSquaredArray = Eigen::ArrayXXd::Zero(imax, jmax);
  Eigen::ArrayXXd hofxSquaredArray = Eigen::ArrayXXd::Zero(imax, jmax);
  Eigen::ArrayXXi superObCount = Eigen::ArrayXXi::Zero(imax, jmax);
  Eigen::ArrayXXi superObLocation = Eigen::ArrayXXi::Zero(imax, jmax);

  // Loop over locations in this record in order to populate the 2D superobbing arrays.
  for (int jloc : locs) {
    const float obsValue = obs[jloc];
    const float hofxValue = hofx[jloc];
    // Only consider locations which have valid observation and background values
    // and are either passing QC or have been marked as passive (i.e. H(x) is computed
    // but the observation is not assimilated).
    if ((flags[jloc] != QCflags::pass && flags[jloc] != QCflags::passive) ||
        obsValue == missing ||
        hofxValue == missing) {
      continue;
    }

    // Determine coordinates of this observation in the 2D superobbing arrays.
    const int rangeIndex =
      std::round((gateRange_[jloc] - scanData.minGateRange +
                  0.5 * scanData.gateWidth) / scanData.gateWidth);
    const int azimIndex =
      std::round((beamAzimuthAngle_[jloc] + 0.5 * scanData.beamWidth) / scanData.beamWidth);

    // Update superobbing arrays at the relevant coordinates.
    obsArray(rangeIndex, azimIndex) = obsValue;
    hofxArray(rangeIndex, azimIndex) = hofxValue;
    innovSquaredArray(rangeIndex, azimIndex) = std::pow(obsValue - hofxValue, 2.0);
    hofxSquaredArray(rangeIndex, azimIndex) = std::pow(hofxValue, 2.0);
    superObCount(rangeIndex, azimIndex) = 1;
    superObLocation(rangeIndex, azimIndex) = jloc;
  }

  // Create superob template.
  const SuperObRadarTemplate SOTemplate(params_.numBeamsInSuperObRegion,
                                        params_.superObRegionRadialExtent,
                                        scanData);

  // Mask existing arrays using the superob template.
  const Eigen::ArrayXXi superObMask = SOTemplate.getMask();
  superObCount *= superObMask;
  superObLocation *= superObMask;
  obsArray *= superObMask.cast<double>();
  hofxArray *= superObMask.cast<double>();
  innovSquaredArray *= superObMask.cast<double>();
  hofxSquaredArray *= superObMask.cast<double>();

  // Compute integral images.
  Eigen::ArrayXXi integratedSuperObCount = superObCount;
  ufo::createIntegralImage(obsArray);
  ufo::createIntegralImage(hofxArray);
  ufo::createIntegralImage(innovSquaredArray);
  ufo::createIntegralImage(hofxSquaredArray);
  ufo::createIntegralImage(integratedSuperObCount);

  // Create superobs.
  const int numBeamsInSuperObRegion = params_.numBeamsInSuperObRegion.value();
  const double superObRegionRadialExtent = params_.superObRegionRadialExtent.value();
  const int superObMinObs = params_.superObMinObs.value();
  const int istep = static_cast<int>(0.5 * (superObRegionRadialExtent / scanData.gateWidth));
  const int jstep = static_cast<int>(0.5 * numBeamsInSuperObRegion);
  const Eigen::ArrayXXd superObDistance = SOTemplate.getDistance();

  // Loop over each superob region.
  for (size_t i = istep + 1; i < scanData.numGates - istep + 1; i += (2 * istep + 1)) {
    for (size_t j = jstep + 1; j < scanData.numBeams - jstep + 1; j += (2 * jstep + 1)) {
      // Count number of observations in this region.
      const int numObsInSuperObRegion = ufo::sumIntegralImagePatch(integratedSuperObCount,
                                                                   i, j, istep, jstep);
      // Do not proceed if there are too few observations.
      if (numObsInSuperObRegion < superObMinObs) {
        continue;
      }

      // Determine superob location by searching for the valid observation with the minimum distance
      // from the centre of the superob template.
      double minObsDistance = std::numeric_limits<double>::max();
      int locSuperob = -1;
      for (size_t ik = i - istep; ik < i + istep + 1; ++ik) {
        for (size_t jk = j - jstep; jk < j + jstep + 1; ++jk) {
          if (superObCount(ik, jk) == 1) {
            usedInSuperOb_[superObLocation(ik, jk)] = true;
            if (superObDistance(ik, jk) < minObsDistance) {
              minObsDistance = superObDistance(ik, jk);
              locSuperob = superObLocation(ik, jk);
            }
          }
        }
      }

      // Do not proceed if none of the valid cells are within the minimum distance.
      if (locSuperob == -1) {
        continue;
      }

      // Compute superob value over the region.
      const double meanObs =
        ufo::sumIntegralImagePatch(obsArray, i, j, istep, jstep) / numObsInSuperObRegion;
      const double meanHofX =
        ufo::sumIntegralImagePatch(hofxArray, i, j, istep, jstep) / numObsInSuperObRegion;
      const double meanInnov = meanObs - meanHofX;
      superobs[locSuperob] = hofx[locSuperob] + meanInnov;

      // This location is un-flagged.
      flagged[locSuperob] = false;

      // Compute total variance and uncertainty of the superob.
      const double meanInnovSquared =
        ufo::sumIntegralImagePatch(innovSquaredArray, i, j, istep, jstep) / numObsInSuperObRegion;
      const double totalVariance = meanInnovSquared - std::pow(meanInnov, 2.0);
      const double totalUncertainty = std::sqrt(std::abs(totalVariance));
      totalUncertainty_[locSuperob] = totalUncertainty;

      // Compute background variance and uncertainty of the superob.
      const double meanHofXSquared =
        ufo::sumIntegralImagePatch(hofxSquaredArray, i, j, istep, jstep) / numObsInSuperObRegion;
      const double backgroundVariance = meanHofXSquared - std::pow(meanHofX, 2.0);
      const double backgroundUncertainty = std::sqrt(std::abs(backgroundVariance));
      backgroundUncertainty_[locSuperob] = backgroundUncertainty;
    }
  }
  oops::Log::trace() << "SuperObRadar::computeSuperOb finished" << std::endl;
}

void SuperObRadar::saveAuxiliaryVariables(const std::string & variableName) const {
  obsdb_.put_db("DiagnosticFlags", "UsedInSuperOb/" + variableName, usedInSuperOb_);
  obsdb_.put_db("TotalUncertainty", variableName, totalUncertainty_);
  obsdb_.put_db("BackgroundUncertainty", variableName, backgroundUncertainty_);
}

}  // namespace ufo
