/*
 * (C) Crown copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/RadarScanEdgeFlag.h"

#include <Eigen/Dense>

#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/DiagnosticFlag.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/IntegralImage.h"

namespace ufo {

static ObsFunctionMaker<RadarScanEdgeFlag>
makerRadarScanEdgeFlag_("RadarScanEdgeFlag");

// -----------------------------------------------------------------------------

RadarScanEdgeFlag::RadarScanEdgeFlag(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Validate and deserialize options
  options_.validateAndDeserialize(conf);

  // Obtain variables associated with the `where` options.
  invars_ += getAllWhereVariables(options_.where);
}

// -----------------------------------------------------------------------------

RadarScanEdgeFlag::~RadarScanEdgeFlag() {}

// -----------------------------------------------------------------------------

void RadarScanEdgeFlag::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<int> & out) const {
  oops::Log::trace() << "RadarScanEdgeFlag::compute started" << std::endl;

  const ioda::ObsSpace & obsdb = in.obsspace();

  // Ensure observations have been grouped into records.
  if (obsdb.obs_group_vars().empty())
    throw eckit::UserError("Group variables configuration is empty", Here());

  // Number of locations.
  const size_t nlocs = obsdb.nlocs();

  // Correspondence between record numbers and locations in the data sample.
  const std::vector<std::size_t> & recnums = obsdb.recidx_all_recnums();

  // Number of records.
  const std::size_t nrecs = recnums.size();

  // Get observed radial velocity.
  std::vector<float> radialVelocity(nlocs);
  obsdb.get_db("ObsValue", "radialVelocity", radialVelocity);

  // Get metadata.
  std::vector<int> numGates(nlocs);
  std::vector<int> numBeams(nlocs);
  std::vector<float> gateWidth(nlocs);
  std::vector<float> beamWidth(nlocs);
  std::vector<float> minGateRange(nlocs);
  std::vector<float> gateRange(nlocs);
  std::vector<float> beamAzimuthAngle(nlocs);
  obsdb.get_db("MetaData", "numberOfGates", numGates);
  obsdb.get_db("MetaData", "numberOfBeams", numBeams);
  obsdb.get_db("MetaData", "gateWidth", gateWidth);
  obsdb.get_db("MetaData", "beamWidth", beamWidth);
  obsdb.get_db("MetaData", "minGateRange", minGateRange);
  obsdb.get_db("MetaData", "gateRange", gateRange);
  obsdb.get_db("MetaData", "beamAzimuthAngle", beamAzimuthAngle);

  // Vector of locations that pass the 'where' clause in the sample
  // (all true if there is no where clause).
  const std::vector<bool> apply = processWhere(options_.where, in, options_.whereOperator);

  // Loop over records (i.e. radar scans).
  for (std::size_t jrec = 0; jrec < nrecs; ++jrec) {
    // Get locations for this record.
    const std::vector<std::size_t> & locs = obsdb.recidx_vector(recnums[jrec]);

    // Obtain metadata for this scan.
    // All scan-related metadata values are the same, so can just take the first location
    // in the record regardless of QC decisions.
    const std::size_t firstLocInScan = locs.front();
    const float scanMinGateRange = minGateRange[firstLocInScan];
    const float scanGateWidth = gateWidth[firstLocInScan];
    const float scanBeamWidth = beamWidth[firstLocInScan];

    // Fill arrays to be processed with the integral image technique.
    // The arrays are extended at the edges.
    // todo(ctgh): Check the values used to extend the arrays here
    // and in the SuperObRadarTemplate class.
    const int imax = numGates[firstLocInScan] + 2;
    const int jmax = numBeams[firstLocInScan] + 3;
    Eigen::ArrayXXd obsArray = Eigen::ArrayXXd::Zero(imax, jmax);
    Eigen::ArrayXXi countArray = Eigen::ArrayXXi::Zero(imax, jmax);
    Eigen::ArrayXXi location = Eigen::ArrayXXi::Zero(imax, jmax);

    for (int jloc : locs) {
      // Do not proceed if this location has been excluded by the `where` clause.
      if (!apply[jloc]) {
        continue;
      }

      // Determine coordinates of this observation in the 2D superobbing arrays.
      const int rangeIndex =
        std::round((gateRange[jloc] - scanMinGateRange +
                    0.5 * scanGateWidth) / scanGateWidth);
      const int azimIndex =
        std::round((beamAzimuthAngle[jloc] + 0.5 * scanBeamWidth) / scanBeamWidth) + 1;

      obsArray(rangeIndex, azimIndex) = radialVelocity[jloc];
      countArray(rangeIndex, azimIndex) = 1;
      location(rangeIndex, azimIndex) = jloc;
    }

    // Wrap arrays in azimuthal direction. j-indices 0, 1 and (jmax - 1) are filled
    // from corresponding entries (jmax - 3), (jmax - 2) and 2.
    for (size_t i = 0; i < imax; ++i) {
      obsArray(i, 0) = obsArray(i, jmax - 3);
      obsArray(i, 1) = obsArray(i, jmax - 2);
      obsArray(i, jmax - 1) = obsArray(i, 2);
      countArray(i, 0) = countArray(i, jmax - 3);
      countArray(i, 1) = countArray(i, jmax - 2);
      // In OPS, the wrapping between obsArray and countArray
      // is not consistent.
      if (options_.opsCompatibilityMode) {
        countArray(i, jmax - 1) = countArray(i, 1);
      } else {
        countArray(i, jmax - 1) = countArray(i, 2);
      }
    }

    // Compute integral image of observations and counts.
    Eigen::ArrayXXd integratedObsArray = obsArray;
    Eigen::ArrayXXi integratedCountArray = countArray;
    ufo::createIntegralImage(integratedObsArray);
    ufo::createIntegralImage(integratedCountArray);

    // Output arrays.
    Eigen::ArrayXXi cleanArray = Eigen::ArrayXXi::Zero(imax, jmax);
    Eigen::ArrayXXi doubleCleanArray = Eigen::ArrayXXi::Zero(imax, jmax);
    Eigen::ArrayXXd LaplaceArray = Eigen::ArrayXXd::Zero(imax, jmax);
    Eigen::ArrayXXd flagArray = Eigen::ArrayXXd::Zero(imax, jmax);

    // Filter parameters.
    const double thresholdLaplaceFilter = options_.thresholdLaplaceFilter.value();
    const int thresholdCleanFilter = options_.thresholdCleanFilter.value();
    const int thresholdDoubleCleanFilter = options_.thresholdDoubleCleanFilter.value();

    // Apply Laplace filter?
    const bool applyLaplaceFilter = thresholdLaplaceFilter > 0.0;
    // Apply double-clean filter?
    const bool applyDoubleCleanFilter = thresholdCleanFilter > 0 ||
      thresholdDoubleCleanFilter > 0;

    for (size_t i = 2; i < imax - 1; ++i) {
      for (size_t j = 2; j < jmax - 1; ++j) {
        // Number of grid boxes containing an observation.
        const int integratedCount = ufo::sumIntegralImagePatch(integratedCountArray,
                                                               i, j, 1, 1);
        // Integrated observation value in the same patch.
        const double integratedObs = ufo::sumIntegralImagePatch(integratedObsArray,
                                                                i, j, 1, 1);
        // Laplace filter array.
        // Addition of 1.0e-6 in the denominator prevents division by zero.
        LaplaceArray(i, j) =
          std::abs((integratedObs - integratedCount * obsArray(i, j)) /
                   (integratedCount + 1.0e-6));
        cleanArray(i, j) = integratedCount;
        // Double integral of values.
        doubleCleanArray(i - 1, j) = ufo::sumIntegralImagePatch(integratedCountArray,
                                                                i, j, 1, 0);
      }
    }

    // Apply flags.
    if (applyLaplaceFilter) {
      for (size_t i = 0; i < imax; ++i) {
        for (size_t j = 0; j < jmax; ++j) {
          if (flagArray(i, j) == 0 &&
              LaplaceArray(i, j) >= thresholdLaplaceFilter) {
            flagArray(i, j) = 1;
          }
        }
      }
    }
    if (applyDoubleCleanFilter) {
      for (size_t i = 0; i < imax; ++i) {
        for (size_t j = 0; j < jmax; ++j) {
          if (flagArray(i, j) == 0 &&
              cleanArray(i, j) < thresholdCleanFilter) {
            flagArray(i, j) = 2;
          } else if (flagArray(i, j) == 0 &&
                     doubleCleanArray(i, j) < thresholdDoubleCleanFilter) {
            flagArray(i, j) = 3;
          }
        }
      }
    }

    // Only flag grid boxes with at least one observation.
    flagArray *= countArray.cast<double>();

    // Transfer decisions to output variable.
    for (size_t i = 2; i < imax - 1; ++i) {
      for (size_t j = 2; j < jmax - 1; ++j) {
        const int loc = location(i, j);
        if (apply[loc] && flagArray(i, j) > 0) {
          out[0][loc] = flagArray(i, j);
        }
      }
    }
  }  // loop over records

  oops::Log::trace() << "RadarScanEdgeFlag::compute finished" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & RadarScanEdgeFlag::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
