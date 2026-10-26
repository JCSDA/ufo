/*
 * (C) Copyright 2017-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/columnretrieval/SatColumnAlgorithm.h"

#include <algorithm>
#include <numeric>
#include <string>
#include <tuple>

#include "eckit/exception/Exceptions.h"
#include "oops/util/missingValues.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {

namespace {

/// Adjust obs and model pressure vertices in place so they share the same
/// top/bottom boundary, depending on the stretch option.
void stretchVertices(StretchVertices stretch,
                     std::vector<double> & pobs, std::vector<double> & pmod) {
  if (stretch == StretchVertices::bottom || stretch == StretchVertices::topbottom) {
    const double surfP = std::max(pmod.back(), pobs.back());
    pmod.back() = surfP;
    pobs.back() = surfP;
  }
  if (stretch == StretchVertices::top || stretch == StretchVertices::topbottom) {
    const double topP = std::min(pmod.front(), pobs.front());
    pmod.front() = topP;
    pobs.front() = topP;
  }
}

}  // namespace

// -----------------------------------------------------------------------------

void checkRetrievalDataSize(const std::string & varName, size_t actual, size_t expected) {
  if (actual != expected) {
    throw eckit::BadValue("RetrievalAncillaryData/" + varName + " has " + std::to_string(actual) +
                          " values but " + std::to_string(expected) + " are expected; check "
                          "that nlayers_retrieval matches the obs file", Here());
  }
}

// -----------------------------------------------------------------------------

void checkSingleSimulatedVariable(size_t nvars) {
  if (nvars != 1) {
    throw eckit::BadValue("ColumnRetrieval simulates a single obs variable but " +
                          std::to_string(nvars) + " were requested", Here());
  }
}

// -----------------------------------------------------------------------------

bool hasMissingRetrievalData(const std::vector<double> & avgKernelObs,
                             const std::vector<double> & prsiObs,
                             double aprioriTerm) {
  const double missing = util::missingValue<double>();
  return aprioriTerm == missing ||
         std::find(avgKernelObs.begin(), avgKernelObs.end(), missing) != avgKernelObs.end() ||
         std::find(prsiObs.begin(), prsiObs.end(), missing) != prsiObs.end();
}

// -----------------------------------------------------------------------------

std::vector<double> columnRetrievalWeights(const std::vector<double> & avgKernelObs,
                                           const std::vector<double> & prsiObs,
                                           const std::vector<double> & prsiModel,
                                           StretchVertices stretch,
                                           double tropopause) {
  std::vector<double> pobs(prsiObs);
  std::vector<double> pmod(prsiModel);
  stretchVertices(stretch, pobs, pmod);

  const double massGrav = Constants::M_dryair * Constants::grav;

  std::vector<double> weights(pmod.size() - 1, 0.0);
  for (size_t k = 0; k < avgKernelObs.size(); ++k) {
    // Averaging kernel is zero above the tropopause
    const double ak = (pobs[k] < tropopause) ? 0.0 : avgKernelObs[k];

    int wiA, wiB;
    double wfA, wfB;
    std::tie(wiA, wfA) = PiecewiseLinearInterpolation::interpolationIndexAndWeight(pmod, pobs[k]);
    std::tie(wiB, wfB) = PiecewiseLinearInterpolation::interpolationIndexAndWeight(pmod,
                                                                                   pobs[k + 1]);

    // A missing obs pressure gives a missing index; callers skip such locations
    // (see hasMissingRetrievalData), so reaching this is a programming error.
    if (wiA < 0 || wiB < 0) {
      throw eckit::BadValue("Error: missing obs pressure vertex in column retrieval", Here());
    }

    // Check monotonically increasing pressure (top-to-bottom convention)
    if (pmod[wiA + 1] < pmod[wiA] || pmod[wiB + 1] < pmod[wiB]) {
      throw eckit::BadValue("Error: inverted pressure coordinate in geovals, "
                            "convention: top->bottom, increasing pressures", Here());
    }
    if (wiA > wiB) {
      throw eckit::BadValue("Error: inverted pressure coordinate in obs, "
                            "convention: top->bottom, increasing pressures", Here());
    }

    if (wiA < wiB) {
      // Multiple model levels span this obs layer
      weights[wiA] += ak * (pmod[wiA + 1] - pmod[wiA]) * wfA / massGrav;
      for (int j = wiA + 1; j < wiB; ++j) {
        weights[j] += ak * (pmod[j + 1] - pmod[j]) / massGrav;
      }
      weights[wiB] += ak * (pmod[wiB + 1] - pmod[wiB]) * (1.0 - wfB) / massGrav;
    } else {
      // wiA == wiB: obs layer is within one model level
      weights[wiA] += ak * (pmod[wiA + 1] - pmod[wiA]) * (wfA - wfB) / massGrav;
    }
  }
  return weights;
}

// -----------------------------------------------------------------------------

double simulateColumnOb(const std::vector<double> & avgKernelObs,
                        const std::vector<double> & prsiObs,
                        const std::vector<double> & prsiModel,
                        const std::vector<double> & profileModel,
                        StretchVertices stretch,
                        double tropopause) {
  const std::vector<double> weights =
      columnRetrievalWeights(avgKernelObs, prsiObs, prsiModel, stretch, tropopause);
  return std::inner_product(weights.begin(), weights.end(), profileModel.begin(), 0.0);
}

// -----------------------------------------------------------------------------

double simulateColumnObTL(const std::vector<double> & avgKernelObs,
                          const std::vector<double> & prsiObs,
                          const std::vector<double> & prsiModel,
                          const std::vector<double> & profileModelTL,
                          StretchVertices stretch,
                          double tropopause) {
  const std::vector<double> weights =
      columnRetrievalWeights(avgKernelObs, prsiObs, prsiModel, stretch, tropopause);
  return std::inner_product(weights.begin(), weights.end(), profileModelTL.begin(), 0.0);
}

// -----------------------------------------------------------------------------

void simulateColumnObAD(const std::vector<double> & avgKernelObs,
                        const std::vector<double> & prsiObs,
                        const std::vector<double> & prsiModel,
                        double hofxAD,
                        std::vector<double> & profileModelAD,
                        StretchVertices stretch,
                        double tropopause) {
  const std::vector<double> weights =
      columnRetrievalWeights(avgKernelObs, prsiObs, prsiModel, stretch, tropopause);
  for (size_t j = 0; j < weights.size(); ++j) {
    profileModelAD[j] += hofxAD * weights[j];
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
