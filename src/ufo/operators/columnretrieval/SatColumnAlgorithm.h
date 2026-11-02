/*
 * (C) Copyright 2017-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "ufo/operators/columnretrieval/ObsColumnRetrievalParameters.h"

namespace ufo {

/// Throw if the RetrievalAncillaryData variable \p varName read from the obs space
/// does not have the \p expected number of values (nlocs x number of levels).
void checkRetrievalDataSize(const std::string & varName, size_t actual, size_t expected);

/// Throw unless the operator simulates exactly one obs variable (\p nvars from the ObsVector).
void checkSingleSimulatedVariable(size_t nvars);

/// True if any averaging kernel value, any obs pressure vertex or the a priori term
/// at one observation location is missing; H(x) is then set to missing at that location.
bool hasMissingRetrievalData(const std::vector<double> & avgKernelObs,
                             const std::vector<double> & prsiObs,
                             double aprioriTerm);

/// Weights of the column retrieval for one observation location: the simulated
/// column is the sum over model layers j of weights[j] * profileModel[j].
///
/// Pressure arrays use the top-to-bottom convention: index 0 is the model top
/// and pressure increases with index.
///
/// \param avgKernelObs Averaging kernel [nlayersObs]
/// \param prsiObs      Obs pressure vertices [nlayersObs+1], top-to-bottom
/// \param prsiModel    Model pressure interfaces [nlayersModel+1], top-to-bottom
/// \param stretch      Vertex alignment option
/// \param tropopause   Tropopause pressure (Pa); layers above set AK=0 for tropospheric obs
/// \return             Weights [nlayersModel]
std::vector<double> columnRetrievalWeights(const std::vector<double> & avgKernelObs,
                                           const std::vector<double> & prsiObs,
                                           const std::vector<double> & prsiModel,
                                           StretchVertices stretch,
                                           double tropopause);

/// Non-linear column retrieval: compute H(x) for one observation location.
///
/// \param profileModel Model tracer profile [nlayersModel] (already unit-converted)
/// \return             Simulated column observation
/// See columnRetrievalWeights for the other arguments.
double simulateColumnOb(const std::vector<double> & avgKernelObs,
                        const std::vector<double> & prsiObs,
                        const std::vector<double> & prsiModel,
                        const std::vector<double> & profileModel,
                        StretchVertices stretch,
                        double tropopause);

/// Tangent linear of simulateColumnOb.
///
/// \param profileModelTL Perturbation of the model tracer profile [nlayersModel]
double simulateColumnObTL(const std::vector<double> & avgKernelObs,
                          const std::vector<double> & prsiObs,
                          const std::vector<double> & prsiModel,
                          const std::vector<double> & profileModelTL,
                          StretchVertices stretch,
                          double tropopause);

/// Adjoint of simulateColumnObTL.
/// Accumulates gradient into profileModelAD (must be pre-initialised by caller).
///
/// \param hofxAD         Adjoint forcing (scalar)
/// \param profileModelAD [inout] gradient w.r.t. model tracer profile [nlayersModel]
void simulateColumnObAD(const std::vector<double> & avgKernelObs,
                        const std::vector<double> & prsiObs,
                        const std::vector<double> & prsiModel,
                        double hofxAD,
                        std::vector<double> & profileModelAD,
                        StretchVertices stretch,
                        double tropopause);

}  // namespace ufo
