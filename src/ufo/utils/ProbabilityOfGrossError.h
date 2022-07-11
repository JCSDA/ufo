/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_PROBABILITYOFGROSSERROR_H_
#define UFO_UTILS_PROBABILITYOFGROSSERROR_H_

#include <algorithm>
#include <string>
#include <vector>

#include "oops/util/missingValues.h"

#include "ufo/utils/metoffice/MetOfficeQCFlags.h"
#include "ufo/utils/ProbabilityOfGrossErrorParameters.h"

namespace ufo {
  /// \brief Bayesian update of probability of gross error (PGE)
  /// \details Update PGE across locations according to (obsVal-bkgVal), obsErr, and bkgErr,
  /// and update flags to say BG check performed, and whether obs rejected or not.
  /// This routine can process both single observations and observations on profile levels.
  /// In the vector case, the variance is assumed to be isotropic
  /// and obsErr and bkgErr give the error in a single vector component.
  ///
  /// \param[in] options: Configurable parameters that govern the operation of this routine.
  /// \param[in] obsVal: Observation values.
  /// \param[in] obsErr: Observation errors.
  /// \param[in] bkgVal: Background values.
  /// \param[in] bkgErr: Background errors.
  /// \param[in] PdBad: Probability density for 'bad' observations.
  /// \param[in] ModelLevels: Have the data been averaged onto model levels?
  /// \param[inout] flags: QC flags.
  /// \param[inout] PGE: Probability of gross error.
  /// \param[in] ErrVarMax: (Optional) Maximum error variance.
  /// \param[in] obsVal2: (Optional) Second component of 2D observation values.
  /// \param[in] bkgVal2: (Optional) Second component of 2D background values.
  /// \param[inout] TotalPd: (Optional) Total (combined) probability distribution.

  void BayesianPGEUpdate(const ProbabilityOfGrossErrorParameters &options,
                         const std::vector<float> &obsVal,
                         const std::vector<float> &obsErr,
                         const std::vector<float> &bkgVal,
                         const std::vector<float> &bkgErr,
                         const std::vector<float> &PdBad,
                         const bool ModelLevels,
                         std::vector<int> &flags,
                         std::vector<float> &PGE,
                         float ErrVarMax = -1,
                         const std::vector<float> *obsVal2 = nullptr,
                         const std::vector<float> *bkgVal2 = nullptr,
                         std::vector<float> *TotalPd = nullptr);
}  // namespace ufo

#endif  // UFO_UTILS_PROBABILITYOFGROSSERROR_H_
