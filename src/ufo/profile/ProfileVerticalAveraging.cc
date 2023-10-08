/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <algorithm>
#include <functional>
#include <utility>

#include "eckit/exception/Exceptions.h"

#include "oops/util/CompareNVectors.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/PropertiesOfNVectors.h"

#include "ufo/profile/ProfileVerticalAveraging.h"

#include "ufo/utils/metoffice/MetOfficeQCFlags.h"

namespace ufo {
  void calculateVerticalAverage(const std::vector <int> &flagsIn,
                                const std::vector <float> &valuesIn,
                                const std::vector <float> &coordIn,
                                const std::vector <float> &bigGap,
                                const std::vector <float> &coordOut,
                                float DZFrac,
                                ProfileAveraging::Method method,
                                std::vector <int> &flagsOut,
                                std::vector <float> &valuesOut,
                                int &numGaps,
                                std::vector <float> *coordMax,
                                std::vector <float> *coordMin)
  {
    const float missingValueFloat = util::missingValue<float>();
    const size_t numRepLev = coordIn.size();
    const size_t numInterp = coordOut.size();
    const size_t numOut = method == ProfileAveraging::Method::Interpolation ?
      numInterp : numInterp - 1;  // Number of output levels.

    // Coordinates in ascending order?
    const bool Ascending = coordOut[0] < coordOut[numInterp - 1];

    // Make local copies of the vertical coordinates in order to allow them to be reversed.
    std::vector <float> ZIn = coordIn;
    std::vector <float> ZOut = coordOut;

    // Multiply coordinates by -1 if ascending.
    if (Ascending) {
      std::transform(ZIn.begin(), ZIn.end(), ZIn.begin(),
                     std::bind(std::multiplies<float>(), std::placeholders::_1, -1));
      std::transform(ZOut.begin(), ZOut.end(), ZOut.begin(),
                     std::bind(std::multiplies<float>(), std::placeholders::_1, -1));
    }

    // Initialise coordMin and coordMax, which record the minimum and maximum
    // coordinates of the values used in the model layer average.
    if (coordMin) coordMin->assign(numOut, missingValueFloat);
    if (coordMax) coordMax->assign(numOut, missingValueFloat);

    // Note observation levels with data that are useful for averaging.

    // Interpolated/averaged values.
    std::vector <float> valuesInterp(numInterp, missingValueFloat);
    // Indices of useful data.
    std::vector <size_t> idxUsefulLevels;
    // True if there is a big gap relative to the previous level.
    std::vector <bool> bigGapWithPreviousLevel(numRepLev + 1, false);
    // Previous value of JLev.
    size_t JLevP = 0;
    for (size_t JLev = 0; JLev < numRepLev; ++JLev) {
      // Ignore levels with missing data or pre-existing rejection flags.
      if (valuesIn[JLev] == missingValueFloat ||
          flagsIn[JLev] & ufo::MetOfficeQCFlags::Elem::FinalRejectFlag)
        continue;
      // Ignore duplicate levels.
      if (idxUsefulLevels.size() > 0 && ZIn[JLevP] == ZIn[JLev]) continue;
      // Set bigGapWithPreviousLevel in two cases:
      if (idxUsefulLevels.size() == 0) {
        // Cannot interpolate before first level.
        bigGapWithPreviousLevel[idxUsefulLevels.size()] = true;
      } else if (ZIn[JLevP] - ZIn[JLev] > std::max(bigGap[JLevP], bigGap[JLev])) {
        numGaps++;
        // Big gap from previous useful level.
        bigGapWithPreviousLevel[idxUsefulLevels.size()] = true;
      }
      // Record index of this level.
      idxUsefulLevels.push_back(JLev);
      JLevP = JLev;
    }

    // Loop over model levels, interpolating from useful reported levels.

    std::vector <int> flagsInterp(numInterp, 0);  // Interpolation flags.
    // Counter over reported levels which is incremented inside the loop over model levels.
    size_t JLev = 0;
    for (size_t MLev = 0; MLev < numInterp; ++MLev) {
      JLevP = JLev;  // Previous value of JLev.
      const double ZMLev = ZOut[MLev];  // Coordinate value at current output level.
      // Increment JLev until the associated coordinate is less than ZMLev
      // or the number of useful levels is reached.
      while (JLev != idxUsefulLevels.size() && ZIn[idxUsefulLevels[JLev]] >= ZMLev)
        ++JLev;

      // JLev is the first reported level above model level MLev.
      // JLev = 0 => ZMLev below bottom level.
      // JLev = idxUsefulLevels.size() => ZMLev above top level.
      if (JLev == 0 || JLev == idxUsefulLevels.size()) {
        // Top or bottom - do not allow extrapolation.
      } else if (!bigGapWithPreviousLevel[JLev]) {
        // Interpolate to model layer boundaries.
        const size_t J1 = idxUsefulLevels[JLev - 1];
        const size_t J2 = idxUsefulLevels[JLev];
        // Weight given to upper layer in interpolation.
        const double Interp_factor = (ZMLev - ZIn[J1]) / (ZIn[J2] - ZIn[J1]);
        // Interpolate reported level values onto model levels.
        valuesInterp[MLev] = valuesIn[J1] + (valuesIn[J2] - valuesIn[J1]) * Interp_factor;
        flagsInterp[MLev] = flagsIn[J1] | flagsIn[J2];
      }

      // Loop if interpolating or at lowest model level.
      if (method == ProfileAveraging::Method::Interpolation || MLev == 0)
        continue;

      // Average over model layers.

      if (JLevP == JLev) {
        // Model layer contains no reported levels.
        if (valuesInterp[MLev - 1] != missingValueFloat &&
            valuesInterp[MLev] != missingValueFloat) {
          // Use mean of values at model layer bounds.
          valuesInterp[MLev - 1] = (valuesInterp[MLev - 1] + valuesInterp[MLev]) * 0.5;
          flagsInterp[MLev - 1] |= flagsInterp[MLev];
          if (coordMax) (*coordMax)[MLev - 1] = ZOut[MLev - 1];
          if (coordMin) (*coordMin)[MLev - 1] = ZOut[MLev];
        }
      } else {
        // Model layer contains at least one reported level.
        // Compute mean of values over the layer,
        // weighting by the difference between coordinates in each segment.

        // Sum of DZ (coordinate differences) within model layer.
        // Used as the denominator in the weighted mean.
        // This is only different from unity if one or both of the upper or lower model bounds
        // are missing.
        float DZSum = 0.0;

        // Contribution from lowest segment, i.e. gap between
        // lower bound of model layer and lower reported level coordinate.
        size_t J2 = idxUsefulLevels[JLevP];
        if (valuesInterp[MLev - 1] != missingValueFloat) {
          const double DZ = ZOut[MLev - 1] - ZIn[J2];
          valuesInterp[MLev - 1] = (valuesInterp[MLev - 1] + valuesIn[J2]) * DZ;
          DZSum += DZ;
          if (coordMax) (*coordMax)[MLev - 1] = ZOut[MLev - 1];
        } else {
          valuesInterp[MLev - 1] = 0.0;
          if (coordMax) (*coordMax)[MLev - 1] = ZIn[J2];
        }

        // Contribution from highest segment, i.e. gap between
        // highest reported level coordinate and upper bound of model layer.
        size_t J1 = idxUsefulLevels[JLev - 1];
        if (valuesInterp[MLev] != missingValueFloat) {
          const double DZ = ZIn[J1] - ZOut[MLev];
          valuesInterp[MLev - 1] += (valuesInterp[MLev] + valuesIn[J1]) * DZ;
          DZSum += DZ;
          flagsInterp[MLev - 1] |= flagsInterp[MLev];
          if (coordMin) (*coordMin)[MLev - 1] = ZOut[MLev];
        } else {
          if (coordMin) (*coordMin)[MLev - 1] = ZIn[J1];
        }

        // Sum contributions from intermediate segments, i.e. gaps between
        // adjacent reported level coordinates.
        for (size_t J = JLevP + 1; J < JLev; ++J) {
          J1 = idxUsefulLevels[J - 1];
          J2 = idxUsefulLevels[J];
          const double DZ = ZIn[J1] - ZIn[J2];
          valuesInterp[MLev - 1] += (valuesIn[J1] + valuesIn[J2]) * DZ;
          DZSum += DZ;
          flagsInterp[MLev - 1] |= flagsIn[J1];
        }
        J1 = idxUsefulLevels[JLev - 1];
        flagsInterp[MLev - 1] |= flagsIn[J1];

        // Calculate layer mean values.

        // Difference in coordinate across model layer.
        const double DZMod = ZOut[MLev - 1] - ZOut[MLev];

        // Only process layers that are full by more than a certain amount.
        if (DZSum > DZFrac * DZMod) {
          // Set partial layer flag if the layer is between DZFrac and 99.5% full.
          if (DZSum <= 0.995 * DZMod)
            flagsInterp[MLev - 1] |= ufo::MetOfficeQCFlags::Profile::PartialLayerFlag;
          // Divide the weighted sum by sum of the weights
          // (the factor of 0.5 ensures correct normalisation).
          valuesInterp[MLev - 1] *= 0.5 / DZSum;
        } else {
          valuesInterp[MLev - 1] = missingValueFloat;
        }
      }

      // Negate and swap coordMax and coordMin for this model level if the coordinates are ascending
      // and the values are not missing.
      if (Ascending && coordMin && coordMax &&
          (*coordMax)[MLev - 1] != missingValueFloat &&
          (*coordMin)[MLev - 1] != missingValueFloat) {
        (*coordMax)[MLev - 1] *= -1;
        (*coordMin)[MLev - 1] *= -1;
        std::swap((*coordMax)[MLev - 1], (*coordMin)[MLev - 1]);
      }
    }

    // Fill output arrays.
    flagsOut = {flagsInterp.begin(), flagsInterp.begin() + numOut};
    valuesOut = {valuesInterp.begin(), valuesInterp.begin() + numOut};
  }
}  // namespace ufo
