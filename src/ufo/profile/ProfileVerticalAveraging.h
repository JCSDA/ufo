/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEVERTICALAVERAGING_H_
#define UFO_PROFILE_PROFILEVERTICALAVERAGING_H_

#include <vector>

namespace ufo {
  namespace ProfileAveraging {
    enum class Method {Interpolation, Averaging};
  }

  /// \brief Profile vertical averaging onto model levels.
  ///
  /// This routine averages the observed values in a profile onto model levels.
  /// The vertical positions of observations are referred to as reported levels.
  ///
  /// Firstly, all valid reported level values are interpolated onto the adjacent model level
  /// boundaries. If there are no reported level values adjacent to a model level then the
  /// model level value is set to missing.
  ///
  /// If the parameter \p method is set to 'Interpolation' then no further action is taken;
  /// otherwise, averaging onto the midpoint of the model layer is performed as follows.
  ///
  /// If a model layer contains no reported levels then the mean of the values at the
  /// upper and lower model layer bounds is taken as the average value.
  /// If the model layer contains at least one reported level then the weighted mean of the
  /// values at the model layer boundaries and the reported levels is calculated.
  /// The distance between adjacent levels (including the model layer boundaries)
  /// is used as the weight in each case.
  /// If either of the model boundary values is missing then the layer is referred to as a
  /// partial layer. In this case a minimum fraction of distance between the upper and lower
  /// model levels must have been included in the weighted mean, otherwise the mean is
  /// set to the missing value.
  /// The minimum fraction is governed by the parameter \p DZFrac.
  ///
  /// The vector \p BigGap contains the thresholds for defining the distance between
  /// reported level coordinates as a 'big gap'.
  /// Interpolation to model levels is not performed for reported levels separted by a big gap.
  ///
  /// The vectors \p coordMin and \p coordMax, which are the coordinates of the highest and lowest
  /// valid values used in the weighted mean calculation, can optionally be recorded.
  /// These values are used in (e.g.) atmospheric profile temperature averaging.
  ///
  /// \param[in] flagsIn: flags on input data.
  /// \param[in] valuesIn: reported-level values.
  /// \param[in] coordIn: reported-level coordinates.
  /// \param[in] bigGap: maximum gap to be filled in.
  /// \param[in] coordOut: model layer boundaries.
  /// \param[in] DZFrac: minimum fraction for partially filled layers to be used.
  /// \param[in] method: interpolate to model layer boundaries or average over model layers?
  /// \param[out] flagsOut: flags on output data.
  /// \param[out] valuesOut: output data (averaged or interpolated).
  /// \param[out] numGaps: number of large gaps in profile.
  /// \param[out] coordMax: maximum coordinate of the values used in the model layer average.
  /// \param[out] coordMin: minimum coordinate of the values used in the model layer average.
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
                                std::vector <float> *coordMax = nullptr,
                                std::vector <float> *coordMin = nullptr);
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEVERTICALAVERAGING_H_
