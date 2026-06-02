/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_MODELBESTFITPRESSURE_H_
#define UFO_FILTERS_MODELBESTFITPRESSURE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit {
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// \brief Parameters controlling the operation of the ModelBestFitPressure filter.
class ModelBestFitPressureParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ModelBestFitPressureParameters, FilterParametersBase)

 public:
  /// Name of the observation pressure variable.
  oops::RequiredParameter<Variable> obs_pressure{"observation pressure", this};
  /// Name of the model pressure variable.
  oops::RequiredParameter<Variable> model_pressure{"model pressure", this};
  /// Minimum allowed pressure region. Model levels will not be considered if their
  /// pressure is below this value.
  oops::Parameter<float> top_pressure{"top pressure", 10000, this};
  /// Winds within this pressure-range of the best-fit pressure will be checked for
  /// consistency within this range of the minimum speed.
  oops::Parameter<float> pressure_band_half_width{"pressure band half-width", 10000, this};
  /// Maximum vector difference allowed, for calculating constraint.
  oops::Parameter<float> upper_vector_diff{"upper vector diff", 4, this};
  /// Minimum vector difference allowed, for calculating constraint.
  oops::Parameter<float> lower_vector_diff{"lower vector diff", 2, this};
  /// Tolerance for vec_diff comparison.
  /// Used when calculating bestfit pressure using parabolic fit.
  oops::Parameter<float> tolerance_vector_diff{"tolerance vector diff", 1.0e-8, this};
  /// Tolerance for pressure comparison. Used for calculating bestfit winds when comparing
  /// pressure with bestfit pressure. Only used if calculate_best_fit_winds is true.
  oops::Parameter<float> tolerance_pressure{"tolerance pressure", 0.01, this};
  /// To calculate bestfit eastward/northward winds by linear interpolation.
  oops::Parameter<bool> calculate_best_fit_winds{"calculate bestfit winds", false, this};
};

/// \brief A filter to calculate the best fit pressure and if the pressure is well constrained.
/// Also can calculate the best fit eastward/northward winds.
///
/// \details The model best-fit pressure is defined as the model pressure (Pa) with the
///  smallest vector difference between the AMV and model background wind, but
///  additionally is not allowed to be above top pressure (top_pressure) (can reasonably
///  expect that AMVs should not be above this level). Vertical interpolation is performed
///  between model levels to find the minimum vector difference.
///
///  Checking if the pressure is well-constrained:
///  1. Remove any winds where the minimum vector difference between the AMV eastward and
///     northward winds and the background column u and v is greater than upper_vector_diff.
///     This check aims to remove cases where there is no good agreement between the AMV
///     and the winds at any level in the background wind column.
///  2. Remove any winds where the vector difference is less than the minimum vector
///     difference + lower_vector_diff outside of a band +/- pressure_band_half_width from the
///     best-fit pressure level. This aims to catch cases where there are secondary minima
///     or very broad minima. In both cases the best-fit pressure is not well constrained.
///  The default parameter values were chosen by eye-balling vector difference profiles and together
///  remove just over half the winds.
///
/// See ModelBestFitPressureParameters for the documentation of
/// the parameters controlling this filter.
class ModelBestFitPressure : public FilterBase,
                        private util::ObjectCounter<ModelBestFitPressure> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ModelBestFitPressureParameters Parameters_;

  static const std::string classname() {return "ufo::ModelBestFitPressure";}

  ModelBestFitPressure(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ModelBestFitPressure();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::pass; }

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_MODELBESTFITPRESSURE_H_
