/*
 * (C) Copyright 2026 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_DUPLICATETHINNING_H_
#define UFO_FILTERS_DUPLICATETHINNING_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

class ObsAccessor;

/// \brief Parameters controlling the DuplicateThinning filter.
class DuplicateThinningParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(DuplicateThinningParameters, FilterParametersBase)

 public:
  /// Vector of names of the variable used to identify duplicates.
  oops::Parameter<std::vector<ufo::Variable>> variableNames{"variable names",
                                           {ufo::Variable("MetaData/latitude")}, this};

  /// Vector of tolerances for considering two observations as duplicates.
  /// Two observations are considered duplicates if the absolute difference
  ///  of all the specified variables are less than or equal to the corresponding tolerance.
  /// Applied exclusively only for variables with float value.
  oops::Parameter<std::vector<float>> tolerance{"tolerance", {0.5}, this};

  oops::OptionalParameter<util::DateTime> analysisTime{"analysis_time",
      "The default value is the start of the time window.", this};

  /// Observations with dateTime within the tolerance of analysis_time are kept
  /// and will be evaluated as candidates for temporal-thinning.
  /// e.g. PT90M for 90 minutes. Default value is PT0H, i.e. no skip window.
  oops::Parameter<util::Duration> analysisTimeTolerance{
      "analysis_time_tolerance", util::Duration("PT0H"), this};

  /// Minimum time gap required between two observations retained in the same
  /// group. Observation is thinned if its dateTime is within
  /// `min_spacing` of any already-retained observation in the same group.
  /// e.g. PT15M for 15 minutes. Default value is set to PT6H.
  oops::Parameter<util::Duration> minSpacing{
      "min_spacing", util::Duration("PT6H"), this};

  /// Tie-breaker for the min-spacing pass: when two observations are
  /// equidistant from `analysis_time`, choose the one whose dateTime is
  /// `"after"` (default) or `"before"` the analysis time.
  oops::Parameter<std::string> equidistantTimeSelection{
      "equidistant_time_selection", "after", this};

 private:
  // This filter doesn't take any extra parameters.
};

/// \brief A filter that, by default, performs the `duplicateThinning`
/// action on observations selected by the variables and tolerances specified
/// i.e. resets the QC flags of these observations to `thinned`
class DuplicateThinning : public FilterBase,
                   private util::ObjectCounter<DuplicateThinning> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef DuplicateThinningParameters Parameters_;

  static const std::string classname() {return "ufo::DuplicateThinning";}

  DuplicateThinning(ioda::ObsSpace &, const Parameters_ &,
                    ioda::ObsDataVector<int> &,
                    ioda::ObsDataVector<float> &);

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;

  int qcFlag() const override { return QCflags::thinned; }
  ObsAccessor createObsAccessor() const;
  Parameters_ parameters_;
  ufo::Variables testVars_;
  std::vector<float> toleranceVec_;
  int64_t timeToleranceSec_;
  int64_t minSpacingSec_;
  bool preferAfter_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_DUPLICATETHINNING_H_
