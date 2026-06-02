/*
 * (C) Copyright 2022 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_PROFILEAVERAGEOBSTOMODLEVELS_H_
#define UFO_FILTERS_PROFILEAVERAGEOBSTOMODLEVELS_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the ProfileAverageObsToModLevels filter.
class ProfileAverageObsToModLevelsParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ProfileAverageObsToModLevelsParameters, FilterParametersBase)

 public:
  /// Model vertical coordinate, e.g. HofX/depthBelowWaterSurface (after VertInterp and
  /// ProfileAverage  have been applied)
  oops::RequiredParameter<Variable> model_vertical_coordinate
    {"model vertical coordinate",
     this};

  /// Observation vertical coordinate, e.g. DerivedObsValue/depthBelowWaterSurface
  oops::RequiredParameter<Variable> observation_vertical_coordinate
    {"observation vertical coordinate",
     this};
};

/// ProfileAverageObsToModLevels: average the observation values on to model levels.
///
/// For each of the filter variables given, this filter computes the model-level average increment:
/// the mean of all where-included, QC-passing, non-missing observation-minus-background values
/// that fall within the range of that model level (bounded by the mid-points between it and the
/// adjacent model level above and below). Each average increment is added to the H(x) value at
/// that model level, and the result is written to the corresponding location of the
/// DerivedObsValue's extended space (its original space has ObsValue copied into it).
///
/// In order for this filter to work correctly the ObsSpace must have been extended as in
/// the following yaml snippet:
///
/// - obs space:
///    extension:
///      allocate companion records with length: 71
///
/// (where 71 can be replaced by the length of the vertical coordinate GeoVaL).
/// AND the model vertical coordinate MUST NOT have all zeros in the extended space - e.g. by
/// making sure to apply `ProfileAverage` obs operator to the observation vertical coordinate.
///
/// Requires the following be specified in .yaml, under
///
/// obs filters:
/// - filter: Average Observations to Model Levels
///
///   * model vertical coordinate: e.g. HofX/depthBelowWaterSurface
///   * observation vertical coordinate: e.g. DerivedObsValue/depthBelowWaterSurface
///

class ProfileAverageObsToModLevels : public FilterBase,
                        private util::ObjectCounter<ProfileAverageObsToModLevels> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ProfileAverageObsToModLevelsParameters Parameters_;

  static const std::string classname() {return "ufo::ProfileAverageObsToModLevels";}

  ProfileAverageObsToModLevels(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ProfileAverageObsToModLevels();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::fguess;}
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PROFILEAVERAGEOBSTOMODLEVELS_H_
