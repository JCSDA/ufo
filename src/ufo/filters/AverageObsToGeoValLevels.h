/*
 * (C) Copyright 2024 NASA GMAO
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_AVERAGEOBSTOGEOVALLEVELS_H_
#define UFO_FILTERS_AVERAGEOBSTOGEOVALLEVELS_H_

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

/// Parameters controlling the operation of the AverageObsToGeoValLevels filter.
class AverageObsToGeoValLevelsParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(AverageObsToGeoValLevelsParameters, FilterParametersBase)

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

/// AverageObsToGeoValLevels: average the observation values on to model levels in the GeoVals.
///
/// For each of the filter variables given, this filter computes the model-level average values
/// that fall within the range of that model level (bounded by the mid-points between it and the
/// adjacent model level above and below). The result is written to the corresponding location
/// of the DerivedObsValue's extended space (its original space has ObsValue copied into it).
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
/// - filter: Average Observations to GeoVals Levels
///
///   * model vertical coordinate: e.g. pressure
///   * observation vertical coordinate: e.g. MetaData/height
///

class AverageObsToGeoValLevels : public FilterBase,
                        private util::ObjectCounter<AverageObsToGeoValLevels> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef AverageObsToGeoValLevelsParameters Parameters_;

  static const std::string classname() {return "ufo::AverageObsToGeoValLevels";}

  AverageObsToGeoValLevels(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~AverageObsToGeoValLevels();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::fguess;}
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_AVERAGEOBSTOGEOVALLEVELS_H_
