/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_MODELOBTHRESHOLD_H_
#define UFO_FILTERS_MODELOBTHRESHOLD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
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

enum class ThresholdType {
  MIN, MAX
};

struct ThresholdTypeParameterTraitsHelper {
  typedef ThresholdType EnumType;
  static constexpr char enumTypeName[] = "ThresholdType";
  static constexpr util::NamedEnumerator<ThresholdType> namedValues[] = {
    { ThresholdType::MIN, "min" },
    { ThresholdType::MAX, "max" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::ThresholdType> :
    public EnumParameterTraits<ufo::ThresholdTypeParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// \brief Parameters controlling the operation of the ModelObThreshold filter.
class ModelObThresholdParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ModelObThresholdParameters, FilterParametersBase)

 public:
  /// Name of the model profile variable (GeoVaLs)
  oops::RequiredParameter<Variable> model_profile{"model profile", this};
  /// Name of the model vertical coordinate variable (GeoVal)
  oops::RequiredParameter<Variable> model_vcoord{"model vertical coordinate", this};
  /// Name of the observation height variable to interpolate to
  oops::RequiredParameter<Variable> obs_height{"observation height", this};
  /// Vector of threshold values
  oops::RequiredParameter<std::vector<double>> thresholds{"thresholds", this};
  /// Vector of vertical coordinates corresponding to vector of thresholds
  oops::RequiredParameter<std::vector<double>> coord_vals{"coordinate values", this};
  /// Threshold type: min or max
  oops::RequiredParameter<ThresholdType> threshold_type{"threshold type", this};
};

// -----------------------------------------------------------------------------

/// \brief A filter that interpolates a model profile (GeoVaL) and a height-dependent
/// threshold to the observation location and flags observations which are outside the
/// specified limit.
///
/// See ModelObThresholdParameters for the documentation of the parameters controlling this filter.
class ModelObThreshold : public FilterBase,
                       private util::ObjectCounter<ModelObThreshold> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ModelObThresholdParameters Parameters_;

  static const std::string classname() {return "ufo::ModelObThreshold";}

  ModelObThreshold(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ModelObThreshold();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::modelobthresh;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_MODELOBTHRESHOLD_H_
