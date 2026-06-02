/*
 * (C) Crown copyright 2025 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_ENSEMBLESTATISTICS_H_
#define UFO_FILTERS_ENSEMBLESTATISTICS_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "ufo/filters/ObsProcessorBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsFilterParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Type of an ensemble statistic.
enum class EnsembleStatistic {
  /// \brief Ensemble mean of H(x) vectors.
  HOFX_MEAN,
  /// \brief Ensemble spread of H(x) vectors.
  HOFX_STDDEV
};

struct EnsembleStatisticParameterTraitsHelper {
  typedef EnsembleStatistic EnumType;
  static constexpr char enumTypeName[] = "EnsembleStatistic";
  static constexpr util::NamedEnumerator<EnsembleStatistic> namedValues[] = {
    { EnsembleStatistic::HOFX_MEAN, "HofXEnsembleMean" },
    { EnsembleStatistic::HOFX_STDDEV, "HofXEnsembleStdDev" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::EnsembleStatistic> :
    public EnumParameterTraits<ufo::EnsembleStatisticParameterTraitsHelper>
{};

}  // namespace oops

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// \brief Parameters controlling the operation of the EnsembleStatistics filter.
class EnsembleStatisticsParameters : public ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(EnsembleStatisticsParameters, ObsFilterParametersBase)

 public:
  /// Simulated variables (and channels) for which ensemble statistics will be calculated.
  /// If not specified, defaults to all simulated variables in the ObsSpace.
  oops::OptionalParameter<std::vector<Variable>> filterVariables{"filter variables", this};

  /// \brief Types of ensemble statistics to calculate.
  ///
  /// Each statistic will be saved in ObsSpace variables from the group with matching name (e.g.
  /// ensemble means of H(x) vectors will be saved in variables from the `HofXEnsembleMean`group).
  oops::Parameter<std::vector<EnsembleStatistic>> statistics{"statistics", {}, this};
};

/// \brief Compute ensemble statistics and store them in ObsSpace variables.
///
/// See EnsembleStatisticsParameters for the documentation of the parameters controlling this
/// filter.
class EnsembleStatistics : public ObsProcessorBase,
                           private util::ObjectCounter<EnsembleStatistics> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef EnsembleStatisticsParameters Parameters_;

  static const std::string classname() {return "ufo::EnsembleStatistics";}

  EnsembleStatistics(ioda::ObsSpace &, const Parameters_ &,
                     std::shared_ptr<ioda::ObsDataVector<int> >,
                     std::shared_ptr<ioda::ObsDataVector<float> >);
  ~EnsembleStatistics() override;

 private:
  void print(std::ostream &) const override;
  void doFilter() override;

  oops::ObsVariables getFilterVariables() const;

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_ENSEMBLESTATISTICS_H_
