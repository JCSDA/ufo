/*
 * (C) Crown copyright 2022, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_OCEANVERTICALSTABILITYCHECK_H_
#define UFO_FILTERS_OCEANVERTICALSTABILITYCHECK_H_

#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/OceanConversions/OceanConversions.interface.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the variables of the Ocean Vertical Stability Check filter.
class OceanVerticalStabilityCheckVarParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OceanVerticalStabilityCheckVarParameters, Parameters)

 public:
  /// Salinity variable, g/kg (from which density is calculated)
  oops::RequiredParameter<Variable> salinity{"salinity", this};

  /// Temperature variable, deg.C (from which density is calculated)
  oops::RequiredParameter<Variable> temperature{"temperature", this};

  /// Pressure variable, dbar (from which density is calculated)
  oops::RequiredParameter<Variable> pressure{"pressure", this};
};

/// Parameters controlling the operation of the Ocean Vertical Stability Check filter.
class OceanVerticalStabilityCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(OceanVerticalStabilityCheckParameters, FilterParametersBase)

 public:
  /// 'variables' parameter class
  oops::RequiredParameter<OceanVerticalStabilityCheckVarParameters> varParams{"variables", this};

  /// If false, do not count spikes. Default: count both spikes and steps.
  oops::Parameter<bool> yesSpikes{"count spikes", true, this};

  /// If false, do not count steps. Default: count both spikes and steps.
  oops::Parameter<bool> yesSteps{"count steps", true, this};

  /// Nominal tolerance for density inversions (kg/m^3).
  oops::Parameter<float> nominal{"nominal tolerance", -0.05, this};

  /// For conditions checking for small spike.
  oops::Parameter<float> threshold{"threshold", 0.25, this,
                                  {oops::exclusiveMinConstraint(0.0f)}};
};

/// \brief OceanVerticalStabilityCheck: Flag observations that are likely erroneous because they
///   constitute a density inversion (density should increase with depth).
/// \details Using the within each record, loop through the observations, checking conditions
///  on level-to-level potential density difference, for whether the obs is:
///  - density spike (drho < tolerance AND |drho[i]+dhro[i+1]| < threshold*|drho[i]-dhro[i+1]|),
///  - density inversion step (drho < tolerance AND NOT spike; flag in both cases)
///
/// Requires the following be specified in .yaml, under...
///
/// obs filters:
/// - filter: Ocean Vertical Stability Check
///   * variables.salinity    # salinity (g/kg),
///   * variables.temperature # temperature (deg.C),
///   * variables.pressure    # pressure (dbar) variables, from which density is calculated.
///
/// May also specify the following optional parameters:
///   * count spikes      # set to false to ignore spikes (default true).
///   * count steps       # set to false to ignore steps (default true).
///   * nominal tolerance # tolerance to density spikes/steps (default -0.05 kg/m^3).
///   * threshold         # use to check spike condition (default 0.25).
///
class OceanVerticalStabilityCheck : public FilterBase,
                        private util::ObjectCounter<OceanVerticalStabilityCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef OceanVerticalStabilityCheckParameters Parameters_;

  static const std::string classname() {return "ufo::OceanVerticalStabilityCheck";}

  OceanVerticalStabilityCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                              std::shared_ptr<ioda::ObsDataVector<int> > flags,
                              std::shared_ptr<ioda::ObsDataVector<float> > obserr);
  ~OceanVerticalStabilityCheck();

 private:
  void print(std::ostream &) const override;

  /// \brief Apply Ocean Vertical Stability filter. Return flagged=true for rejected obs.
  void applyFilter(const std::vector<bool> & apply, const Variables & filtervars,
                   std::vector<std::vector<bool>> & flagged) const override;

  /// \brief Return track flag for observations rejected by Ocean Vertical Stability check.
  int qcFlag() const override {return QCflags::track;}

  /// \brief Given salinity (g/kg), temperature (deg.C) and pressure (dbar),
  ///  set density difference, drho, for the given record (group).
  std::vector<float> calculateDensityDiff(const std::vector<float> &salinity,
                                          const std::vector<float> &temperature,
                                          const std::vector<float> &pressure,
                                          const std::vector<size_t> &obs_indices) const;

  /// \brief Go through the obs in the record (group), finding which are spikes and steps,
  ///  according to given conditions, flag as appropriate.
  void identifyDensityInversions(std::vector<bool> &isThinned,
                                   std::vector<bool> &spikeFlag,
                                   std::vector<bool> &stepFlag,
                                   const std::vector<float> &densityDiff,
                                   const std::vector<size_t> &obs_indices) const;

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OCEANVERTICALSTABILITYCHECK_H_
