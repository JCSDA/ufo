/*
 * (C) Crown copyright 2022, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_SPIKEANDSTEPCHECK_H_
#define UFO_FILTERS_SPIKEANDSTEPCHECK_H_

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
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the tolerance options of the Spike and Step Check filter.
class SpikeAndStepToleranceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SpikeAndStepToleranceParameters, Parameters)

 public:
  /// Nominal tolerance for spikes of y against x, i.e. tolerance at x = 0.
  oops::RequiredParameter<float> nominal{"nominal value", this,
                                        {oops::exclusiveMinConstraint(0.0f)}};

  /// For conditions checking for small spike or consistent large gradient.
  oops::Parameter<float> threshold{"threshold", 0.5, this,
                                  {oops::exclusiveMinConstraint(0.0f)}};

  /// For conditions checking for small spike.
  oops::Parameter<float> smallThreshold{"small spike threshold", 0.25, this,
                                       {oops::exclusiveMinConstraint(0.0f)}};

  /// dy/dx tolerance, for checking for small spikes.
  oops::Parameter<float> gradientTolerance{"gradient", std::numeric_limits<float>::max(),
                                           this, {oops::exclusiveMinConstraint(0.0f)}};

  /// Resolution of dx for calculating dy/dx.
  oops::Parameter<float> dxResolution{"gradient x resolution",
                                      std::numeric_limits<float>::epsilon(), this,
                                      {oops::exclusiveMinConstraint(0.0f)}};

  /// Factors to multiply nominal tolerance to get tolerances in different x domains.
  ///  If not given, then nominal tolerance applies across whole x domain.
  oops::OptionalParameter<std::vector<float>> toleranceFactors{"factors", this};

  /// Boundary points of x domains with different tolerance values.
  ///  If not given, then nominal tolerance applies across whole x domain.
  oops::OptionalParameter<std::vector<float>> toleranceBoundaries{"x boundaries", this};
};

/// Parameters controlling the boundary layer options of the Spike and Step Check filter.
class SpikeAndStepBoundaryParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SpikeAndStepBoundaryParameters, Parameters)

 public:
  /// If x in this domain, the tolerance to a step in y is modified (see step tolerance range).
  oops::Parameter<std::vector<float>> boundaryRange{"x range", {0, 0}, this};

  /// If x within boundary layer range, then dy in this range does not count as a step.
  oops::Parameter<std::vector<float>> stepTolRange{"step tolerance range", {0, 0}, this};

  /// Ignore level if dx greater than this.
  oops::Parameter<std::vector<float>> maxDx{"maximum x interval",
                                           {std::numeric_limits<float>::max(),
                                            std::numeric_limits<float>::max()},
                                           this};
};

/// Parameters controlling the operation of the Spike and Step Check filter.
class SpikeAndStepCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(SpikeAndStepCheckParameters, FilterParametersBase)

 public:
  /// Independent variable, i.e. x of dy/dx
  oops::RequiredParameter<Variable> xVar{"independent", this};

  /// Dependent variable, i.e. y of dy/dx
  oops::RequiredParameter<Variable> yVar{"dependent", this};

  /// If false, do not count spikes. Default: count both spikes and steps.
  oops::Parameter<bool> yesSpikes{"count spikes", true, this};

  /// If false, do not count steps. Default: count both spikes and steps.
  oops::Parameter<bool> yesSteps{"count steps", true, this};

  /// Tolerance options (see class SpikeAndStepToleranceParameters).
  oops::RequiredParameter<SpikeAndStepToleranceParameters> toleranceOptions{"tolerance", this};

  /// Boundary layer options (see class SpikeAndStepBoundaryParameters).
  oops::OptionalParameter<SpikeAndStepBoundaryParameters> boundaryOptions{"boundary layer", this};
};

/// \brief SpikeAndStepCheck: Flag observations that are likely erroneous because they constitute
/// a spike or step within a record, e.g. profiles of ocean temperature vs. depth.
/// \details Using the series difference dy and dx computed from the dependent and independent
///  variables y and x within each record, loop through the observations, checking conditions
///  for whether the obs is:
///  - large spike (|dy| > tolerance(x) and opposite sign either side of the obs),
///  - small spike (|dy| > 0.5*tolerance(x) and dy/dx > gradient tolerance),
///  - step (|dy| > tolerance(x) and not spike; flag in these 3 cases),
///  - part of a large consistent gradient (i.e. not a one-off; accept),
///  - expected trend in boundary layer e.g. thermocline - see boundary layer options (accept),
///  - or too far from adjacent obs (dx too large; ignore).
///
/// Requires the following be specified in .yaml, under...
///
/// obs filters:
/// - filter: Spike and Step Check
///   * filter variables    # which variables should be flagged
///   * dependent           # dependent (y-)variable checked for spikes/steps
///   * independent         # independent (x-)variable against which y varies
///   * tolerance:          # tolerance options for spike/step conditions of y(x)
///     - nominal value     # tolerance value at x = 0
///
/// May also specify the following optional parameters:
///   * count spikes # set to false to ignore spikes (default true).
///   * count steps  # set to false to ignore steps (default true).
///   * tolerance:
///     - gradient  # if dy/dx greater, could be a spike.
///                 # (Don't check this condition if parameter not given.)
///     - gradient x resolution # precision to which dx can be known (otherwise dy/dx undefined).
///     - threshold  # use in combination with nominal tolerance for checking
///                  # spike/step conditions
///     - small spike threshold  # like 'threshold' but for small spikes
///     - factors      # tolerance (boundaries, factors) define variation of tolerance(x),
///     - x boundaries #  such that tolerance(x) = nominal * tolerance factor at x.
///   * boundary layer:
///     - x range  # when bounded by these two x values, you're in the boundary layer...
///     - step tolerance range # ...so relax tolerance factor for steps to within this range...
///     - maximum x interval   # ...and ignore an obs if its dx is greater than this when...
///                            # ...[within, outside] the range of the boundary layer.
///
class SpikeAndStepCheck : public FilterBase,
                        private util::ObjectCounter<SpikeAndStepCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef SpikeAndStepCheckParameters Parameters_;

  static const std::string classname() {return "ufo::SpikeAndStepCheck";}

  SpikeAndStepCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                  std::shared_ptr<ioda::ObsDataVector<int> > flags,
                  std::shared_ptr<ioda::ObsDataVector<float> > obserr);
  ~SpikeAndStepCheck();

  /// \brief Vectors of dependent and independent variables, their differences, and gradient
  ///  (obs locations per record).
  struct xyStruct {
    /// \brief \p yrec y-values in this record
    std::vector<float> yrec;
    /// \brief \p xrec x-values in this record
    std::vector<float> xrec;
    /// \brief \p ydiff consecutive differences of y-values
    std::vector<float> ydiff;
    /// \brief \p xdiff consecutive differences of x-values
    std::vector<float> xdiff;
    /// \brief \p dydx gradient (dy/dx) values
    std::vector<float> dydx;
  };

 private:
  void print(std::ostream &) const override;

  /// \brief Apply spike and step check filter. Return flagged=true for rejected obs.
  void applyFilter(const std::vector<bool> & apply, const Variables & filtervars,
                   std::vector<std::vector<bool>> & flagged) const override;

  /// \brief Return track flag for observations rejected by spike and step check.
  int qcFlag() const override {return QCflags::track;}

  /// \brief Given x (independent variable) and y (dependent variable),
  ///  set x, y, dx, dy and dy/dx for the given record (group).
  void set_xyrec(const std::vector<float> &x,
                 const std::vector<float> &y,
                 const std::vector<size_t> &obs_indices,
                 xyStruct &xy,
                 const Parameters_ &parameters_) const;

  /// \brief Given the nominal tolerance, tolerance factors and x-boundaries,
  ///  set the vector of tolerance values for every obs in the record (group).
  std::vector<float> set_tolerances(const std::vector<float> &xrec,
                                    const std::vector<size_t> &obs_indices,
                                    const Parameters_ &parameters_) const;

  /// \brief Go through the obs in the record (group), finding which are spikes and steps,
  ///  according to given conditions, flag as appropriate.
  void identifyThinnedObservations(std::vector<bool> &isThinned,
                                   std::vector<bool> &spikeFlag,
                                   std::vector<bool> &stepFlag,
                                   xyStruct &xy,
                                   const std::vector<float> &tolerances,
                                   const std::vector<size_t> &obs_indices,
                                   const Parameters_ &parameters_) const;

  /// Reimplemented to detect incompatible options.
  void validateParameters() const;

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_SPIKEANDSTEPCHECK_H_
