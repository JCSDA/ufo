/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_PROBABILITYGROSSERRORWHOLEREPORT_H_
#define UFO_FILTERS_PROBABILITYGROSSERRORWHOLEREPORT_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/ProbabilityOfGrossErrorParameters.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the ProbabilityGrossErrorWholeReport filter.
class ProbabilityGrossErrorWholeReportParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ProbabilityGrossErrorWholeReportParameters, FilterParametersBase)

 public:
  /// Parameters related to PGE calculations. The value of \p PGECrit
  /// is obtained from here.
  ProbabilityOfGrossErrorParameters PGEParameters{this};
};

/// \brief This filter calculates the probability that an entire report is affected by gross error.
/// \details Synoptic stations typically provide reports at regular intervals. A report is a
///   combination of variables observed by different sensors at a single location. Reports may
///   include some, but not necessarily all, of pressure, temperature, dew point and wind speed
///   and direction.
///
///   The probability that the whole report is affected by gross error is calculated
///   through the Bayesian combination of the probability of gross error of individual
///   observations. This is based on the logic that if multiple observations within a report
///   appear dubious based on a Bayesian Background check, it is likely that the whole report
///   is affected by, for example, location error. This filter should be called after the
///   Bayesian Background Check.
///
///   Variables which are to have their probability of gross error updated should be specified
///   using the "filter variables" YAML option. All variables included in "filter variables" will
///   be used to calculate the probability that the whole report is affected by gross error unless
///   the option \c not_used_in_whole_report is set to true for that variable.
///
///   Variables can be either scalar or vector (with two Cartesian components, such as the eastward
///   and northward wind components). In the latter case the two components need to specified one
///   after the other in the "filter variables" list, with the second component having the
///   \c second_component_of_two option set to true. For each variable, the option
///   \c Probability_Density_Bad is used to set the prior probability density of that variable being
///   "bad". The filter can also apply a specific prior probability density of bad observations for
///   the following Met Office SubTypes:
///   * Bogus
///   * Synop (SynopManual, SynopAuto, MetarManual, MetarAuto, SynopMob, SynopBufr, WOW)
///
/// Example:
///
/// \code{.yaml}
///- filter: Bayesian Whole Report
///  filter variables:
///  - name: pressure_at_model_surface
///    options:
///      Probability_Density_Bad: 0.1
///      Bogus_Probability_Density_Bad: 0.1
///  - name: air_temperature_at_2m
///    options:
///      Probability_Density_Bad: 0.1
///  - name: eastward_wind
///    options:
///      Probability_Density_Bad: 0.1
///      Synop_Probability_Density_Bad: 0.1
///      Bogus_Probability_Density_Bad: 0.1
///  - name: northward_wind
///    options:
///      not_used_in_whole_report: true
///      second_component_of_two: true
///  - name: relative_humidity_at_2m
///    options:
///      not_used_in_whole_report: true
///      Probability_Density_Bad: 0.1
///  PGE threshold: 0.15
/// \endcode
///

class ProbabilityGrossErrorWholeReport : public FilterBase,
                        private util::ObjectCounter<ProbabilityGrossErrorWholeReport> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ProbabilityGrossErrorWholeReportParameters Parameters_;

  static const std::string classname() {return "ufo::ProbabilityGrossErrorWholeReport";}

  ProbabilityGrossErrorWholeReport(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int>>,
                  std::shared_ptr<ioda::ObsDataVector<float>>);
  ~ProbabilityGrossErrorWholeReport();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::bayesianQC;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PROBABILITYGROSSERRORWHOLEREPORT_H_
