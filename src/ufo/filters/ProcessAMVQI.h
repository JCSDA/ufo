/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_PROCESSAMVQI_H_
#define UFO_FILTERS_PROCESSAMVQI_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/generic/ObsFilterParametersBase.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/ObsProcessorBase.h"

namespace ufo {

/// \brief Parameters controlling the operation of the ProcessAMVQI filter.
class ProcessAMVQIParameters : public oops::ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ProcessAMVQIParameters, ObsFilterParametersBase)

 public:
  /// How many generating applications to search for.
  oops::RequiredParameter<size_t> number_of_apps {
    "number of generating apps",
    this
  };
};

/// \brief A filter to convert new BUFR formatted data into variables with names
/// corressponding to the wind generating application.
///
/// \details This filter will convert variables of "windPercentConfidence<number>" and
/// "windGeneratingApplication<number>" to variables named corresponding to the
/// generating application number (see Table 1).
///
/// Table 1:
/// Generating application number = QI type
/// 1 = Full weighted mixture of individual quality tests
/// 2 = Weighted mixture of individual tests, but excluding forecast comparison
/// 3 = Recursive filter function
/// 4 = Common quality index (QI) without forecast
/// 5 = QI without forecast
/// 6 = QI with forecast
/// 7 = Estimated Error (EE) in m/s converted to a percent confidence
///
/// See ProcessAMVQIParameters for the documentation of
/// the parameters controlling this filter.
class ProcessAMVQI : public ObsProcessorBase,
                     private util::ObjectCounter<ProcessAMVQI> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ProcessAMVQIParameters Parameters_;

  static const std::string classname() {return "ufo::ProcessAMVQI";}

  ProcessAMVQI(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
               std::shared_ptr<ioda::ObsDataVector<int>> flags,
               std::shared_ptr<ioda::ObsDataVector<float>> obserr);
  ~ProcessAMVQI() override;

 private:
  void print(std::ostream &) const override;
  void doFilter() override;

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PROCESSAMVQI_H_
