/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ProcessAMVQI.h"

#include <algorithm>
#include <string>
#include <vector>

#include "ioda/distribution/Accumulator.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

#include "ufo/filters/processWhere.h"

namespace ufo {

// -----------------------------------------------------------------------------

ProcessAMVQI::ProcessAMVQI(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                           std::shared_ptr<ioda::ObsDataVector<int>> flags,
                           std::shared_ptr<ioda::ObsDataVector<float>> obserr)
  : ObsProcessorBase(obsdb, false /*deferToPost?*/, flags, obserr),
    parameters_(parameters)
{
  oops::Log::trace() << "ProcessAMVQI contructor starting" << std::endl;
}

// -----------------------------------------------------------------------------

ProcessAMVQI::~ProcessAMVQI() {
  oops::Log::trace() << "ProcessAMVQI destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/*! \brief A filter to convert new BUFR formatted data into variables with names
 *  corresponding to the generating application.
 *
 *  Example:
 *  \code{.unparsed}
 *  obs filter:
 *  - filter: Process AMV QI
 *    number of generating apps: 4
 *  \endcode
 *
 *  \author A.Martins (Met Office)
 *
 *  \date 02/08/2021: Created
 */
void ProcessAMVQI::doFilter() {
  oops::Log::trace() << "ProcessAMVQI doFilter" << std::endl;

  const float missing = util::missingValue<float>();
  const int int_missing = util::missingValue<int>();
  const size_t nlocs = obsdb_.nlocs();
  const size_t number_of_apps = parameters_.number_of_apps.value();

  // vectors to store BUFR data
  std::vector<float> percent_confidence(nlocs);
  std::vector<int> wind_generating_application(nlocs);

  // Create vector of strings for windPercentConfidence<number> to new variable
  // Wind generating application number = QI type
  // Table 1:
  // 1 = Full weighted mixture of individual quality tests
  // 2 = Weighted mixture of individual tests, but excluding forecast comparison
  // 3 = Recursive filter function
  // 4 = Common quality index (QI) without forecast
  // 5 = QI without forecast
  // 6 = QI with forecast
  // 7 = Estimated Error (EE) in m/s converted to a percent confidence
  std::vector<std::string> new_variables = {
    "qiWeightedMixtureFull",
    "qiWeightedMixtureWithoutForecast",
    "qiRecursiveFilterFunction",
    "qiCommon",
    "qiWithoutForecast",
    "qiWithForecast",
    "qiEstimatedError" };
  const size_t total_variables = new_variables.size();

  // create variable to write to obsdb
  std::vector<std::vector<float>> new_percent_confidence(total_variables,
                                                         std::vector<float>(nlocs, missing));

  // diagnostic variables to be summed over all processors
  std::unique_ptr<ioda::Accumulator<size_t>> count_bad_app_accumulator =
    obsdb_.distribution()->createAccumulator<size_t>();

  // Get BUFR data
  for (size_t iapp = 0; iapp < number_of_apps; ++iapp) {
    // names of variables
    std::string percent_confidence_name = "windPercentConfidence";
    std::string wind_generating_application_name = "windGeneratingApplication";

    obsdb_.get_db("MetaData",
                  percent_confidence_name.append(std::to_string(iapp + 1)),
                  percent_confidence);

    obsdb_.get_db("MetaData",
                  wind_generating_application_name.append(std::to_string(iapp + 1)),
                  wind_generating_application);

    for (size_t idata = 0; idata < nlocs; ++idata) {
      if (wind_generating_application[idata] != int_missing) {
        // check index is in range before filling new percent confidence values
        if (wind_generating_application[idata] >= 1 &&
            wind_generating_application[idata] <= total_variables) {
          new_percent_confidence[wind_generating_application[idata] - 1][idata] =
              percent_confidence[idata];
        } else {
          count_bad_app_accumulator->addTerm(idata, 1);
        }
      }
    }
  }

  // sum number of generating application values that are out of range
  const std::size_t count_bad_app = count_bad_app_accumulator->computeResult();
  if (count_bad_app > 0) {
    oops::Log::warning() << "Process AMV QI: " << count_bad_app
      << " instances of generating application out of range" << std::endl;
  }

  // Need to check database for the named variables and add to them if they exist.
  for (size_t inum = 0; inum < total_variables; ++inum) {
    if (std::any_of(new_percent_confidence[inum].begin(),
                    new_percent_confidence[inum].end(),
                    [missing] (float elem) {return elem != missing;})) {
      // Check if new_variable already exists in obsdb
      std::vector<float> existing_new_percent_confidence(nlocs, missing);
      if (obsdb_.has("MetaData",
                     new_variables[inum])) {
        obsdb_.get_db("MetaData",
                      new_variables[inum],
                      existing_new_percent_confidence);
      } else {
        oops::Log::trace() << "ProcessAMIQI: New variable: " << new_variables[inum] <<
                              " not present in input data" << std::endl;
      }

      // Check if existing_new_percent_confidence has non-missing values
      // and update new_percent_confidence.

      for (size_t idata = 0; idata < nlocs; ++idata) {
        if (existing_new_percent_confidence[idata] != missing) {
          if (new_percent_confidence[inum][idata] == missing) {
            new_percent_confidence[inum][idata] =
                existing_new_percent_confidence[idata];
          } else {
            oops::Log::trace() << "[WARN] New variable (created from new BUFR data)" << "\n"
                                  "[WARN] already has value in old BUFR format" << "\n"
                                  "[WARN] new BUFR: " <<
                                  new_percent_confidence[inum][idata] << "\n"
                                  "[WARN] old BUFR: " <<
                                  existing_new_percent_confidence[idata] << "\n";
          }
        }
      }
      // write to new or existing QI vectors
      obsdb_.put_db("MetaData",
                    new_variables[inum],
                    new_percent_confidence[inum]);
    }
  }
}

// -----------------------------------------------------------------------------

void ProcessAMVQI::print(std::ostream & os) const {
  os << "ProcessAMVQI filter" << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
