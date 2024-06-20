/*
 * (C) Crown Copyright 2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_PRINTFILTERDATA_H_
#define UFO_FILTERS_PRINTFILTERDATA_H_

#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "boost/variant.hpp"

#include "ioda/core/ParameterTraitsObsDtype.h"
#include "oops/generic/ObsFilterParametersBase.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/ObsProcessorBase.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters related to the variables to be printed.
class VariablePrintParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VariablePrintParameters, Parameters)

 public:
  /// Variable information (name and optional channels).
  /// This can either be a string containing the variable name, e.g.
  ///
  ///       - variable: ObsValue/airTemperature
  ///
  /// or list both the name and channels to print, e.g.
  ///
  ///       - variable:
  ///           name: ObsValue/brightnessTemperature
  ///           channels: 1-7, 15-28
  ///
  /// See ufo/utils/parameters/ParameterTraitsVariable.h for further information.
  oops::RequiredParameter<Variable>
  variable{"variable",
            "Variable information (name and optional channels)",
            this};

  /// Set of levels to be printed for multi-level data.
  oops::Parameter<std::set<int>>
  levels{"levels",
         "Set of levels to be printed for multi-level data",
         {},
         this};
};

/// Parameters controlling the operation of the PrintFilterData filter.
class PrintFilterDataParameters : public oops::ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(PrintFilterDataParameters, ObsFilterParametersBase)

 public:
  /// Message to print at the start of the output.
  oops::OptionalParameter<std::string> message{"message", this};

  /// Print summary of ObsFilterData.
  oops::Parameter<bool> summary{"summary", true, this};

  /// List of variables to print.
  oops::Parameter<std::vector<VariablePrintParameters>> variables{"variables", {}, this};

  /// The filter data will be printed for all locations whose global unique ObsSpace indices
  /// are greater than or equal to this value.
  /// This option can be used in combination with `maximum location` to control which
  /// locations are printed.
  oops::Parameter<int> locmin{"minimum location", 0, this, {oops::minConstraint(0)}};

  /// The filter data will be printed for all locations whose global unique ObsSpace indices
  /// are less than this value.
  /// This option can be used in combination with `minimum location` to control which
  /// locations are printed.
  /// If this is set to zero then it will be redefined as the number of locations - 1.
  /// If the maximum location is less than the minimum location then an exception will be thrown.
  oops::Parameter<int> locmax{"maximum location", 0, this, {oops::minConstraint(0)}};

  /// Only get and print data from rank 0. The precise consequence of this being set to `true`
  /// depends on the ObsSpace distribution that is used (round robin, inefficient, etc.).
  /// If this option is `false`, data from all ranks are gathered and printed on rank 0.
  oops::Parameter<bool> printRank0{"print only rank 0", false, this};

  /// The maximum width (in characters) of the output text.
  oops::Parameter<int> maxTextWidth{"maxmimum text width", 120, this, {oops::minConstraint(0)}};

  /// The width (in characters) of the columns in the output table.
  oops::Parameter<int> columnWidth{"column width", 20, this, {oops::minConstraint(0)}};

  /// The precision of floating-point numbers.
  oops::Parameter<int> floatPrecision{"float precision", 3, this, {oops::minConstraint(0)}};

  /// Whether or not to use scientific notation. If false, fixed format is used.
  oops::Parameter<bool> scientificNotation {"scientific notation", false, this};

  /// Conditions used to select locations at which the filter data should be printed.
  /// If not specified, printing will be performed at all required locations.
  oops::Parameter<std::vector<WhereParameters>> where{"where", {}, this};

  /// Operator used to combine the results of successive `where` options at the same location.
  /// The available operators are `and` and `or`.
  oops::Parameter<WhereOperator> whereOperator{"where operator", WhereOperator::AND, this};

  /// If set to true, variable assignment will be done after the obs operator has been invoked
  /// (even if the filter doesn't require any variables from the GeoVaLs or HofX groups).
  oops::Parameter<bool> deferToPost{"defer to post", false, this};

  /// If this option is true, retrieval of a particular group name from the filter data will only
  /// consider that group. If this option is false, the retrieval will first check for the same
  /// group prefixed with "Derived"; if such a group is present then the data from that will be
  /// retrieved. If the Derived group is not present, data from the original group will then be
  /// retrieved.
  oops::Parameter<bool> skipDerived{"skip derived", true, this};

  /// If set to true, output will appear in the oops `test` stream.
  oops::Parameter<bool> outputToTest{"output to test stream", false, this};
};

/// Type identifier used in calls to getData.
template <typename VariableType>
struct identifier{};

/// \brief Prints requested filter data whenever it is invoked.
/// \details This filter can be placed at any point in a configuration file.
/// When run, it will print out the values stored in a list of variables requested by the user.
/// Any variable stored in the ObsFilterData can be printed.
/// Model-level data such as GeoVaLs are printed out on individual (user-configured) levels.
class PrintFilterData : public ObsProcessorBase,
                        private util::ObjectCounter<PrintFilterData> {
 public:
  typedef PrintFilterDataParameters Parameters_;

  static const std::string classname() {return "ufo::PrintFilterData";}

  PrintFilterData(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                  std::shared_ptr<ioda::ObsDataVector<int> > flags,
                  std::shared_ptr<ioda::ObsDataVector<float> > obserr);

 private:  // variables
  Parameters_ parameters_;

  /// Collection of filter data.
  mutable std::unordered_map<std::string, boost::variant
                             <std::vector <int>,
                              std::vector <float>,
                              std::vector <std::string>,
                              std::vector <util::DateTime>,
                              std::vector <bool>>> filterData_;

  /// Output stream.
  std::ostream & os_;

 private:  // functions
  void print(std::ostream &) const override;
  void doFilter() override;

  /// Get the name of multi-level data at a particular level.
  std::string getVariableNameAtLevel(const std::string & varname, const int level) const;

  /// Get the name of a variable with an associated channel.
  std::string getVariableNameWithChannel(const Variable & variable, const int channel) const;

  /// Determine the maximum variable name length.
  int getMaxVariableNameLength() const;

  /// Determine whether a variable can exist on multiple levels at each location (e.g. GeoVaLs).
  bool isMultiLevelData(Variable variable) const;

  /// Get all requested data.
  void getAllData() const;

  /// Get the requested variable.
  /// Specialisation to bool is required because ioda::Distribution does not handle that type.
  template <typename VariableType>
  void getData(const Variable & variableName) const {
    this->getData(variableName, identifier<VariableType>());
  }
  template <typename VariableType>
  void getData(const Variable & variableName, identifier<VariableType>) const;
  void getData(const Variable & variableName, identifier<bool>) const;

  /// Get the requested multi-level data (in the GeoVaLs, ObsDiag, or ObsBiasTerm groups).
  void getMultiLevelData(const Variable & variableName,
                         const std::vector<int> & levels) const;

  /// Print all data.
  void printAllData() const;

  /// Print the requested variable.
  /// Specialisation to bool is required because ioda::Distribution does not handle that type.
  template <typename VariableType>
  void printVariable
  (const std::string & varname, const int loc) const {
    this->printVariable(varname, loc, identifier<VariableType>());
  }
  template <typename VariableType>
  void printVariable
  (const std::string & varname, const int loc, identifier<VariableType>) const;
  void printVariable
  (const std::string & varname, const int loc, identifier<bool>) const;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PRINTFILTERDATA_H_
