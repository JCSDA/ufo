/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <utility>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/PrintFilterData.h"

namespace ufo {

PrintFilterData::PrintFilterData(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : ObsProcessorBase(obsdb, parameters.deferToPost, std::move(flags), std::move(obserr)),
    parameters_(parameters)
{
  oops::Log::trace() << "PrintFilterData constructor" << std::endl;
  allvars_ += getAllWhereVariables(parameters.where);

  // Add all variables to be printed to `allvars_`.
  for (const VariablePrintParameters & variableParams : parameters_.variables.value()) {
    const Variable variable(variableParams.variable);
    if (data_.has(variable))
      allvars_ += variable;
  }
}

template <typename VariableType>
void PrintFilterData::getData(const Variable & variable, identifier<VariableType>) const {
  ioda::ObsDataVector<VariableType> variableData(obsdb_, variable.toOopsVariables());

  // Treatment depends on whether channels are present.
  if (variable.channels().size() == 0) {
    // If channels are not present use data_.get().
    data_.get(variable, variableData, parameters_.skipDerived.value());
    std::vector<VariableType> globalVariableData = variableData[0];
    if (!parameters_.printRank0)
      obsdb_.distribution()->allGatherv(globalVariableData);
    filterData_[getVariableNameWithChannel(variable, 0)] = std::move(globalVariableData);
  } else {
    // If channels are present use obsdb_.get_db().
    // A try-catch block is used here because obsdb_.has() does not take channels into account.
    for (size_t ich = 0; ich < variable.size(); ++ich) {
      const std::string variableWithChannel = variable.variable(ich);
      try {
        obsdb_.get_db(variable.group(), variableWithChannel, variableData[ich],
                      {}, parameters_.skipDerived);
      } catch (...) {
        oops::Log::info() << getVariableNameWithChannel(variable, ich)
                          << " not present in filter data" << std::endl;
        continue;
      }
      std::vector<VariableType> globalVariableData = variableData[ich];
      if (!parameters_.printRank0)
        obsdb_.distribution()->allGatherv(globalVariableData);
      filterData_[getVariableNameWithChannel(variable, ich)] = std::move(globalVariableData);
    }
  }
}

void PrintFilterData::getData(const Variable & variable, identifier<bool>) const {
  // todo(ctgh): Check whether this specialisation can be removed by adding
  // a bool implementation of allGatherv.
  ioda::ObsDataVector<bool> variableData(obsdb_, variable.toOopsVariables());

  // Treatment depends on whether channels are present.
  if (variable.channels().size() == 0) {
    // If channels are not present use data_.get().
    data_.get(variable, variableData, parameters_.skipDerived.value());
    // Note conversion to int from bool.
    std::vector<int> globalVariableData(variableData[0].begin(), variableData[0].end());
    if (!parameters_.printRank0)
      obsdb_.distribution()->allGatherv(globalVariableData);
    filterData_[getVariableNameWithChannel(variable, 0)] = std::move(globalVariableData);
  } else {
    // If channels are present use obsdb_.get_db().
    // A try-catch block is used here because obsdb_.has() does not take channels into account.
    for (size_t ich = 0; ich < variable.size(); ++ich) {
      const std::string variableWithChannel = variable.variable(ich);
      try {
        obsdb_.get_db(variable.group(), variableWithChannel, variableData[ich],
                      {}, parameters_.skipDerived);
      } catch (...) {
        oops::Log::info() << getVariableNameWithChannel(variable, ich)
                          << " not present in filter data" << std::endl;
        continue;
      }
      std::vector<int> globalVariableData(variableData[ich].begin(), variableData[ich].end());
      if (!parameters_.printRank0)
        obsdb_.distribution()->allGatherv(globalVariableData);
      filterData_[getVariableNameWithChannel(variable, ich)] = std::move(globalVariableData);
    }
  }
}

void PrintFilterData::getMultiLevelData(const Variable & variable,
                                        const std::vector<int> & levels) const {
  std::vector<float> variableData(obsdb_.nlocs());
  for (const int level : levels) {
    const std::string MultiLevelVariableName =
      this->getVariableNameAtLevel(variable.fullName(), level);
    // Ensure the level is not out of bounds.
    if (level < 0 || level >= data_.nlevs(variable)) {
      oops::Log::info() << MultiLevelVariableName << " not present in filter data" << std::endl;
      continue;
    }
    data_.get(variable, level, variableData);
    filterData_[MultiLevelVariableName] = std::move(variableData);
  }
}

template <typename VariableType>
void PrintFilterData::printVariable
(const std::string & varname, const int loc, const std::vector<int> & apply,
 identifier<VariableType>) const {
  if (!apply[loc]) {
    oops::Log::info() << std::right << std::setw(parameters_.columnWidth) << "masked by where";
  } else {
    const VariableType value = boost::get<std::vector<VariableType>>(filterData_[varname])[loc];
    if (value == util::missingValue<VariableType>())
      oops::Log::info() << std::right << std::setw(parameters_.columnWidth) << "missing";
    else
      oops::Log::info() << std::right << std::setw(parameters_.columnWidth) << value;
  }
}

void PrintFilterData::printVariable
(const std::string & varname, const int loc, const std::vector<int> & apply,
 identifier<bool>) const {
  if (!apply[loc]) {
    oops::Log::info() << std::right << std::setw(parameters_.columnWidth) << "masked by where";
  } else {
    const int value = boost::get<std::vector<int>>(filterData_[varname])[loc];
    oops::Log::info() << std::right << std::setw(parameters_.columnWidth) << value;
    // There is not currently a missing boolean value.
  }
}

std::string PrintFilterData::getVariableNameAtLevel(const std::string & varname,
                                                    const int level) const {
  std::stringstream ssvarname;
  ssvarname << varname << " (level " << level << ")";
  return ssvarname.str();
}

std::string PrintFilterData::getVariableNameWithChannel(const Variable & variable,
                                                        const int channel) const {
  std::stringstream ssvarname;
  ssvarname << variable.group() << "/" << variable.variable(channel);
  return ssvarname.str();
}

bool PrintFilterData::isMultiLevelData(Variable variable) const {
  const std::string groupName = variable.group();
  return groupName == "GeoVaLs" || groupName == "ObsDiag" || groupName == "ObsBiasTerm";
}

void PrintFilterData::getAllData() const {
  // Populate all requested vectors of filter data.
  for (const VariablePrintParameters & variableParams : parameters_.variables.value()) {
    const Variable variable = variableParams.variable;
    if (data_.has(variable)) {
      switch (data_.dtype(variable)) {
      case ioda::ObsDtype::Integer:
        this->getData<int>(variable);
        break;
      case ioda::ObsDtype::Float:
        if (isMultiLevelData(variable)) {
          const std::set<int> levelset = variableParams.levels;
          const std::vector<int> levels(levelset.begin(), levelset.end());
          this->getMultiLevelData(variable, levels);
        } else {
          this->getData<float>(variable);
        }
        break;
      case ioda::ObsDtype::String:
        this->getData<std::string>(variable);
        break;
      case ioda::ObsDtype::DateTime:
        this->getData<util::DateTime>(variable);
        break;
      case ioda::ObsDtype::Bool:
        this->getData<bool>(variable);
        break;
      default:
        throw eckit::BadParameter("Invalid variable type for printing", Here());
      }
    } else {
      oops::Log::info() << variable.fullName() << " not present in filter data" << std::endl;
    }
  }
}

int PrintFilterData::getMaxVariableNameLength() const {
  // Find maximum string length by cycling through all of the variables to print.
  int maxVariableNameLength = std::string("Location").length();
  for (const VariablePrintParameters & variableParams : parameters_.variables.value()) {
    const Variable variable = variableParams.variable;
    if (isMultiLevelData(variable)) {
      const std::set<int> levelset = variableParams.levels;
      const std::vector<int> levels(levelset.begin(), levelset.end());
      for (int level : levels) {
        const std::string MultiLevelVariableName =
          this->getVariableNameAtLevel(variable.fullName(), level);
        if (filterData_.find(MultiLevelVariableName) == filterData_.end())
          continue;
        maxVariableNameLength = MultiLevelVariableName.length() > maxVariableNameLength ?
          MultiLevelVariableName.length() :
          maxVariableNameLength;
      }
    } else {
      if (variable.channels().size() == 0) {
        maxVariableNameLength = variable.fullName().length() > maxVariableNameLength ?
          variable.fullName().length() :
          maxVariableNameLength;
      } else {
        for (size_t ich = 0; ich < variable.size(); ++ich) {
          const std::string varname = this->getVariableNameWithChannel(variable, ich);
          maxVariableNameLength = varname.length() > maxVariableNameLength ?
            varname.length() :
            maxVariableNameLength;
        }
      }
    }
  }
  return maxVariableNameLength;
}

void PrintFilterData::printAllData() const {
  // Set up values that govern the appearance of the output.
  const int maxVariableNameLength = this->getMaxVariableNameLength();
  const int nlocs = parameters_.printRank0 ? obsdb_.nlocs() : obsdb_.globalNumLocs();
  const int locmin = parameters_.locmin >= nlocs ? nlocs - 1 : parameters_.locmin;
  const int locmax = parameters_.locmax == 0 ? nlocs : parameters_.locmax;
  if (locmin > locmax)
    throw eckit::BadValue("Minimum location cannot be larger than maximum location", Here());
  const int columnWidth = parameters_.columnWidth;
  // Number of locations to print per table row
  // (adding 3 accounts for " | " appended to each column).
  const int nlocsPerRow =
    std::max((parameters_.maxTextWidth - maxVariableNameLength) / (columnWidth + 3), 1);

  // Select locations at which the filter will be applied.
  const std::vector<bool> apply = processWhere(parameters_.where, data_, parameters_.whereOperator);
  std::vector<int> globalApply(apply.begin(), apply.end());
  if (!parameters_.printRank0)
    obsdb_.distribution()->allGatherv(globalApply);

  // Loop over each group of locations and print the contents of each variable.
  for (int locgroup = locmin; locgroup < locmax; locgroup += nlocsPerRow) {
    // Print table header.
    oops::Log::info() << std::setw(maxVariableNameLength) << "Location" << " | ";
    for (int loc = locgroup; loc < locgroup + nlocsPerRow && loc < locmax; ++loc)
      oops::Log::info() << std::setw(columnWidth) << loc << " | ";
    oops::Log::info() << std::endl;
    // Print division bar below header.
    oops::Log::info() << std::string(maxVariableNameLength, '-') << "-+-";
    for (int loc = locgroup; loc < locgroup + nlocsPerRow && loc < locmax; ++loc)
      oops::Log::info() << std::string(columnWidth, '-') << "-+-";
    oops::Log::info() << std::endl;
    // Print each variable in turn.
    for (const VariablePrintParameters & variableParams : parameters_.variables.value()) {
      const Variable variable = variableParams.variable;
      if (!data_.has(variable))
        continue;
      if (isMultiLevelData(variable)) {
        const std::set<int> levels = variableParams.levels;
        for (int level : levels)  {
          const std::string MultiLevelVariableName =
            this->getVariableNameAtLevel(variable.fullName(), level);
          if (filterData_.find(MultiLevelVariableName) == filterData_.end())
            continue;
          oops::Log::info() << std::setw(maxVariableNameLength) << MultiLevelVariableName << " | ";
          for (int loc = locgroup; loc < locgroup + nlocsPerRow && loc < locmax; ++loc) {
            this->printVariable<float>(MultiLevelVariableName, loc, globalApply);
            oops::Log::info() << " | ";
          }
          oops::Log::info() << std::endl;
        }
      } else {
        for (size_t ich = 0; ich < variable.size(); ++ich) {
          const std::string varname = this->getVariableNameWithChannel(variable, ich);
          if (filterData_.find(varname) == filterData_.end())
            continue;
          oops::Log::info() << std::setw(maxVariableNameLength) << varname << " | ";
          for (int loc = locgroup; loc < locgroup + nlocsPerRow && loc < locmax; ++loc) {
            switch (data_.dtype(variable)) {
            case ioda::ObsDtype::Integer:
              this->printVariable<int>(varname, loc, globalApply);
              break;
            case ioda::ObsDtype::Float:
              this->printVariable<float>(varname, loc, globalApply);
              break;
            case ioda::ObsDtype::String:
              this->printVariable<std::string>(varname, loc, globalApply);
              break;
            case ioda::ObsDtype::DateTime:
              this->printVariable<util::DateTime>(varname, loc, globalApply);
              break;
            case ioda::ObsDtype::Bool:
              this->printVariable<bool>(varname, loc, globalApply);
              break;
            default:
              break;
            }
            oops::Log::info() << " | ";
          }
          oops::Log::info() << std::endl;
        }
      }
    }
    oops::Log::info() << std::endl;
  }
}

void PrintFilterData::doFilter() const {
  oops::Log::trace() << "PrintFilterData doFilter started" << std::endl;
  oops::Log::debug() << *this;

  // Print welcome message.
  oops::Log::info() << std::endl;
  oops::Log::info() << "############################" << std::endl;
  oops::Log::info() << "### Printing filter data ###" << std::endl;
  oops::Log::info() << "############################" << std::endl;
  oops::Log::info() << std::endl << std::endl;

  if (parameters_.message.value() != boost::none)
    oops::Log::info() << *parameters_.message.value() << std::endl << std::endl;

  if (parameters_.summary)
    oops::Log::info() << data_;

  this->getAllData();
  this->printAllData();

  oops::Log::trace() << "PrintFilterData doFilter finished" << std::endl;
}

void PrintFilterData::print(std::ostream & os) const {
  os << "PrintFilterData: config = " << parameters_ << std::endl;
}

}  // namespace ufo
