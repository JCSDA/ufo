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
    parameters_(parameters),
    os_(parameters.outputToTest ? oops::Log::test() : oops::Log::info())
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
  ioda::ObsDataVector<VariableType> variableData(obsdb_, variable.toOopsObsVariables());

  // Treatment depends on whether channels are present.
  if (variable.channels().size() == 0) {
    // If channels are not present use data_.get().
    data_.get(variable, variableData, parameters_.skipDerived.value());
    std::vector<VariableType> globalVariableData = variableData[0];
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
        os_ << getVariableNameWithChannel(variable, ich)
            << " not present in filter data" << std::endl;
        continue;
      }
      std::vector<VariableType> globalVariableData = variableData[ich];
      obsdb_.distribution()->allGatherv(globalVariableData);
      filterData_[getVariableNameWithChannel(variable, ich)] = std::move(globalVariableData);
    }
  }
}

void PrintFilterData::getData(const Variable & variable, identifier<bool>) const {
  ioda::ObsDataVector<bool> variableData(obsdb_, variable.toOopsObsVariables());

  // Treatment depends on whether channels are present.
  if (variable.channels().size() == 0) {
    // If channels are not present use data_.get().
    data_.get(variable, variableData, parameters_.skipDerived.value());
    // Note conversion to int from bool.
    std::vector<int> globalVariableData(variableData[0].begin(), variableData[0].end());
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
        os_ << getVariableNameWithChannel(variable, ich)
            << " not present in filter data" << std::endl;
        continue;
      }
      std::vector<int> globalVariableData(variableData[ich].begin(), variableData[ich].end());
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
      os_ << MultiLevelVariableName << " not present in filter data" << std::endl;
      continue;
    }
    data_.get(variable, level, variableData);
    std::vector<float> globalVariableData(variableData);
    obsdb_.distribution()->allGatherv(globalVariableData);
    filterData_[MultiLevelVariableName] = std::move(globalVariableData);
  }
}

template <typename VariableType>
void PrintFilterData::printVariable
(const std::string & varname, const int loc, identifier<VariableType>) const {
  const VariableType value = boost::get<std::vector<VariableType>>(filterData_[varname])[loc];
  if (value == util::missingValue<VariableType>()) {
    os_ << std::right << std::setw(parameters_.columnWidth) << "missing";
  } else {
    os_ << std::right
        << std::setw(parameters_.columnWidth)
        << (parameters_.scientificNotation ? std::scientific : std::fixed)
        << std::setprecision(parameters_.floatPrecision)
        << value;
  }
}

void PrintFilterData::printVariable
(const std::string & varname, const int loc, identifier<bool>) const {
  const int value = boost::get<std::vector<int>>(filterData_[varname])[loc];
  os_ << std::right << std::setw(parameters_.columnWidth) << value;
  // There is not currently a missing boolean value.
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
      os_ << variable.fullName() << " not present in filter data" << std::endl;
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
  const int nlocs = obsdb_.globalNumLocs();
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
  if (parameters_.printRank0 && obsdb_.comm().rank() != 0) {
    std::fill(globalApply.begin(), globalApply.end(), 0);
  }
  obsdb_.distribution()->allGatherv(globalApply);
  // Obtain global indices of each location in the ObsSpace.
  const std::vector<std::size_t> index = obsdb_.index();
  std::vector<int> globalIndex(index.begin(), index.end());
  obsdb_.distribution()->allGatherv(globalIndex);

  // Determine which locations to print.
  std::vector<int> locsToPrint;
  // Also record the indices of these locations in the global location vector.
  std::vector<int> indicesToPrint;

  for (size_t idx = 0; idx < globalIndex.size(); ++idx) {
    const int loc = std::distance(globalIndex.begin(),
                                  std::find(globalIndex.begin(), globalIndex.end(), idx));
    if (idx >= locmin && idx < locmax && globalApply[loc]) {
      indicesToPrint.push_back(idx);
      locsToPrint.push_back(loc);
    }
  }

  // Rows of locations (and corresponding indices in the global location vector)
  // to print in the output table.
  std::vector<std::vector<int>> rowOfLocsToPrint;
  std::vector<std::vector<int>> rowOfIndicesToPrint;
  int count = 0;
  for (int i = 0; i < locsToPrint.size(); ++i) {
    const int idx = indicesToPrint[i];
    const int loc = locsToPrint[i];
    if (count % nlocsPerRow == 0) {
      rowOfLocsToPrint.push_back({});
      rowOfIndicesToPrint.push_back({});
    }
    rowOfLocsToPrint.back().push_back(loc);
    rowOfIndicesToPrint.back().push_back(idx);
    count++;
  }

  // Print each row in turn.
  for (int i = 0; i < rowOfLocsToPrint.size(); ++i) {
    const auto locGroup = rowOfLocsToPrint[i];
    const auto idxGroup = rowOfIndicesToPrint[i];
    // Print table header.
    os_ << std::setw(maxVariableNameLength) << "Location" << " | ";
    for (int idx : idxGroup)
      os_ << std::setw(columnWidth) << idx << " | ";
    os_ << std::endl;
    // Print division bar below header.
    os_ << std::string(maxVariableNameLength, '-') << "-+-";
    for (int loc : locGroup)
      os_ << std::string(columnWidth, '-') << "-+-";
    os_ << std::endl;
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
          os_ << std::setw(maxVariableNameLength) << MultiLevelVariableName << " | ";
          for (int loc : locGroup) {
            this->printVariable<float>(MultiLevelVariableName, loc);
            os_ << " | ";
          }
          os_ << std::endl;
        }
      } else {
        for (size_t ich = 0; ich < variable.size(); ++ich) {
          const std::string varname = this->getVariableNameWithChannel(variable, ich);
          if (filterData_.find(varname) == filterData_.end())
            continue;
          os_ << std::setw(maxVariableNameLength) << varname << " | ";
          for (int loc : locGroup) {
            switch (data_.dtype(variable)) {
            case ioda::ObsDtype::Integer:
              this->printVariable<int>(varname, loc);
              break;
            case ioda::ObsDtype::Float:
              this->printVariable<float>(varname, loc);
              break;
            case ioda::ObsDtype::String:
              this->printVariable<std::string>(varname, loc);
              break;
            case ioda::ObsDtype::DateTime:
              this->printVariable<util::DateTime>(varname, loc);
              break;
            case ioda::ObsDtype::Bool:
              this->printVariable<bool>(varname, loc);
              break;
            default:
              break;
            }
            os_ << " | ";
          }
          os_ << std::endl;
        }
      }
    }
    os_ << std::endl;
  }
}

void PrintFilterData::doFilter() {
  oops::Log::trace() << "PrintFilterData doFilter started" << std::endl;
  oops::Log::debug() << *this;

  // Print welcome message.
  os_ << std::endl;
  os_ << "############################" << std::endl;
  os_ << "### Printing filter data ###" << std::endl;
  os_ << "############################" << std::endl;
  os_ << std::endl << std::endl;

  if (parameters_.message.value() != boost::none)
    os_ << *parameters_.message.value() << std::endl << std::endl;

  if (parameters_.summary)
    os_ << data_;

  this->getAllData();
  this->printAllData();

  oops::Log::trace() << "PrintFilterData doFilter finished" << std::endl;
}

void PrintFilterData::print(std::ostream & os) const {
  os << "PrintFilterData: config = " << parameters_ << std::endl;
}

}  // namespace ufo
