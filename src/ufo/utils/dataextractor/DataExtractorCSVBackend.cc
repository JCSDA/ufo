/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <unordered_map>
#include <utility>  // for move
#include <vector>

#include <boost/multi_array.hpp>
#include <boost/variant.hpp>

#include "eckit/exception/Exceptions.h"
#include "eckit/parser/CSVParser.h"
#include "eckit/utils/StringTools.h"

#include "ioda/Misc/StringFuncs.h"   // for convertV1PathToV2Path
#include "ioda/ObsSpace.h"  // for ObsDtype

#include "oops/util/Logger.h"

#include "ufo/utils/dataextractor/DataExtractorCSVBackend.h"
#include "ufo/utils/dataextractor/DataExtractorInput.h"

namespace ufo
{

namespace {

/// Number of header rows in CSV files.
const size_t numHeaderRows = 2;

/// Representation of missing values in CSV files (same as in NetCDF's CDL).
const char *missingValuePlaceholder = "_";

/// Visitor whose operator() takes a vector and appends to it the value passed to the constructor
/// (treating "_" as a placeholder for "missing value").
class AppendValueVisitor : public boost::static_visitor<void> {
 public:
  explicit AppendValueVisitor(const eckit::Value &value) :
    value_(value)
  {}

  void operator()(std::vector<int> &values) const {
    int value;
    if (value_.as<std::string>() == missingValuePlaceholder)
      value = util::missingValue<int>();
    else
      // NOLINTNEXTLINE(runtime/int): It's not our fault that eckit uses the 'long long' type...
      value = static_cast<int>(value_.as<long long>());
    values.push_back(value);
  }

  void operator()(std::vector<float> &values) const {
    float value;
    if (value_.as<std::string>() == missingValuePlaceholder)
      value = util::missingValue<float>();
    else
      value = static_cast<float>(value_.as<double>());
    values.push_back(value);
  }

  void operator()(std::vector<std::string> &values) const {
    std::string value = value_.as<std::string>();
    if (value == missingValuePlaceholder)
      value = util::missingValue<std::string>();
    values.push_back(value);
  }

 private:
  const eckit::Value &value_;
};

template <typename Source, typename Destination>
void convertVectorToColumnArray(const std::vector<Source> &source,
                                    boost::multi_array<Destination, 3> &destination) {
  const Source missingSource = util::missingValue<Source>();
  const Destination missingDestination = util::missingValue<Destination>();
  destination.resize(boost::extents[source.size()][1][1]);
  for (size_t i = 0; i < source.size(); ++i)
    if (source[i] != missingSource)
      destination[i][0][0] = static_cast<Destination>(source[i]);
    else
      destination[i][0][0] = missingDestination;
}

/// Visitor that converts an std::vector to a boost::multi_array with one column.
template <typename ExtractedValue>
class ConvertToBoostMultiArrayVisitor : public boost::static_visitor<void> {
 public:
  explicit ConvertToBoostMultiArrayVisitor(boost::multi_array<ExtractedValue, 3> &output) :
    output_(output)
  {}

  template <typename T,
            typename std::enable_if<std::is_convertible<T, ExtractedValue>::value, bool>::type
            = true>
  void operator()(const std::vector<T> &values) {
    output_.resize(boost::extents[values.size()][1][1]);
    for (size_t i = 0; i < values.size(); ++i)
      output_[i][0][0] = values[i];
  }

  template <typename T,
            typename std::enable_if<!std::is_convertible<T, ExtractedValue>::value, bool>::type
            = true>
  void operator()(const std::vector<T> &) {
    // Should never be called
    throw eckit::NotImplemented(Here());
  }

 private:
  boost::multi_array<ExtractedValue, 3> &output_;
};

/// \brief Find the index of the column whose name ends with `@` followed by `payloadGroup`
/// or begins with `payloadGroup` followed by `/`.
///
/// Throw an exception if there's no such column or there's more than one.
size_t findPayloadColumn(const std::vector<std::string> &columnNames,
                         const std::string &payloadGroup) {
  const std::string prefix = payloadGroup + '/';
  const std::string suffix = '@' + payloadGroup;
  auto isInPayloadGroup = [&prefix, &suffix](const std::string &name) {
    return eckit::StringTools::beginsWith(name, prefix) ||
           eckit::StringTools::endsWith(name, suffix);
  };
  auto payloadColumnIt = std::find_if(columnNames.begin(), columnNames.end(), isInPayloadGroup);
  if (payloadColumnIt == columnNames.end())
    throw eckit::UserError("No payload column found: no column name begins with '" + prefix +
                           "' or ends with '" + suffix + "'",
                           Here());
  if (std::any_of(payloadColumnIt + 1, columnNames.end(), isInPayloadGroup))
    throw eckit::UserError("Multiple payload candidates found: "
                           "more than one column name begins with '" + prefix +
                           "' or ends with '" + suffix + "'", Here());
  return payloadColumnIt - columnNames.begin();
}

template <typename T>
std::vector<T> createColumn(size_t numValues) {
  std::vector<T> values;
  values.reserve(numValues);
  return values;
}

/// \brief Throw an exception if contents of columns of type `type` can't be converted to values
/// of type `ExtractedValue`.
template <typename ExtractedValue>
void checkPayloadColumnType(const std::string &type);

template <>
void checkPayloadColumnType<float>(const std::string &type) {
  if (type != "float" && type != "int")
    throw eckit::UserError("The payload column must contain numeric data", Here());
}

template <>
void checkPayloadColumnType<int>(const std::string &type) {
  if (type != "float" && type != "int")
    throw eckit::UserError("The payload column must contain numeric data", Here());
}

template <>
void checkPayloadColumnType<std::string>(const std::string &type) {
  if (type != "string" && type != "datetime")
    throw eckit::UserError("The payload column must contain strings or datetimes", Here());
}

}  // namespace

template <typename ExtractedValue>
DataExtractorCSVBackend<ExtractedValue>::DataExtractorCSVBackend(const std::string &filepath)
  : filepath_(filepath)
{}

template <typename ExtractedValue>
DataExtractorInput<ExtractedValue> DataExtractorCSVBackend<ExtractedValue>::loadData(
    const std::string &interpolatedArrayGroup) const {
  DataExtractorInput<ExtractedValue> result;

  const eckit::Value contents = eckit::CSVParser::decodeFile(filepath_, false /* hasHeader? */);
  const size_t numRows = contents.size();
  // Ensure we have at least three lines:
  // * column names
  // * data types
  // * one row of values.
  if (numRows <= numHeaderRows)
    throw eckit::UserError("No data could be loaded from the file '" + filepath_ + "'", Here());
  const size_t numValues = numRows - numHeaderRows;

  // Read column names from the first line
  const eckit::Value nameHeader = contents[0];
  const size_t numColumns = nameHeader.size();
  std::vector<std::string> columnNames(numColumns);
  columnNames.reserve(numColumns);
  for (size_t column = 0; column < numColumns; ++column)
    columnNames[column] = nameHeader[column].as<std::string>();

  const size_t payloadColumnIndex = findPayloadColumn(columnNames, interpolatedArrayGroup);

  // Now that we won't need to include column names in any further error messages, convert
  // them to the ioda-v2 convention (Group/var rather than var@Group)
  for (std::string &columnName : columnNames)
    columnName = ioda::convertV1PathToV2Path(columnName);

  // Read data types from the second line
  const eckit::Value typeHeader = contents[1];
  if (typeHeader.size() != numColumns)
    throw eckit::UserError("The number of columns in line 2 differs from that in line 1", Here());

  // Allocate vectors for values to be loaded from subsequent lines
  std::vector<DataExtractorInputBase::Coordinate> columns(numColumns);
  for (size_t column = 0; column < numColumns; ++column) {
    const std::string type = typeHeader[column];
    if (column == payloadColumnIndex)
      checkPayloadColumnType<ExtractedValue>(type);
    if (type == "string" || type == "datetime") {
      columns[column] = createColumn<std::string>(numValues);
    } else if (type == "int" || type == "integer") {
      columns[column] = createColumn<int>(numValues);
    } else if (type == "float") {
      columns[column] = createColumn<float>(numValues);
    } else {
      throw eckit::UserError("Unsupported data type '" + type + "'", Here());
    }
  }

  // Load values from the rest of the CSV file
  for (size_t row = numHeaderRows; row < numRows; ++row) {
    const eckit::Value rowContents = contents[row];
    if (rowContents.size() == 1 && rowContents[0] == "")
      continue;  // empty line
    if (rowContents.size() != numColumns)
      throw eckit::UserError("The number of columns in line " + std::to_string(1 + row) +
                             " differs from that in line 1", Here());
    for (size_t column = 0; column < numColumns; ++column)
      boost::apply_visitor(AppendValueVisitor(rowContents[column]), columns[column]);
  }

  // Store the loaded data in the result object
  const int firstDim = 0;
  result.dim2CoordMapping.resize(1);
  for (size_t column = 0; column < numColumns; ++column) {
    if (column == payloadColumnIndex) {
      ConvertToBoostMultiArrayVisitor<ExtractedValue> visitor(result.payloadArray);
      boost::apply_visitor(visitor, columns[column]);
    } else {
      result.coordsVals[columnNames[column]] = std::move(columns[column]);
      result.coord2DimMapping[columnNames[column]] = std::vector<int> {firstDim};
      result.coordNDims[columnNames[column]] = 1;
      result.dim2CoordMapping[firstDim].push_back(columnNames[column]);
    }
  }

  if (result.payloadArray.shape()[0] == 0)
    throw eckit::UserError("No data could be loaded from the file '" + filepath_ + "'", Here());

  return result;
}

// Explicit instantiations
template class DataExtractorCSVBackend<float>;
template class DataExtractorCSVBackend<int>;
template class DataExtractorCSVBackend<std::string>;

}  // namespace ufo
