/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <list>
#include <string>
#include <unordered_map>
#include <utility>             // pair
#include <vector>

#include <boost/variant.hpp>
#include "Eigen/Core"
#include "unsupported/Eigen/CXX11/Tensor"

#include "eckit/exception/Exceptions.h"
#include "eckit/utils/StringTools.h"

#include "ioda/Engines/HH.h"
#include "ioda/Group.h"
#include "ioda/Misc/SFuncs.h"   // for convertV1PathToV2Path

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/utils/dataextractor/DataExtractorInput.h"
#include "ufo/utils/dataextractor/DataExtractorNetCDFBackend.h"

namespace ufo
{

namespace {

/// \brief Helper function for determining the dimension mapping names for a variable
std::vector<std::string> fetchDimNameMapping(
    const ioda::Variable &variable,
    const std::string &varName,
    const std::list<std::pair<std::string, ioda::Variable>> &coordinates) {
  std::vector<std::string> dimnames;
  if (variable.isDimensionScale()) {
    // Variable corresponds to the dimension name so maps to itself (it's a dimension coordinate).
    dimnames = {varName};
  } else {
    // Lookup the dimension mapping for this variable.
    // Loop over the dimension scales associated with a particular dimension.
    auto dsm = variable.getDimensionScaleMappings(coordinates);
    for (const auto &item : dsm) {
      // There shouldn't be more than 1 dimension scale per dimension.
      if (item.size() != 1)
        throw eckit::Exception("Variable '" + varName + "' has a multi-dimensional coordinate. "
                               "This is not supported.", Here());
      dimnames.emplace_back(item[0].first);
    }
  }
  return dimnames;
}

/// \brief Add variable `var` of type `T` to `coordsVals` under key `key`.
template<typename T>
void updateVariable(const std::string &key, ioda::Variable var,
                    DataExtractorInput::Coordinates &coordsVals) {
  // Read the variable from the input file
  std::vector<T> values = var.readAsVector<T>();

  // Replace source fill values with corresponding missing marks
  if (var.hasFillValue()) {
    ioda::detail::FillValueData_t sourceFvData = var.getFillValue();
    const T sourceFillValue = ioda::detail::getFillValue<T>(sourceFvData);
    // TODO(someone): This won't call the correct overload for datetimes (which are treated
    // as strings by ioda-engines).
    const T oopsFillValue = util::missingValue(oopsFillValue);
    if (oopsFillValue != sourceFillValue)
      std::replace(values.begin(), values.end(), sourceFillValue, oopsFillValue);
  }

  // Store the variable in coordsVals_
  coordsVals.emplace(key, std::move(values));
}

/// \brief Add variable `var` to `coordsVals` under key `key`.
void update(const std::string &key, ioda::Variable var,
            DataExtractorInput::Coordinates &coordsVals) {
  if (var.isA<int>()) {
    updateVariable<int>(key, var, coordsVals);
  } else if (var.isA<float>()) {
    updateVariable<float>(key, var, coordsVals);
  } else if (var.isA<std::string>()) {
    updateVariable<std::string>(key, var, coordsVals);
  } else {
    throw eckit::Exception("Data type not yet supported.", Here());
  }
}

/// \brief Helper function to determine dimension index from a dimension name
std::vector<int> getDimMapping(const std::vector<std::string> &dimnames,
                               const std::vector<std::string> &interpolatedArrayDimnames) {
  std::vector<int> dimIndex;
  // Loop over dimension names provided
  for (const std::string & dimname : dimnames) {
    int ind = -1;
    for (const std::string &obsDimname : interpolatedArrayDimnames) {
      ind += 1;
      if (dimname == obsDimname) {
        dimIndex.emplace_back(ind);
      }
    }
  }
  if (dimIndex.size() != dimnames.size()) {
    oops::Log::debug() << "One or more of the dimension names provided do not map to the "
                          "array to be interpolated: ";
    for (std::string dimname : dimnames) oops::Log::debug() << dimname << " ";
      oops::Log::debug() << std::endl;
    throw eckit::Exception(
      "Dimension mappings not associated with the array to be interpolated are not "
      "currently supported.", Here());
  }
  return dimIndex;
}

/// \brief Return the name of the unique variable whose name ends with `@` followed by
/// `payloadGroup` or begins with `payloadGroup` followed by '/'.
///
/// Throw an exception if there's no such variable or there's more than one.
const std::string &findPayloadVariable(const std::vector<std::string> &varNames,
                                       const std::string &payloadGroup) {
  const std::string prefix = payloadGroup + '/';
  const std::string suffix = '@' + payloadGroup;
  auto isInPayloadGroup = [&prefix, &suffix](const std::string &name) {
    return eckit::StringTools::beginsWith(name, prefix) ||
           eckit::StringTools::endsWith(name, suffix);
  };
  auto payloadVarIt = std::find_if(varNames.begin(), varNames.end(), isInPayloadGroup);
  if (payloadVarIt == varNames.end())
    throw eckit::UserError("No payload column found: no column name begins with '" + prefix +
                           "' or ends with '" + suffix + "'",
                           Here());
  if (std::any_of(payloadVarIt + 1, varNames.end(), isInPayloadGroup))
    throw eckit::UserError("Multiple payload candidates found: "
                           "more than one column name begins with '" + prefix +
                           "' or ends with '" + suffix + "'", Here());
  return *payloadVarIt;
}

}  // namespace

DataExtractorNetCDFBackend::DataExtractorNetCDFBackend(const std::string &filepath)
  : filepath_(filepath)
{}

DataExtractorInput DataExtractorNetCDFBackend::loadData(
    const std::string &interpolatedArrayGroup) const {
  DataExtractorInput result;

  // Open the input file
  const ioda::Group group = ioda::Engines::HH::openFile(
        filepath_, ioda::Engines::BackendOpenModes::Read_Only);

  // All our coords from the file.
  std::vector<std::string> vars = group.listObjects(ioda::ObjectType::Variable, true)
                                       .at(ioda::ObjectType::Variable);

  // Find the array to be interpolated
  std::string interpolatedArrayName = findPayloadVariable(vars, interpolatedArrayGroup);

  // Loop over ALL coords and fetch the corresponding Variable object.
  std::unordered_map<std::string, ioda::Variable> coords;
  for (const auto &coord_name : vars) {
    ioda::Variable cvar = group.vars[coord_name];
    coords.emplace(coord_name, cvar);
  }
  // Turn this into a list of pairs - ioda engines works with this.
  std::list<std::pair<std::string, ioda::Variable>> lcoords(coords.begin(), coords.end());

  // Process the metadata of the array to be interpolated
  // --------------------------------
  const ioda::Variable &interpolatedArrayCoord = coords[interpolatedArrayName];
  // The variable to be interpolated shouldn't be a dimension mapping
  ASSERT(interpolatedArrayCoord.isDimensionScale() == false);
  // Fetch the dimension mapping names for this array
  std::vector<std::string> interpolatedArrayDimnames = fetchDimNameMapping(
        interpolatedArrayCoord, interpolatedArrayName, lcoords);
  const ioda::Dimensions dimdim = interpolatedArrayCoord.getDimensions();

  // Does the array to be interpolated represent a full covariance matrix (or a stack of them)?
  // If so, we'll need to extract diagonals.
  // --------------------------------
  std::string isCovariant {"false"};
  if (interpolatedArrayCoord.atts.exists("full")) {
    interpolatedArrayCoord.atts["full"].read(isCovariant);
  }

  // Remove the array to be interpolated from our list of coords
  vars.erase(std::remove(vars.begin(), vars.end(), interpolatedArrayName), vars.end());

  // --------------------------------
  // Fetch dimension mappings for our NetCDF variables. Store them initially in more generic
  // containers than those used by result.coord2DimMapping and result.dim2CoordMapping; later
  // we'll post-process them and make sure they can be stored in the corresponding member variables
  // of 'result'.

  // Maps coordinate names to dimensions (0 or 1) of the payload array
  std::unordered_map<std::string, std::vector<int>> coord2DimMapping;
  // Maps dimensions of the payload array (0 or 1) to coordinate names
  std::unordered_map<int, std::vector<std::string>> dim2CoordMapping;
  // Coordinate values
  DataExtractorInput::Coordinates coordsVals;

  for (const std::string& varName : vars) {
    const ioda::Variable &variable = coords[varName];
    std::vector<std::string> dimnames = fetchDimNameMapping(variable, varName, lcoords);

    // Load our variables data and keep hold of it along with relevant information.
    update(varName, variable, coordsVals);

    // Create a mapping between name and array dimension and vice versa.
    std::vector<int> ddim = getDimMapping(dimnames, interpolatedArrayDimnames);
    for (int dim : ddim) {
      dim2CoordMapping[dim].emplace_back(varName);
    }
    coord2DimMapping[varName] = ddim;
  }

  // Load the array to be interpolated
  // --------------------------------
  // NOTE:
  // - We might want to eventually put together a more generalised approach (a collapse
  //   method taking the dimension as argument)
  // - We might also want to not assume that the array is ordered? - by using
  //   coordinate info. though I'm not sure why it wouldn't be...
  if (isCovariant == "true") {
    // This is a full matrix or a stack of full matrices - pull out the diagonals
    // --------------------------------
    int dimCollapse = 0;  // Dimension to collapse
    ASSERT(dimdim.dimsCur[0] == dimdim.dimsCur[1]);  // Sanity check the matrix is square.

    if (dimdim.dimensionality == 3) {
      // ioda::Variable requires us to define the size of the tensor beforehand
      // We use row major ordering so that it resemblence handling the NetCDF data structure.
      Eigen::Tensor<float, 3, Eigen::RowMajor> data(dimdim.dimsCur[0], dimdim.dimsCur[1],
                                                    dimdim.dimsCur[2]);
      interpolatedArrayCoord.readWithEigenTensor(data);
      result.payloadArray.resize(dimdim.dimsCur[1], dimdim.dimsCur[2]);

      for (int i = 0; i < data.dimension(0); i++) {
        for (int k = 0; k < data.dimension(2); k++) {
          result.payloadArray(i, k) = data(i, i, k);
        }
      }
    } else if (dimdim.dimensionality == 2) {
      Eigen::ArrayXXf data;
      interpolatedArrayCoord.readWithEigenRegular(data);
      result.payloadArray.resize(dimdim.dimsCur[1], 1);

      for (int i = 0; i < data.rows(); i++) {
        result.payloadArray(i, 0) = data(i, i);
      }
    } else {
      throw eckit::Exception("Expecting 3D or 2D array for error covariance.", Here());
    }

    // Update our dimension mapping after collapsing a dimension.
    // --------------------------------
    // Remove the outermost dimension mapping of our data array.
    interpolatedArrayDimnames.erase(interpolatedArrayDimnames.begin());
    // Now update all affected coordinates (those who vary with dim0).
    for (const std::string &coord : dim2CoordMapping.at(dimCollapse))
      coord2DimMapping.erase(coord);
    dim2CoordMapping.clear();
    for (auto &coord : coord2DimMapping) {
      // Shift remaining dimension mappings since the collapse
      for (int &dim : coord.second) {
        dim--;
        dim2CoordMapping[dim].emplace_back(coord.first);  // update our reverse lookup
      }
    }
  } else if ((dimdim.dimensionality == 2) || (dimdim.dimensionality == 1)) {
    // This is just the diagonals - read directly into our container
    // --------------------------------
    interpolatedArrayCoord.readWithEigenRegular(result.payloadArray);
  } else {
    throw eckit::Exception("The array to be interpolated has an unsupported number of dimensions.",
                           Here());
  }

  // Store the coordinate values and dimension mappings in 'result', converting any ioda-v1-style
  // variable names (var@Group) to ioda-v2 style (Group/var)

  for (auto &coordAndVals : coordsVals)
    result.coordsVals[ioda::convertV1PathToV2Path(coordAndVals.first)] = coordAndVals.second;

  for (const auto &coordAndDims : coord2DimMapping) {
    if (coordAndDims.second.size() != 1)
      throw eckit::Exception("Coordinate '" + coordAndDims.first +
                             "'is used to index multiple dimensions of the payload array", Here());
    result.coord2DimMapping[ioda::convertV1PathToV2Path(coordAndDims.first)] =
        coordAndDims.second.front();
  }

  if (!dim2CoordMapping.empty()) {
    const int maxDimIndex = std::max_element(dim2CoordMapping.begin(),
                                             dim2CoordMapping.end())->first;
    result.dim2CoordMapping.resize(maxDimIndex + 1);
    for (const auto &dimAndCoords : dim2CoordMapping) {
      result.dim2CoordMapping[dimAndCoords.first] = dimAndCoords.second;
      for (std::string &varName : result.dim2CoordMapping[dimAndCoords.first])
        varName = ioda::convertV1PathToV2Path(varName);
    }
  }

  return result;
}

}  // namespace ufo
