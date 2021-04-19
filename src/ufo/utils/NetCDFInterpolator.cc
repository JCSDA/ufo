/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/NetCDFInterpolator.h"

#include "ioda/Engines/Factory.h"
#include "ioda/Group.h"

#include "unsupported/Eigen/CXX11/Tensor"  // Eigen Tensors

/// \brief Boost visitor which allows us to sort a vector.
class SortUpdateVisitor : public boost::static_visitor<void> {
 public:
  explicit SortUpdateVisitor(ufo::RecursiveSplitter &splitter) : splitter(splitter) {}

  template <typename T>
  void operator()(const std::vector<T> &coord) {
    splitter.groupBy(coord);
  }

  void operator()(const std::vector<float> &coord) {
    splitter.sortGroupsBy(
      [&coord](int indexA, int indexB) {
        return coord[static_cast<size_t>(indexA)] < coord[static_cast<size_t>(indexB)];
      });
  }

  ufo::RecursiveSplitter &splitter;
};


/// \brief Boost visitor which allows us to sort a vector.
class SortVisitor : public boost::static_visitor<void> {
 public:
  explicit SortVisitor(const ufo::RecursiveSplitter &splitter) : splitter(splitter) {}

  template <typename T>
  void operator()(std::vector<T> &coord) {
    std::vector<T> newCoord;
    newCoord.reserve(coord.size());
    for (const auto &group : splitter.groups()) {
      for (const auto &index : group) {
        newCoord.push_back(coord[index]);
      }
    }
    // Replace the coordinate with the sorted one.
    coord = std::move(newCoord);
  }

  const ufo::RecursiveSplitter &splitter;
};


namespace ufo {

NetCDFInterpolator::NetCDFInterpolator(const std::string &filepath, const std::string &group) {
  // Read the data from the file
  load(filepath, group);
  // Start by constraining to the full range of our data
  resetExtract();
  // Initialise splitter for both dimensions
  splitter_.emplace_back(ufo::RecursiveSplitter(static_cast<size_t>(interpolatedArray2D_.rows())));
  splitter_.emplace_back(ufo::RecursiveSplitter(static_cast<size_t>(interpolatedArray2D_.cols())));
}


void NetCDFInterpolator::sort() {
  Eigen::ArrayXXf sortedArray = interpolatedArray2D_;
  extractMapIter_ = extractMap_.begin();

  for (auto &dimCoord : dim2CoordMapping_) {
    // Reorder coordinates
    for (auto &coord : dimCoord.second) {
      auto &coordVal = coordsVals_[coord];
      SortVisitor visitor(splitter_[static_cast<size_t>(dimCoord.first)]);
      boost::apply_visitor(visitor, coordVal);
    }
    // Reorder the array to be interpolated
    if (dimCoord.first == 0) {
      int ind = -1;
      for (const auto &group : splitter_[static_cast<size_t>(dimCoord.first)].groups()) {
        for (const auto &index : group) {
          ind++;
          oops::Log::debug() << "Sort index dim0; index-from: " << ind << " index-to: " <<
            index << std::endl;
          for (Eigen::Index j = 0; j < interpolatedArray2D_.cols(); j++) {
            sortedArray(ind, j) = interpolatedArray2D_(static_cast<Eigen::Index>(index), j);
          }
        }
      }
      // Replace the unsorted array with the sorted one.
      interpolatedArray2D_ = sortedArray;
    } else if (dimCoord.first == 1) {
      int ind = -1;
      for (const auto &group : splitter_[static_cast<size_t>(dimCoord.first)].groups()) {
        for (const auto &index : group) {
          ind++;
          oops::Log::debug() << "Sort index dim1; index-from: " << ind << " index-to: " <<
            index << std::endl;
          for (Eigen::Index i = 0; i < interpolatedArray2D_.rows(); i++) {
            sortedArray(i, ind) = interpolatedArray2D_(i, static_cast<Eigen::Index>(index));
          }
        }
      }
      // Replace the unsorted array with the sorted one.
      interpolatedArray2D_ = sortedArray;
    } else {
        throw eckit::Exception("Unable to reorder the array to be interpolated: "
                               "it has more than 2 dimensions.", Here());
    }
  }
}


void NetCDFInterpolator::scheduleSort(const std::string &varName, const InterpMethod &method) {
  auto &coordVal = coordsVals_[varName];
  int dimIndex;
  std::vector<int> dimIndices = coord2DimMapping_.at(varName);
  if (dimIndices.size() != 1) {
    throw eckit::Exception(
      "Expecting a coordinate mapping to a single dimension",
      Here());
  }
  dimIndex = dimIndices[0];

  SortUpdateVisitor visitor(splitter_[static_cast<size_t>(dimIndex)]);
  boost::apply_visitor(visitor, coordVal);

  // Update our map between coordinate (variable) and interpolation/extract method
  extractMap_.emplace_back(std::pair<std::string, InterpMethod> {varName, method});
}


void NetCDFInterpolator::resetExtract() {
  constrainedRanges_[0].begin = 0;
  constrainedRanges_[0].end = static_cast<int>(interpolatedArray2D_.rows());
  constrainedRanges_[1].begin = 0;
  constrainedRanges_[1].end = static_cast<int>(interpolatedArray2D_.cols());
  resultSet_ = false;
  extractMapIter_ = extractMap_.begin();
}


float NetCDFInterpolator::getResult() {
  // Fetch the result
  if (resultSet_) {
    // This was derived from linear interpolation so return it.
    resetExtract();
    return result_;
  }

  int sizeDim0 = constrainedRanges_[0].end - constrainedRanges_[0].begin;
  int sizeDim1 = constrainedRanges_[1].end - constrainedRanges_[1].begin;
  if (sizeDim0 != 1 || sizeDim1 != 1) {
    throw eckit::Exception("Previous calls to extract() have failed to identify "
                           "a single value to return.", Here());
  }
  float res = interpolatedArray2D_(constrainedRanges_[0].begin, constrainedRanges_[1].begin);
  resetExtract();
  return res;
}


void NetCDFInterpolator::update(const std::string &key, ioda::Variable var) {
  if (var.isA<int>()) {
    coordsVals_.emplace(key, var.readAsVector<int>());
  } else if (var.isA<float>()) {
    coordsVals_.emplace(key, var.readAsVector<float>());
  } else if (var.isA<std::string>()) {
    coordsVals_.emplace(key, var.readAsVector<std::string>());
  } else {
    throw eckit::Exception("Data type not yet supported.", Here());
  }
}


std::vector<std::string> NetCDFInterpolator::fetchDimNameMapping(
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
      ASSERT(item.size() == 1);
      dimnames.emplace_back(item[0].first);
    }
  }
  return dimnames;
}


void NetCDFInterpolator::load(const std::string &filepath,
                              const std::string &interpolatedArrayGroup) {
  // Create backend using an hdf5 file for reading.
  //   see 05b-ObsGroupAppend.cpp for reference
  ioda::Engines::BackendNames backendName;
  backendName = ioda::Engines::BackendNames::Hdf5File;
  ioda::Engines::BackendCreationParameters backendParams;
  backendParams.fileName = filepath;
  backendParams.action = ioda::Engines::BackendFileActions::Open;
  backendParams.openMode = ioda::Engines::BackendOpenModes::Read_Only;
  ioda::Group group = ioda::Engines::constructBackend(backendName, backendParams);

  // All our coords from the file.
  std::vector<std::string> vars = group.vars.list();

  // Find the array to be interpolated
  const std::string suffix = "@" + interpolatedArrayGroup;
  auto it = std::find_if(
    vars.begin(), vars.end(),
    [&suffix](const std::string& str) { return str.find(suffix) != std::string::npos; });
  if (it == vars.end()) {
    throw eckit::Exception(
      "No variable called ..." + suffix + " found within the NetCDF file", Here());
  }
  interpolatedArrayName_ = *(it);

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
  const ioda::Variable &interpolatedArrayCoord = coords[interpolatedArrayName_];
  // The variable to be interpolated shouldn't be a dimension mapping
  ASSERT(interpolatedArrayCoord.isDimensionScale() == false);
  // Fetch the dimension mapping names for this array
  interpolatedArrayDimnames_ = fetchDimNameMapping(interpolatedArrayCoord, interpolatedArrayName_,
                                                   lcoords);
  const ioda::Dimensions &dimdim = interpolatedArrayCoord.getDimensions();

  // Does the array to be interpolated represent a full covariance matrix (or a stack of them)?
  // If so, we'll need to extract diagonals.
  // --------------------------------
  std::string isCovariant {"false"};
  if (interpolatedArrayCoord.atts.exists("full")) {
    interpolatedArrayCoord.atts["full"].read(isCovariant);
  }

  // Fetch dimension mappings for our NetCDF variables.
  // --------------------------------
  // Remove the array to be interpolated from our list of coords
  vars.erase(std::remove(vars.begin(), vars.end(), interpolatedArrayName_), vars.end());
  for (const std::string& varName : vars) {
    const ioda::Variable &variable = coords[varName];
    std::vector<std::string> dimnames = fetchDimNameMapping(variable, varName, lcoords);

    // Load our variables data and keep hold of it along with relevant information.
    update(varName, variable);

    // Create a mapping between name and array dimension and vise versa.
    std::vector<int> ddim = getDimMapping(dimnames);
    for (int dim : ddim) {
      dim2CoordMapping_[dim].emplace_back(varName);
    }
    coord2DimMapping_[varName] = ddim;
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
      interpolatedArray2D_.resize(dimdim.dimsCur[1], dimdim.dimsCur[2]);

      for (int i = 0; i < data.dimension(0); i++) {
        for (int k = 0; k < data.dimension(2); k++) {
          interpolatedArray2D_(i, k) = data(i, i, k);
        }
      }
    } else if (dimdim.dimensionality == 2) {
      Eigen::ArrayXXf data;
      interpolatedArrayCoord.readWithEigenRegular(data);
      interpolatedArray2D_.resize(dimdim.dimsCur[1], 1);

      for (int i = 0; i < data.rows(); i++) {
        interpolatedArray2D_(i, 0) = data(i, i);
      }
    } else {
      throw eckit::Exception("Expecting 3D or 2D array for error covariance.", Here());
    }

    // Update our dimension mapping after collapsing a dimension.
    // --------------------------------
    // Remove the outermost dimension mapping of our data array.
    interpolatedArrayDimnames_.erase(interpolatedArrayDimnames_.begin());
    // Now update all effected coordinates (those who vary with dim0).
    std::vector<std::string> coordsDim0 = dim2CoordMapping_[dimCollapse];
    if (coordsDim0.size() > 1) {
      throw eckit::Exception(
        "More than one variable varies over the outermost dimension that we are collapsing.",
        Here());
    }
    coord2DimMapping_.erase(coordsDim0[0]);
    dim2CoordMapping_.clear();
    for (auto &coord : coord2DimMapping_) {
      // Shift remaining dimension mappings since the collapse
      for (int &dim : coord.second) {
        dim--;
        dim2CoordMapping_[dim].emplace_back(coord.first);  // update our reverse lookup
      }
    }
  } else if ((dimdim.dimensionality == 2) || (dimdim.dimensionality == 1)) {
    // This is just the diagonals - read directly into our container
    // --------------------------------
    interpolatedArrayCoord.readWithEigenRegular(interpolatedArray2D_);
  } else {
    throw eckit::Exception("The array to be interpolated has an unsupported number of dimensions.",
                           Here());
  }
}


std::vector<int> NetCDFInterpolator::getDimMapping(const std::vector<std::string> &dimnames) {
  std::vector<int> dimIndex;
  // Loop over dimension names provided
  for (const std::string & dimname : dimnames) {
    int ind = -1;
    for (std::string obsDimname : interpolatedArrayDimnames_) {
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


}  // namespace ufo
