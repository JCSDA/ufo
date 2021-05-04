/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <boost/make_unique.hpp>

#include "eckit/utils/StringTools.h"

#include "ioda/Misc/SFuncs.h"

#include "ufo/utils/dataextractor/DataExtractorCSVBackend.h"
#include "ufo/utils/dataextractor/DataExtractorInput.h"
#include "ufo/utils/dataextractor/DataExtractorNetCDFBackend.h"
#include "ufo/utils/NetCDFInterpolator.h"

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


void NetCDFInterpolator::load(const std::string &filepath,
                              const std::string &interpolatedArrayGroup) {
  std::unique_ptr<DataExtractorBackend> backend = createBackendFor(filepath);
  DataExtractorInput input = backend->loadData(interpolatedArrayGroup);
  coord2DimMapping_ = std::move(input.coord2DimMapping);
  dim2CoordMapping_ = std::move(input.dim2CoordMapping);
  coordsVals_ = std::move(input.coordsVals);
  interpolatedArray2D_ = std::move(input.payloadArray);
}

std::unique_ptr<DataExtractorBackend> NetCDFInterpolator::createBackendFor(
    const std::string &filepath) {
  const std::string lowercasePath = eckit::StringTools::lower(filepath);
  if (eckit::StringTools::endsWith(lowercasePath, ".nc") ||
      eckit::StringTools::endsWith(lowercasePath, ".nc4"))
    return boost::make_unique<DataExtractorNetCDFBackend>(filepath);
  else if (eckit::StringTools::endsWith(lowercasePath, ".csv"))
    return boost::make_unique<DataExtractorCSVBackend>(filepath);
  else
    throw eckit::BadValue("File '" + filepath + "' has an unrecognized extension. "
                          "The supported extensions are .nc, .nc4 and .csv", Here());
}


void NetCDFInterpolator::sort() {
  Eigen::ArrayXXf sortedArray = interpolatedArray2D_;
  nextCoordToExtractBy_ = coordsToExtractBy_.begin();

  for (size_t dim = 0; dim < dim2CoordMapping_.size(); ++dim) {
    // Reorder coordinates
    for (auto &coord : dim2CoordMapping_[dim]) {
      auto &coordVal = coordsVals_[coord];
      SortVisitor visitor(splitter_[dim]);
      boost::apply_visitor(visitor, coordVal);
    }
    // Reorder the array to be interpolated
    if (dim == 0) {
      int ind = -1;
      for (const auto &group : splitter_[dim].groups()) {
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
    } else if (dim == 1) {
      int ind = -1;
      for (const auto &group : splitter_[dim].groups()) {
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
  // Map any names of the form var@Group to Group/var
  const std::string canonicalVarName = ioda::convertV1PathToV2Path(varName);

  const CoordinateValues &coordVal = coordsVals_.at(canonicalVarName);
  const int dimIndex = coord2DimMapping_.at(canonicalVarName);

  SortUpdateVisitor visitor(splitter_[static_cast<size_t>(dimIndex)]);
  boost::apply_visitor(visitor, coordVal);

  // Update our map between coordinate (variable) and interpolation/extract method
  coordsToExtractBy_.emplace_back(Coordinate{varName, coordVal, method, dimIndex});
}


void NetCDFInterpolator::resetExtract() {
  constrainedRanges_[0].begin = 0;
  constrainedRanges_[0].end = static_cast<int>(interpolatedArray2D_.rows());
  constrainedRanges_[1].begin = 0;
  constrainedRanges_[1].end = static_cast<int>(interpolatedArray2D_.cols());
  resultSet_ = false;
  nextCoordToExtractBy_ = coordsToExtractBy_.begin();
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


}  // namespace ufo
