/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_OBSERRORNETCDFHANDLER_H_
#define UFO_UTILS_OBSERRORNETCDFHANDLER_H_

#include <algorithm>           // sort
#include <limits>              // std::numeric_limits
#include <list>                // list
#include <string>
#include <unordered_map>
#include <utility>             // pair
#include <vector>

#include "boost/variant.hpp"
#include "Eigen/Dense"         // Eigen Arrays and Matrices

#include "eckit/exception/Exceptions.h"
#include "ioda/Engines/Factory.h"
#include "ioda/Group.h"
#include "ioda/ObsGroup.h"
#include "oops/util/Logger.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {


enum class InterpMethod {
  EXACT, NEAREST, LINEAR
};


/// \brief Observation error NetCDF handler, performing extract and interpolation to derive the
/// observation error (error standard deviation).
/// \details The NetCDF format is as follows:
///   - The observation error data variable from the file will be identified by a name
///     with @ErrorVariance suffix in the varname.
///   - Any number of variables are supported of types: int; float and string (inc. DateTime).
///   - Dimensionality of the observation error is currently limited to 1D, 2D or 3D.
///     - 1D or 2D corresponding to the diagonal only variance.
///     - 2D or 3D corresponding to the error covariance matrix defined.
///   - Variables (coordinates) are expected to be 1D only.
///   - An error covariance matrix is defined by the presence of a local attribute, full = "true" on
///     the error variance variable, otherwise the file is assumed to describe diagonals only.
///     - For an error covariance, the first dimension is collapsed as a consequence of extracting
///       the diagonals.
///
/// Here is an example NetCDF metadata describing variance (diagonals only) with multi-variant:
/// \code
/// dimensions:
///   air_pressure@MetaData = 10 ;
///   index = 5 ;
/// variables:
///   int air_temperature@ErrorVariance(air_pressure@MetaData, index) ;
///     air_temperature@ErrorVariance:coordinates = "observation_type@MetaData index" ;
///   int air_pressure@MetaData(air_pressure@MetaData) ;
///   int index(index) ;
///   int observation_type@MetaData(index) ;
/// \endcode
///
/// Here is an example NetCDF metadata describing a multi-variant error covariance matrix:
/// \code
/// dimensions:
///   channel_number = 10 ;
///   channel_number@MetaData = 10 ;
///   index = 8 ;
/// variables:
///   float air_temperature@ErrorVariance(channel_number, channel_number@MetaData, index) ;
///     air_temperature@ErrorVariance:coordinates = "latitude_band@MetaData \
///         processing_center@MetaData satellite_id@MetaData" ;
///   int index(index) ;
///   int channel_number(channel_number) ;
///   int channel_number@MetaData(channel_number@MetaData) ;
///   int latitude_band@MetaData(index) ;
///   int processing_center@MetaData(index) ;
///   int satellite_id@MetaData(index) ;
/// \endcode
///
/// Notice how a channel_number describes the outsermost two dimensions.  Since the very outermost
/// dimension is collapsed, it's name does not have to match anything from our ObsSpace as it is
/// discarded after extracting the diagonals.
///
/// Here is a summary of particulars to the interpolation algorithms available:
/// - Linear interpolation takes place between values of variance.  We then determine
///   its square root in order to return the standard deviation observation error (@ObsError).
/// - Nearest neighbour 'interpolation' chooses the **first** nearest value to be found in the case
///   of equidistant neighbours of different values. The algorithm will then return the one or more
///   locations corresponding to this nearest value.  Let's illustrate by extracting the nearest
///   neighbours to 1.5:
///
///   [1, 1, 2, 3, 4, 5]
///
///   Both 1 and 2 are equidtstant, but we take the first found equidistant neighbour (1), then
///   return all indices matching this neighbour (indices 0 and 1 in this case).
class ObsErrorNetCDFHandler
{
 public:
  /// \brief Create an object representing a NetCDF handle to an observation error variance.
  /// \details This object is capable of sorting the data from this file; extracting the relevant
  /// values for a given observation as well as performing linear interpolation to finally derive
  /// our observation error value.  See class documentation for further details on the expected
  /// structure of the NetCDF.
  /// \param[in] filepath is the file-path to the observation variance file.
  explicit ObsErrorNetCDFHandler(const std::string &filepath);

  /// \brief Update the instruction on how to sort the data for the provided variable name.
  /// \details This works iteratively by further splitting the RecurssiveSplitter sub-groups
  /// according to the variable name provided.  By this, it is possible to sort the data in
  /// such a way to ensure that extraction always results in 1 contiguous chunk.
  /// No data or coordinates are actually physically sorted yet.  Special treatment for float
  /// type variables, where this is used to sort each of the sub-groups.
  /// \param[in] varName is the variable name (NetCDF coordinate name) that we will be sorting.
  /// \internal This member function call corresponds to a RecursiveSplitter.groupBy call, useful
  /// to sort according to nearesr/exact match variables.  In the special case of float type,
  /// RecursiveSplitter.sortGroupsBy is used.
  void scheduleSort(const std::string &varName, const InterpMethod &method);

  /// \brief Finalise the sort, sorting each of the coordinates describing the variance, as well
  /// as the variance itself.
  /// \details Utilising the instructions provided by the user calling the scheduleSort() member
  /// function, we now physically sort the variance array itself along with all coordinates which
  /// describe it.
  /// \internal Applies the RecurssiveSplitter object and neccessarily creates copies to achieve
  /// this sort.
  void sort();

  /// \brief Perform extract, given an observation value for the coordinate associated with this
  /// extract iteration.
  /// \details Calls the relevant extract methood (linear, nearest or exact), corresponding to the
  /// coordinate associated with this extract iteration (along with the associated interpolation
  /// method).  The extract then utilises the observation value ('obVal') to perform this
  /// extraction.
  /// \param[in] obVal is the observation value used for the extract operation.
  template<typename T>
  void extract(const T &obVal) {
    if (extractMapIter_ == extractMap_.end()) {
      throw eckit::UserError("Too many extract() calls made for the expected number of variables.",
                             Here());
    }

    // Fetch corresponding NetCDF var data and the dimension mapping
    const std::string &varName = extractMapIter_->first;
    switch (extractMapIter_->second) {
      case InterpMethod::EXACT:
        exactMatch(varName, obVal);
        break;
      case InterpMethod::LINEAR:
        obErrorVal_ = getObsErrorValue(varName, obVal);
        obErrorValSet_ = true;
        break;
      case InterpMethod::NEAREST:
        nearestMatch(varName, obVal);
        break;
      default:
        throw eckit::UserError("Only Linear, nearest and exact interpolation methods supported.",
                               Here());
    }
    ++extractMapIter_;
  }

  /// \brief Fetch the observation error value.
  /// \details This will only be succesful if previous calls to extract() have resulted in a single
  /// value to return.
  float getObsErrorValue();

 private:
  /// \brief Fetch the observation error variance at a 'location'.
  /// \details Fetching the observation error variance at a 'location' (obVal) is to linearly
  /// interpolate to the location along the coordinate of specified name which describes the
  /// observation variance.
  /// \param[in] varName is the variable name (NetCDF coordinate name) that we will be using to
  /// perform the interpolation.
  /// \param[in] obVal is the interpolation location.
  template<typename T>
  float getObsErrorValue(const std::string &varName, const T &obVal) {
    // Fetch corresponding NetCDF var data and the dimension mapping
    int dimIndex;
    const std::vector<T> &ncdfVal = getVarDim<T>(varName, dimIndex);

    // Sanity check constraint
    int sizeDim0 = constrainedRanges_[0].end - constrainedRanges_[0].begin;
    int sizeDim1 = constrainedRanges_[1].end - constrainedRanges_[1].begin;
    if ((dimIndex == 1 && !(sizeDim1 > 1 && sizeDim0 == 1)) ||
        !(sizeDim0 > 1 && sizeDim1 == 1)) {
      throw eckit::Exception("Linear interpolation failed - data must be 1D.", Here());
    }

    // Constrain our index range in the relevant dimension.
    const Range &range = constrainedRanges_[static_cast<size_t>(dimIndex)];

    if ((obVal > ncdfVal[range.end - 1]) || (obVal < ncdfVal[range.begin])) {
      throw eckit::Exception("Linear interpolation failed, value is beyond grid extent."
                             "No extrapolation supported.",
                             Here());
    }
    // Find first index of ncdfVal >= obVal
    int nnIndex = std::lower_bound(ncdfVal.begin() + range.begin,
                                     ncdfVal.begin() + range.end,
                                     obVal) - ncdfVal.begin();

    // Determine upper or lower indices from this
    if (ncdfVal[nnIndex] == obVal) {
      // No interpolation required (is equal)
      float res;
      if (dimIndex == 1) {
        res = static_cast<float>(obserrData2D_(constrainedRanges_[0].begin, nnIndex));
      } else {
        res = static_cast<float>(obserrData2D_(nnIndex, constrainedRanges_[1].begin));
      }
      return res;
    }
    // Linearly interpolate between these two indices.
    auto zUpper = *(obserrData2D_.data());
    auto zLower = *(obserrData2D_.data());
    if (dimIndex == 1) {
      zUpper = obserrData2D_(constrainedRanges_[0].begin, nnIndex);
      zLower = obserrData2D_(constrainedRanges_[0].begin, nnIndex-1);
    } else {
      zUpper = obserrData2D_(nnIndex, constrainedRanges_[1].begin);
      zLower = obserrData2D_(nnIndex-1, constrainedRanges_[1].begin);
    }
    float res = ((static_cast<float>(obVal - ncdfVal[nnIndex-1]) /
                  static_cast<float>(ncdfVal[nnIndex] - ncdfVal[nnIndex-1])) *
                 (zUpper - zLower)) + zLower;
    return res;
  }

  float getObsErrorValue(const std::string &varName, const std::string &obVal) {
    throw eckit::UserError("VarName: " + varName +
                           " - linear interpolation not compatible with string type.", Here());
  }

  /// \brief Update our extract constraint based on a nearest match against the specified
  /// coordinate associated with the observation variance.
  /// \details
  ///
  /// Method:
  /// - Find **first** discovered nearest value in our loop.
  /// - Determine which indices match this nearest value.
  ///   (more than one index could have this one value).
  ///
  ///   [1, 1, 2, 3, 4, 5]
  ///
  /// Nearest neighbour extraction of “1”, has more than one neighbour.
  /// That is, more than one index with the same value have the same distance:
  ///
  ///   [1, 1]  i.e. range=(0, 2)
  ///
  /// - Note that an alternative implementation could consider equidistant
  ///   values, though it was decided this was not desirable behaviour:
  ///
  ///   [1, 1, 2, 3, 4, 5]
  ///
  /// Nearest neighbour extraction of “1.5” could be then considered to have 3
  /// equidistant neighbours (1, 1, 2).  That is, two different values with the
  /// same distance.
  ///
  /// [1, 1, 2] i.e. range=(0, 3)
  ///
  /// \param[in] varName is the variable name (NetCDF coordinate name) that we match against.
  /// \param[in] obVal is the value to match, against the NetCDF coordinate of name 'varName'.
  template<typename T>
  void nearestMatch(const std::string &varName, const T &obVal) {
    // Fetch corresponding NetCDF var data and the dimension mapping
    int dimIndex;
    const std::vector<T> &ncdfVal = getVarDim<T>(varName, dimIndex);

    // Constrain our index range in the relevant dimension.
    Range &range = constrainedRanges_[static_cast<size_t>(dimIndex)];

    // Find first index of ncdfVal >= obVal
    int nnIndex = std::lower_bound(ncdfVal.begin() + range.begin,
                                   ncdfVal.begin() + range.end,
                                   obVal) - ncdfVal.begin();
    if (nnIndex >= range.end) {
      nnIndex = range.end - 1;
    }

    // Now fetch the nearest neighbour index (lower index prioritised for different values with
    // same distance)
    T dist = std::abs(ncdfVal[nnIndex] - obVal);
    if ((ncdfVal[nnIndex] > obVal) && (nnIndex > range.begin) &&
        (std::abs(ncdfVal[nnIndex - 1] - obVal) <= dist))
      nnIndex--;

    // Now find **same value** equidistant neighbours
    auto bounds = std::equal_range(ncdfVal.begin() + range.begin,
                                   ncdfVal.begin() + range.end,
                                   ncdfVal[nnIndex]);
    range = {static_cast<int>(bounds.first - ncdfVal.begin()),
             static_cast<int>(bounds.second - ncdfVal.begin())};
    oops::Log::debug() << "Nearest match; name: " << varName << " range: " <<
      range.begin << "," << range.end << std::endl;
  }

  void nearestMatch(const std::string &varName, const std::string &obVal) {
    throw eckit::UserError("Nearest match not compatible with string type.", Here());
  }

  /// \brief Update our extract constraint based on an exact match against the specified coordinate
  /// associated with the observation variance.
  /// \param[in] varName is the variable name (NetCDF coordinate name) that we match against.
  /// \param[in] obVal is the value to match, against the NetCDF coordinate of name 'varName'.
  template<typename T>
  void exactMatch(const std::string &varName, const T &obVal) {
    // Fetch corresponding NetCDF var data and the dimension mapping
    int dimIndex;
    const std::vector<T> &ncdfVal = getVarDim<T>(varName, dimIndex);

    // Constrain our index range in the relevant dimension.
    Range &range = constrainedRanges_[static_cast<size_t>(dimIndex)];

    // Find the first and last matching index
    auto bounds = std::equal_range(ncdfVal.begin() + range.begin,
                                   ncdfVal.begin() + range.end,
                                   obVal);
    range = {static_cast<int>(bounds.first - ncdfVal.begin()),
             static_cast<int>(bounds.second - ncdfVal.begin())};

    if (range.begin == range.end) {
      throw eckit::Exception(
        "No match found for exact match extraction of " + varName,
        Here());
    }
    oops::Log::debug() << "Exact match; name: " << varName << " range: " <<
      range.begin << "," << range.end << std::endl;
  }

  /// \brief Reset the extraction range for this object.
  /// \details Each time an exactMatch or nearestMatch call is made for one or more variable,
  /// the extraction range is further constrained to match our updated match conditions.  After
  /// the final 'extract' is made (i.e. an observation error is derived) it is desirable to reset
  /// the extraction range by calling this method.
  /// \internal This is called by the getObsErrorValue member functions just before returning the
  /// observation error.
  void resetExtract();

  /// \brief Fetch the NetCDF coordinate data of specified name.
  /// \details We cache this data so as not to require performing multiple loads from disk.
  template <typename T>
  T& get(const std::string &key) {
    try {
      return boost::get<T> (coordsVals_[key]);
    } catch (boost::bad_get) {
      throw eckit::BadParameter("Unable to find coordinate with this type", Here());
    }
  }

  /// \brief Add the specified variable (from NetCDF) to our cache for later retrieval.
  void update(const std::string &key, ioda::Variable var);

  /// \brief Convenience helper function to return us the netCDF data along with the dimension
  /// mapping
  template<typename T>
  const std::vector<T> & getVarDim(const std::string &varName, int &dimIndex) {
    const auto &ncdfVal = get<std::vector<T>>(varName);
    std::vector<int> dimIndices = coord2DimMapping_.at(varName);
    if (dimIndices.size() != 1) {
      throw eckit::Exception(
        "Expecting a coordinate mapping to a single dimension",
        Here());
    }
    dimIndex = dimIndices[0];
    return ncdfVal;
  }

  /// \brief Helper function to determine dimension index from a dimension name
  std::vector<int> getDimMapping(const std::vector<std::string> &dimnames);

  /// \brief Helper function for determining the dimension mapping names for a variable
  std::vector<std::string> fetchDimNameMapping(
      const ioda::Variable &variable,
      const std::string &varName,
      const std::list<std::pair<std::string, ioda::Variable>> &coordinates);

  /// \brief Perform any load from the NetCDF (the variance data, coordinates etc.)
  void load();

  std::string filepath_;
  ioda::Group group_;
  std::string obsErrorName_;

  // Dimension names of our observation error variance.
  std::vector<std::string> obErrorDimnames_;
  // Object represent the extraction range in both dimensions.
  struct Range {int begin, end;};
  std::array<Range, 2> constrainedRanges_;

  // Container for holding our coordinates of numerous types from NetCDF (cache).
  std::unordered_map<std::string,
                     boost::variant<std::vector<int>,
                                    std::vector<float>,
                                    std::vector<std::string>
                                   >
                    > coordsVals_;
  // Our variance
  Eigen::ArrayXXf obserrData2D_;
  float obErrorVal_;
  bool obErrorValSet_;
  // Container for re-ordering our data
  std::vector<ufo::RecursiveSplitter> splitter_;

  // Dimension to variable mapping (+vice versa)
  std::unordered_map<std::string, std::vector<int>> coord2DimMapping_;
  std::unordered_map<int, std::vector<std::string>> dim2CoordMapping_;

  // Expected coord name and interpolation method mapping.
  std::vector<std::pair<std::string, InterpMethod>> extractMap_;
  std::vector<std::pair<std::string, InterpMethod>>::iterator extractMapIter_;
};

}  // namespace ufo

#endif  // UFO_UTILS_OBSERRORNETCDFHANDLER_H_
