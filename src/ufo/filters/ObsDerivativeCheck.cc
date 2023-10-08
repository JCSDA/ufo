/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ObsDerivativeCheck.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Sphere.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsDerivativeCheck::ObsDerivativeCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                               std::shared_ptr<ioda::ObsDataVector<int> > flags,
                               std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr),
    parameters_(parameters)
{
  oops::Log::debug() << "ObsDerivativeCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsDerivativeCheck::~ObsDerivativeCheck() {}

// -----------------------------------------------------------------------------

void get_locs(const std::vector<std::size_t> & rSort, const size_t & i1, const size_t & i2,
              const size_t & ilocs, size_t & ii1, size_t & ii2)
// rSort - the vector of integers of indices of the sorted ObsDataVector
//         used to get the size of the sorted record
// i1 and i2 are by default 0 but can be defined in YAML to assign fixed indices
//         for the derivative computed for each point in the record
//         if i1 and i2 are both zero then a local derivative will be computed and indices
//         will be computed by this function
// ii1 and ii2 are the output indices from the function to use in computation of derivatives
//         ii1=i1 and ii2=i2 for when i1 and i2 != 0
//         ii1 and ii2 are computed below for all other cases
{
  if ( rSort.size() == 1 ) {
    // set both ii1 and ii2 to 0 which will be used to indicate derivative is zero
    ii1 = 0;
    ii2 = 0;
  } else if ( rSort.size() == 2 ) {
    // set ii1 to 0 and ii2 to 1 because there are only 2 points in this group
    ii1 = 0;
    ii2 = 1;
  } else {
    if ( (i1 == 0) && (i2 == 0) ) {
      // this means local derivatives will be computed for each point
      if ( ilocs == 0 ) {
        // special case for boundary condition 1
        ii1 = 0;
        ii2 = 1;
      } else if ( ilocs == rSort.size()-1 ) {
        // special case for boundary condition 2
        ii1 = ilocs-1;
        ii2 = ilocs;
      } else {
        // define local derivative normally as
        // (y(i+1) - y(i-1)) / (x(i+1) - x(i-1))
        ii1 = ilocs-1;
        ii2 = ilocs+1;
      }
    } else {
      // compute derivative for each point to be from indices provided in YAML
      // note this may need to be more complex later to compute the indices desired
      // set index to the last value if -1 defined in YAML
      (i1 == -1) ? ii1 = rSort.size()-1 : ii1 = i1;
      (i2 == -1) ? ii2 = rSort.size()-1 : ii2 = i2;
      // set index to the last value if YAML greater than indices
      (i1 > rSort.size()-1) ? ii1 = rSort.size()-1 : ii1 = i1;
      (i2 > rSort.size()-1) ? ii2 = rSort.size()-1 : ii2 = i2;
    }
  }
}

// -----------------------------------------------------------------------------

void ObsDerivativeCheck::applyFilter(const std::vector<bool> & apply,
                                 const Variables & filtervars,
                                 std::vector<std::vector<bool>> & flagged) const {
  const float missing = util::missingValue<float>();
  const double radiusEarth = Constants::mean_earth_rad*1000.0;

  // first we want to get the config of the two vars to use in computing the derivative
  const std::string strInd_ = parameters_.independent;
  const std::string strDep_ = parameters_.dependent;
  // specified indices to use if not local derivatives
  const size_t i1 = parameters_.i1;
  const size_t i2 = parameters_.i2;
  size_t ii1 = 0;  // indices that will be used to compute derivativs
  size_t ii2 = 0;
  // min/max value setup
  const float minddx = parameters_.minvalue;
  const float maxddx = parameters_.maxvalue;

  oops::Log::debug() << "ObsDerivativeCheck: Independent Var = " << strInd_ << std::endl;
  oops::Log::debug() << "ObsDerivativeCheck: Dependent Var = " << strDep_ << std::endl;

  // now grab the vars and compute the derivative
  const size_t nlocs_ = obsdb_.nlocs();
  std::vector<float> dydx(nlocs_);

  // different options depending on the independent variable
  // the generic case:
  //    dy/dx = y(ii2)-y(ii1) / x(ii2) - x(ii1)
  // when x is datetime:
  //    time is computed in seconds and denominator is handled differently
  //    when y is distance:
  //       longitude and latitude are used to compute great circle distances
  //       for numerator
  // when x is distance:
  //       longitude and latitude are used to compute great circle distances
  //       for denominator
  if ( strInd_ == "dateTime" ) {  // special case for datetime
      std::vector<util::DateTime> varIndT_(nlocs_);
      obsdb_.get_db("MetaData", strInd_, varIndT_);
      ioda::ObsSpace::RecIdxIter irec;
      if ( strDep_ == "distance" ) {  // special case for moving obs
        ioda::ObsDataVector<float> varDepX_(obsdb_, "longitude", "MetaData");
        ioda::ObsDataVector<float> varDepY_(obsdb_, "latitude", "MetaData");
        for ( irec = obsdb_.recidx_begin(); irec != obsdb_.recidx_end(); ++irec ) {
          std:: size_t rNum = obsdb_.recidx_recnum(irec);
          std::vector<std::size_t> rSort = obsdb_.recidx_vector(irec);
          for (size_t ilocs = 0; ilocs < rSort.size(); ++ilocs) {
            get_locs(rSort, i1, i2, ilocs, ii1, ii2);
            if ( ii1 == ii2 ) {
              // no derivative if the indices are the same
              dydx[rSort[ilocs]] = 0.0;
            } else {
              // compute distance for each lat lon pair
              eckit::geometry::Point2 Ob1(varDepX_["longitude"][rSort[ii2]],
                                          varDepY_["latitude"][rSort[ii2]]);
              eckit::geometry::Point2 Ob2(varDepX_["longitude"][rSort[ii1]],
                                          varDepY_["latitude"][rSort[ii1]]);
              float dist_m = static_cast<float>(eckit::geometry::Sphere::distance(radiusEarth,
                                                                                 Ob1, Ob2));
              // compute derivative
              dydx[rSort[ilocs]] = dist_m /
                                   (varIndT_[rSort[ii2]] - varIndT_[rSort[ii1]]).toSeconds();
            }
          }
        }
      } else {
        ioda::ObsDataVector<float> varDep_(obsdb_, strDep_, "MetaData");
        for ( irec = obsdb_.recidx_begin(); irec != obsdb_.recidx_end(); ++irec ) {
          std:: size_t rNum = obsdb_.recidx_recnum(irec);
          std::vector<std::size_t> rSort = obsdb_.recidx_vector(irec);
          for (size_t ilocs = 0; ilocs < rSort.size(); ++ilocs) {
            get_locs(rSort, i1, i2, ilocs, ii1, ii2);
            if ( ii1 == ii2 ) {
              // no derivative if the indices are the same
              dydx[rSort[ilocs]] = 0.0;
            } else {
              // compute derivative
              dydx[rSort[ilocs]] = (varDep_[strDep_][rSort[ii2]] -
                                    varDep_[strDep_][rSort[ii1]]) /
                                   (varIndT_[rSort[ii2]] -
                                    varIndT_[rSort[ii1]]).toSeconds();
            }
          }
        }
     }
  // end if strInd_ == "datetime"
  } else if ( strInd_ == "distance" ) {  // convert lat/lon change to distance
      ioda::ObsDataVector<float> varDep_(obsdb_, strDep_, "MetaData");
      ioda::ObsDataVector<float> varIndX_(obsdb_, "longitude", "MetaData");
      ioda::ObsDataVector<float> varIndY_(obsdb_, "latitude", "MetaData");
      ioda::ObsSpace::RecIdxIter irec;
      for ( irec = obsdb_.recidx_begin(); irec != obsdb_.recidx_end(); ++irec ) {
        std:: size_t rNum = obsdb_.recidx_recnum(irec);
        std::vector<std::size_t> rSort = obsdb_.recidx_vector(irec);
        for (size_t ilocs = 0; ilocs < rSort.size(); ++ilocs) {
          get_locs(rSort, i1, i2, ilocs, ii1, ii2);
          if ( ii1 == ii2 ) {
            // no derivative if the indices are the same
            dydx[rSort[ilocs]] = 0.0;
          } else {
            // compute distance for each lat lon pair
            eckit::geometry::Point2 Ob1(varIndX_["longitude"][rSort[ii2]],
                                        varIndY_["latitude"][rSort[ii2]]);
            eckit::geometry::Point2 Ob2(varIndX_["longitude"][rSort[ii1]],
                                        varIndY_["latitude"][rSort[ii1]]);
            float dist_m = static_cast<float>(eckit::geometry::Sphere::distance(radiusEarth,
                                                                               Ob1, Ob2));
            dydx[rSort[ilocs]] = (varDep_[strDep_][rSort[ii2]] - varDep_[strDep_][rSort[ii1]]) /
                                 dist_m;
          }
        }
      }
  // end if strInd_ == "distance"
  } else {  // standard case where independent var is not datetime or distance
      ioda::ObsDataVector<float> varDep_(obsdb_, strDep_, "MetaData");
      ioda::ObsDataVector<float> varInd_(obsdb_, strInd_, "MetaData");
      ioda::ObsSpace::RecIdxIter irec;
      for ( irec = obsdb_.recidx_begin(); irec != obsdb_.recidx_end(); ++irec ) {
        std:: size_t rNum = obsdb_.recidx_recnum(irec);
        std::vector<std::size_t> rSort = obsdb_.recidx_vector(irec);
        for (size_t ilocs = 0; ilocs < rSort.size(); ++ilocs) {
          get_locs(rSort, i1, i2, ilocs, ii1, ii2);
          if ( ii1 == ii2 ) {
            dydx[rSort[ilocs]] = 0.0;
          } else {
            dydx[rSort[ilocs]] = (varDep_[strDep_][rSort[ii2]] - varDep_[strDep_][rSort[ii1]]) /
                                 (varInd_[strInd_][rSort[ii2]] - varInd_[strInd_][rSort[ii1]]);
          }
        }
      }
  }

  // determine if the derivative is outside the specified range
  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (apply[jobs]) {
        if (dydx[jobs] > maxddx && maxddx != missing) flagged[jv][jobs] = flagged[jv][jobs] = true;
        if (dydx[jobs] < minddx && minddx != missing) flagged[jv][jobs] = flagged[jv][jobs] = true;
      }
    }
  }
}

// -----------------------------------------------------------------------------

void ObsDerivativeCheck::print(std::ostream & os) const {
  os << "ObsDerivativeCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
