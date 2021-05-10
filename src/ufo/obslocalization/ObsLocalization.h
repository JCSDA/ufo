/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSLOCALIZATION_H_
#define UFO_OBSLOCALIZATION_OBSLOCALIZATION_H_

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/util/Earth.h"

#include "eckit/config/Configuration.h"
#include "eckit/container/KDTree.h"
#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"
#include "eckit/geometry/UnitSphere.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/Printable.h"

#include "ufo/obslocalization/ObsLocParameters.h"

namespace ufo {

/// Horizontal Box car observation space localization
template<class MODEL>
class ObsLocalization: public util::Printable {
  typedef typename MODEL::GeometryIterator   GeometryIterator_;

 public:
  struct TreeTrait {
    typedef eckit::geometry::Point3 Point;
    typedef double                  Payload;
  };
  typedef eckit::KDTreeMemory<TreeTrait> KDTree;
  ObsLocalization(const eckit::Configuration &, const ioda::ObsSpace &);

  /// compute localization and save localization values in \p obsvector and
  /// localization flags (1: outside of localization; 0: inside localization area)
  /// in \p outside
  virtual void computeLocalization(const GeometryIterator_ &, ioda::ObsDataVector<int> & outside,
                           ioda::ObsVector & obsvector) const;

  const std::vector<int> & localobs() const {return localobs_;}
  const std::vector<double> & horizontalObsdist() const {return obsdist_;}
  const ObsLocParameters & localizationOptions() const {return options_;}

 private:
  ObsLocParameters options_;
  mutable std::vector<double> obsdist_;
  mutable std::vector<int> localobs_;

  void print(std::ostream &) const override;

  /// KD-tree for searching for local obs
  std::unique_ptr<KDTree> kd_;

  std::vector<float> lats_;
  std::vector<float> lons_;

  /// TODO(travis) distribution name is needed for temporary fix, should be removed eventually
  std::string distName_;
};

// -----------------------------------------------------------------------------

/*!
 * \details Creates a KDTree class member that can be used for searching for local obs
 */
template<typename MODEL>
ObsLocalization<MODEL>::ObsLocalization(const eckit::Configuration & config,
                                        const ioda::ObsSpace & obsspace)
  : lats_(obsspace.nlocs()), lons_(obsspace.nlocs())
{
  options_.deserialize(config);

  // check that this distribution supports local obs space
  // TODO(travis) this has been moved to computeLocalization as a quick fix for a bug.
  distName_ = obsspace.distribution().name();

  const size_t nlocs = obsspace.nlocs();
  // Get latitudes and longitudes of all observations.
  obsspace.get_db("MetaData", "longitude", lons_);
  obsspace.get_db("MetaData", "latitude", lats_);

  if (options_.searchMethod == SearchMethod::KDTREE) {
    kd_ = std::unique_ptr<KDTree> ( new KDTree() );
    // Define points list from lat/lon values
    typedef typename KDTree::PointType Point;
    std::vector<typename KDTree::Value> points;
    for (unsigned int i = 0; i < nlocs; i++) {
      eckit::geometry::Point2 lonlat(lons_[i], lats_[i]);
      Point xyz = Point();
      // FIXME: get geometry from yaml, for now assume spherical Earth radius.
      atlas::util::Earth::convertSphericalToCartesian(lonlat, xyz);
      double index = static_cast<double>(i);
      typename KDTree::Value v(xyz, index);
      points.push_back(v);
    }
    // Create KDTree class member from points list.
    kd_->build(points.begin(), points.end());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocalization<MODEL>::computeLocalization(const GeometryIterator_ & i,
                                            ioda::ObsDataVector<int> & outside,
                                            ioda::ObsVector & locvector) const {
  oops::Log::trace() << "ObsLocalization::computeLocalization" << std::endl;

  // check that this distribution supports local obs space
  // TODO(travis) this should be in the constructor, but currently
  //  breaks LETKF when using a split observer/solver
  if ( distName_ != "Halo" && distName_ != "InefficientDistribution" ) {
    std::string message = "Can not use ObsLocalization with distribution=" + distName_;
    throw eckit::BadParameter(message);
  }

  // clear arrays before proceeding
  localobs_.clear();
  obsdist_.clear();

  eckit::geometry::Point2 refPoint = *i;
  size_t nlocs = lons_.size();
  if ( options_.searchMethod == SearchMethod::BRUTEFORCE ) {
    oops::Log::trace() << "Local obs searching via brute force." << std::endl;

    for (unsigned int jj = 0; jj < nlocs; ++jj) {
      eckit::geometry::Point2 searchPoint(lons_[jj], lats_[jj]);
      double localDist = options_.distance(refPoint, searchPoint);
      if ( localDist < options_.lengthscale ) {
        localobs_.push_back(jj);
        obsdist_.push_back(localDist);
      }
    }
    const boost::optional<int> & maxnobs = options_.maxnobs;
    if ( (maxnobs != boost::none) && (localobs_.size() > *maxnobs ) ) {
      for (unsigned int jj = 0; jj < localobs_.size(); ++jj) {
          oops::Log::debug() << "Before sort [i, d]: " << localobs_[jj]
              << " , " << obsdist_[jj] << std::endl;
      }
      // Construct a temporary paired vector to do the sorting
      std::vector<std::pair<std::size_t, double>> localObsIndDistPair;
      for (unsigned int jj = 0; jj < obsdist_.size(); ++jj) {
        localObsIndDistPair.push_back(std::make_pair(localobs_[jj], obsdist_[jj]));
      }

      // Use a lambda function to implement an ascending sort.
      sort(localObsIndDistPair.begin(), localObsIndDistPair.end(),
           [](const std::pair<std::size_t, double> & p1,
              const std::pair<std::size_t, double> & p2){
                return(p1.second < p2.second);
              });

      // Unpair the sorted pair vector
      for (unsigned int jj = 0; jj < obsdist_.size(); ++jj) {
        localobs_[jj] = localObsIndDistPair[jj].first;
        obsdist_[jj] = localObsIndDistPair[jj].second;
      }

      // Truncate to maxNobs length
      localobs_.resize(*maxnobs);
      obsdist_.resize(*maxnobs);
    }
  } else if (nlocs > 0) {
    // Check (nlocs > 0) is needed,
    // otherwise, it will cause ASERT check fail in kdtree.findInSphere, and hang.

    oops::Log::trace() << "Local obs searching via KDTree" << std::endl;

    if ( options_.distanceType == DistanceType::CARTESIAN)
     ABORT("ObsLocalization:: search method must be 'brute_force' when using 'cartesian' distance");

    // Using the radius of the earth
    eckit::geometry::Point3 refPoint3D;
    atlas::util::Earth::convertSphericalToCartesian(refPoint, refPoint3D);
    double alpha =  (options_.lengthscale / options_.radius_earth)/ 2.0;  // angle in radians
    double chordLength = 2.0*options_.radius_earth * sin(alpha);  // search radius in 3D space

    auto closePoints = kd_->findInSphere(refPoint3D, chordLength);

    // put closePoints back into localobs_ and obsdist
    for (unsigned int jloc = 0; jloc < closePoints.size(); ++jloc) {
       localobs_.push_back(closePoints[jloc].payload());  // observation
       obsdist_.push_back(closePoints[jloc].distance());  // distance
    }

    // The obs are sorted in the kdtree call
    const boost::optional<int> & maxnobs = options_.maxnobs;
    if ( (maxnobs != boost::none) && (localobs_.size() > *maxnobs ) ) {
      // Truncate to maxNobs length
      localobs_.resize(*maxnobs);
      obsdist_.resize(*maxnobs);
    }
  }
  for (size_t jloc = 0; jloc < outside.nlocs(); ++jloc) {
    for (size_t jvar = 0; jvar < outside.nvars(); ++jvar) {
      outside[jvar][jloc] = 1;
    }
  }
  const size_t nvars = locvector.nvars();
  for (size_t jlocal = 0; jlocal < localobs_.size(); ++jlocal) {
    // obsdist is calculated at each location; need to update R for each variable
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      outside[jvar][localobs_[jlocal]] = 0;
      locvector[jvar + localobs_[jlocal] * nvars] = 1.0;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocalization<MODEL>::print(std::ostream & os) const {
  os << "ObsLocalization (box car) horizontal localization with " << options_.lengthscale
     << " lengthscale" << std::endl;
}

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSLOCALIZATION_H_
