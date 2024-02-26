/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSHORLOCALIZATION_H_
#define UFO_OBSLOCALIZATION_OBSHORLOCALIZATION_H_

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

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/ObsLocalizationBase.h"
#include "oops/util/missingValues.h"

#include "ufo/obslocalization/ObsHorLocParameters.h"
#include "ufo/ObsTraits.h"

namespace ufo {

/// Horizontal Box car observation space localization
template<class MODEL>
class ObsHorLocalization: public oops::ObsLocalizationBase<MODEL, ObsTraits> {
  typedef typename MODEL::GeometryIterator   GeometryIterator_;

 public:
  ObsHorLocalization(const eckit::Configuration &, const ioda::ObsSpace &);

  /// Compute localization and save localization values in \p locvector.
  /// Missing values indicate that observation is outside of localization.
  /// The lengthscale from ObsHorLocParameters is used.
  void computeLocalization(const GeometryIterator_ &,
                           ioda::ObsVector & locvector) const override;

 protected:
  struct LocalObs {
    /// The list of indexes for ObsVector pointing to the valid local obs.
    std::vector<int> index;

    /// The horizontal distance of each local ob from the search point.
    std::vector<double> distance;

    /// The maximum search distance that was used for this local obs search.
    double lengthscale;
  };

  /// For a given distance, returns the local observations and their distances.
  /// Intended to be called by \c computeLocalization() .
  const LocalObs getLocalObs(const GeometryIterator_ &, double lengthscale) const;

  /// Compute box car localization using the set of \p localobs and save
  /// localization values in \p locvector. Mmissing values are set for obs
  /// outside of localization.
  /// Intended to be called by \c computeLocalization() .
  virtual void localizeLocalObs(const GeometryIterator_ &,
                                ioda::ObsVector & locvector,
                                const LocalObs & localobs) const;

  /// Get the lengthscale specified in the parameters.
  double lengthscale() const {return options_.lengthscale;}

 private:
  ObsHorLocParameters options_;

  void print(std::ostream &) const override;

  /// KD-tree for searching for local obs
  struct TreeTrait {
    typedef eckit::geometry::Point3 Point;
    typedef double                  Payload;
  };
  typedef eckit::KDTreeMemory<TreeTrait> KDTree;
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
ObsHorLocalization<MODEL>::ObsHorLocalization(const eckit::Configuration & config,
                                              const ioda::ObsSpace & obsspace)
  : options_(), lats_(obsspace.nlocs()), lons_(obsspace.nlocs())
{
  options_.deserialize(config);
  // check that this distribution supports local obs space
  // TODO(travis) this has been moved to computeLocalization as a quick fix for a bug.
  distName_ = obsspace.distribution()->name();

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
void ObsHorLocalization<MODEL>::computeLocalization(const GeometryIterator_ & i,
                                                 ioda::ObsVector & locvector) const {
  oops::Log::trace() << "ObsHorLocalization::computeLocalization" << std::endl;

  // get the set of local observations using the lengthscale given in the
  // config file options.
  const LocalObs & localobs = getLocalObs(i, options_.lengthscale);

  // compute localization of those local obs. Note that since this is
  // a virtual method, it could be overriden by dervied classes
  localizeLocalObs(i, locvector, localobs);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsHorLocalization<MODEL>::localizeLocalObs(const GeometryIterator_ & i,
                                              ioda::ObsVector & locvector,
                                              const LocalObs & localobs) const {
  oops::Log::trace() << "ObsHorLocalization::computeLocalization(lengthscale)" << std::endl;

  // set all to missing (outside of localization distance)
  const double missing = util::missingValue<double>();
  for (size_t jj = 0; jj < locvector.size(); ++jj) {
    locvector[jj] = missing;
  }

  // set localization for the obs inside localization distance to 1.0
  const size_t nvars = locvector.nvars();
  const size_t nlocal = localobs.index.size();
  for (size_t jlocal = 0; jlocal < nlocal; ++jlocal) {
    // obsdist is calculated at each location; need to update R for each variable
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      locvector[jvar + localobs.index[jlocal] * nvars] = 1.0;
    }
  }
}


template<typename MODEL>
const typename ObsHorLocalization<MODEL>::LocalObs
ObsHorLocalization<MODEL>::getLocalObs(const GeometryIterator_ & i,
                                    double lengthscale) const {
  oops::Log::trace() << "ObsHorLocalization::getLocalObs" << std::endl;

  if ( lengthscale <= 0.0 ) {
    throw eckit::BadParameter("lengthscale parameter should be >= 0.0");
  }

  // check that this distribution supports local obs space
  // TODO(travis) this should be in the constructor, but currently
  //  breaks LETKF when using a split observer/solver
  if ( distName_ != "Halo" && distName_ != "InefficientDistribution" ) {
    std::string message = "Can not use ObsHorLocalization with distribution=" + distName_;
    throw eckit::BadParameter(message);
  }

  LocalObs localobs;
  localobs.lengthscale = lengthscale;
  eckit::geometry::Point3 refPoint = *i;
  eckit::geometry::Point2 refPoint2(refPoint[0], refPoint[1]);
  size_t nlocs = lons_.size();
  if ( options_.searchMethod == SearchMethod::BRUTEFORCE ) {
    oops::Log::trace() << "Local obs searching via brute force." << std::endl;

    for (unsigned int jj = 0; jj < nlocs; ++jj) {
      eckit::geometry::Point3 searchPoint(lons_[jj], lats_[jj], 0.0);
      double localDist = options_.distance(refPoint, searchPoint);
      if ( localDist < lengthscale ) {
        localobs.index.push_back(jj);
        localobs.distance.push_back(localDist);
      }
    }
    const boost::optional<int> & maxnobs = options_.maxnobs;
    if ( (maxnobs != boost::none) && (localobs.index.size() > *maxnobs ) ) {
      for (unsigned int jj = 0; jj < localobs.index.size(); ++jj) {
          oops::Log::debug() << "Before sort [i, d]: " << localobs.index[jj]
              << " , " << localobs.distance[jj] << std::endl;
      }
      // Construct a temporary paired vector to do the sorting
      std::vector<std::pair<std::size_t, double>> localObsIndDistPair;
      for (unsigned int jj = 0; jj < localobs.distance.size(); ++jj) {
        localObsIndDistPair.push_back(std::make_pair(localobs.index[jj], localobs.distance[jj]));
      }

      // Use a lambda function to implement an ascending sort.
      sort(localObsIndDistPair.begin(), localObsIndDistPair.end(),
           [](const std::pair<std::size_t, double> & p1,
              const std::pair<std::size_t, double> & p2){
                return(p1.second < p2.second);
              });

      // Unpair the sorted pair vector
      for (unsigned int jj = 0; jj < localobs.distance.size(); ++jj) {
        localobs.index[jj] = localObsIndDistPair[jj].first;
        localobs.distance[jj] = localObsIndDistPair[jj].second;
      }

      // Truncate to maxNobs length
      localobs.index.resize(*maxnobs);
      localobs.distance.resize(*maxnobs);
    }
  } else if (nlocs > 0) {
    // Check (nlocs > 0) is needed,
    // otherwise, it will cause ASERT check fail in kdtree.findInSphere, and hang.

    oops::Log::trace() << "Local obs searching via KDTree" << std::endl;

    if ( options_.distanceType == DistanceType::CARTESIAN)
      ABORT("ObsHorLocalization:: search method must be 'brute_force' when using"
            " 'cartesian' distance");

    // Using the radius of the earth
    eckit::geometry::Point3 refPoint3DTemp;
    atlas::util::Earth::convertSphericalToCartesian(refPoint2, refPoint3DTemp);
    double alpha =  (lengthscale / options_.radius_earth)/ 2.0;  // angle in radians
    double chordLength = 2.0*options_.radius_earth * sin(alpha);  // search radius in 3D space

    auto closePoints = kd_->findInSphere(refPoint3DTemp, chordLength);

    // put closePoints back into localobs and obsdist
    localobs.index.reserve(closePoints.size());
    localobs.distance.reserve(closePoints.size());
    for (unsigned int jloc = 0; jloc < closePoints.size(); ++jloc) {
       localobs.index.push_back(closePoints[jloc].payload());  // observation
       localobs.distance.push_back(closePoints[jloc].distance());  // distance
    }

    // The obs are sorted in the kdtree call
    const boost::optional<int> & maxnobs = options_.maxnobs;
    if ( (maxnobs != boost::none) && (localobs.index.size() > *maxnobs ) ) {
      // Truncate to maxNobs length
      localobs.index.resize(*maxnobs);
      localobs.distance.resize(*maxnobs);
    }
  }

  return localobs;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsHorLocalization<MODEL>::print(std::ostream & os) const {
  os << "ObsHorLocalization (box car) horizontal localization with " << options_.lengthscale
     << " lengthscale" << std::endl;
}

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSHORLOCALIZATION_H_
