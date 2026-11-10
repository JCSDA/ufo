/*
 * (C) Copyright 2026 NCAR
 * (C) Copyright 2026 University of Oklahoma / Cooperative Institute for
 *                    Severe and High-Impact Weather Research and Operations (CIWRO)
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Authors:
 *   Jun Park (NCAR)
 *   Yunheng Wang (CIWRO / University of Oklahoma) -- contact
 */

#ifndef UFO_OBSLOCALIZATION_OBS3DLOCALIZATION_H_
#define UFO_OBSLOCALIZATION_OBS3DLOCALIZATION_H_

#include <algorithm>
#include <cfloat>
#include <cmath>
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

#include "ioda/distribution/Halo.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/generic/gc99.h"
#include "oops/util/missingValues.h"

#include "ufo/obslocalization/Obs3DLocParameters.h"
#include "ufo/obslocalization/ObsLocalizationBase.h"

namespace ufo {

/// \brief 3D observation-space localization using a combined normalized distance:
///   combined = sqrt((r_H / L_H)^2 + (r_V / L_V)^2)
/// where r_H is horizontal geodesic distance, r_V is vertical coordinate distance,
/// and L_H, L_V are the respective lengthscales.
/// The GC99 localization function is applied to `combined`.
/// Observations with combined >= 1 are excluded (set to missing).
///
/// \note Rationale. This is a true 3D ("ellipsoid") distance, in contrast to
/// applying separate horizontal and vertical GC99 tapers and multiplying them
/// (an "H x V" localization). The multiplicative form is order-independent only
/// in principle; in UFO the horizontal taper initialises the localization vector
/// while the vertical taper multiplies it, so only the "H x V" ordering is valid.
/// A combined 3D distance removes that ordering ambiguity and is expected to be
/// more important for convective-scale data assimilation, where small horizontal
/// and vertical localization lengthscales make the difference between the ellipsoid
/// and the product of two 1D tapers more significant. It also reduces the LETKF
/// measurement-update cost.
///
/// \note This observation-space localization is used by algorithms that perform
/// vertical localization in observation space (e.g. LETKF), which require a
/// vertical coordinate (e.g. `height`) in the IODA file. It is not needed for
/// GETKF, which performs model-space vertical localization; GETKF only applies
/// horizontal localization in the YAML and therefore does not require a vertical
/// coordinate in the IODA files.
template<class ITERATOR>
class Obs3DLocalization : public ObsLocalizationBase<ITERATOR> {
 public:
  Obs3DLocalization(const eckit::Configuration &, const ioda::ObsSpace &);

  void computeLocalization(const ITERATOR &,
                           ioda::ObsVector & locvector) const override;

  /// \brief Compute the 3D (ellipsoid) Gaspari-Cohn localization weight between
  /// two points. \p p1 and \p p2 carry longitude, latitude and the vertical
  /// coordinate in their x, y and z components respectively. Returns
  /// gc99(sqrt((r_H/L_H)^2 + (r_V/L_V)^2)), or 0 outside the unit ellipsoid.
  double computeLocalization(const eckit::geometry::Point3 &,
                             const eckit::geometry::Point3 &) const override;

 private:
  Obs3DLocParameters options_;

  /// Vertical coordinates of all observations
  std::vector<float> vCoord_;

  /// Horizontal coordinates of all observations
  std::vector<float> lats_;
  std::vector<float> lons_;

  std::string distName_;
  double haloSize_;

  /// Localization lengthscales, read and validated once in the constructor
  double horLengthscale_;
  double vertLengthscale_;

  void print(std::ostream &) const override;

  // KD-tree for horizontal pre-screening
  struct TreeTrait {
    typedef eckit::geometry::Point3 Point;
    typedef double                  Payload;
  };
  typedef eckit::KDTreeMemory<TreeTrait> KDTree;
  std::unique_ptr<KDTree> kd_;
};

// -----------------------------------------------------------------------------

template<typename ITERATOR>
Obs3DLocalization<ITERATOR>::Obs3DLocalization(const eckit::Configuration & config,
                                               const ioda::ObsSpace & obsspace)
  : options_(), lats_(obsspace.nlocs()), lons_(obsspace.nlocs())
{
  options_.validateAndDeserialize(config);

  distName_ = obsspace.distribution()->name();
  if (distName_ == "Halo") {
    auto haloDist = std::dynamic_pointer_cast<const ioda::Halo>(obsspace.distribution());
    haloSize_ = haloDist->haloSize();
  } else {
    haloSize_ = 1.e+20;
  }

  // Read and validate the localization lengthscales once here, rather than on
  // every grid point in computeLocalization().
  horLengthscale_  = options_.horLengthscale;
  vertLengthscale_ = options_.vertLengthscale;
  if (horLengthscale_  <= 0.0)
    throw eckit::BadParameter("horizontal lengthscale must be > 0");
  if (vertLengthscale_ <= 0.0)
    throw eckit::BadParameter("vertical lengthscale must be > 0");

  // The horizontal pre-screening is bounded by the halo; a horizontal
  // lengthscale larger than the halo would silently miss observations
  // (cf. ObsHorLocalization::getLocalObs).
  if (distName_ == "Halo" && horLengthscale_ > haloSize_) {
    throw eckit::BadParameter("Can not use Obs3DLocalization with distribution=" + distName_ +
                              " and horizontal lengthscale=" + std::to_string(horLengthscale_) +
                              " greater than haloSize=" + std::to_string(haloSize_));
  }

  // Get horizontal coordinates
  obsspace.get_db("MetaData", "longitude", lons_);
  obsspace.get_db("MetaData", "latitude", lats_);

  // Get vertical coordinates
  const size_t nlocs = obsspace.nlocs();
  if (options_.logTransform.value()) {
    vCoord_.resize(nlocs);
    obsspace.get_db(options_.iodaVerticalCoordinateGroup,
                    options_.iodaVerticalCoordinate, vCoord_);
    for (unsigned int jj = 0; jj < nlocs; ++jj) {
      if (vCoord_[jj] == 0) { vCoord_[jj] = FLT_EPSILON; }
      vCoord_[jj] = log(vCoord_[jj]);
    }
  } else {
    vCoord_.resize(nlocs);
    obsspace.get_db(options_.iodaVerticalCoordinateGroup,
                    options_.iodaVerticalCoordinate, vCoord_);
  }

  // Build KD-tree for horizontal pre-screening
  if (options_.searchMethod == SearchMethod::KDTREE) {
    kd_ = std::unique_ptr<KDTree>(new KDTree());
    typedef typename KDTree::PointType Point;
    std::vector<typename KDTree::Value> points;
    points.reserve(nlocs);
    for (unsigned int i = 0; i < nlocs; ++i) {
      eckit::geometry::Point2 lonlat(lons_[i], lats_[i]);
      Point xyz = Point();
      atlas::util::Earth::convertSphericalToCartesian(lonlat, xyz);
      points.push_back(typename KDTree::Value(xyz, static_cast<double>(i)));
    }
    kd_->build(points.begin(), points.end());
  }
}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
void Obs3DLocalization<ITERATOR>::computeLocalization(const ITERATOR & i,
                                                      ioda::ObsVector & locvector) const {
  oops::Log::trace() << "Obs3DLocalization::computeLocalization" << std::endl;

  if ( distName_ != "Halo" && distName_ != "InefficientDistribution" ) {
    throw eckit::BadParameter("Can not use Obs3DLocalization with distribution=" + distName_);
  }

  const double horLS  = horLengthscale_;
  const double vertLS = vertLengthscale_;

  const double missing = util::missingValue<double>();
  const size_t nvars   = locvector.nvars();
  const size_t nlocs   = lons_.size();

  // Set all to missing initially
  for (size_t jj = 0; jj < locvector.size(); ++jj) {
    locvector[jj] = missing;
  }

  // Reference point
  eckit::geometry::Point3 refPoint3D = *i;
  eckit::geometry::Point2 refPoint2D(refPoint3D[0], refPoint3D[1]);
  double vCoordRef = refPoint3D[2];
  if (options_.logTransform.value()) {
    if (vCoordRef == 0) { vCoordRef = FLT_EPSILON; }
    vCoordRef = log(vCoordRef);
  }

  // --- Step 1: gather horizontal candidate indices ---
  std::vector<size_t> candidates;

  if (options_.searchMethod == SearchMethod::KDTREE && nlocs > 0) {
    // Use full horizontal lengthscale as KDTree search radius.
    // The ellipsoid boundary along the horizontal axis is exactly horLS.
    eckit::geometry::Point3 refPoint3DCart;
    atlas::util::Earth::convertSphericalToCartesian(refPoint2D, refPoint3DCart);
    double alpha       = (horLS / options_.radius_earth) / 2.0;
    double chordLength = 2.0 * options_.radius_earth * sin(alpha);

    auto closePoints = kd_->findInSphere(refPoint3DCart, chordLength);
    candidates.reserve(closePoints.size());
    for (const auto & cp : closePoints) {
      candidates.push_back(static_cast<size_t>(cp.payload()));
    }
  } else {
    // Brute force: all obs are candidates; 3D filter below will cull them
    candidates.resize(nlocs);
    for (size_t jj = 0; jj < nlocs; ++jj) candidates[jj] = jj;
  }

  // --- Step 2: compute combined 3D normalized distance and apply localization ---
  // Collect (index, combined_distance) pairs for sorting if maxnobs is set
  std::vector<std::pair<size_t, double>> localObsPairs;
  localObsPairs.reserve(candidates.size());

  for (size_t jj : candidates) {
    // Horizontal geodesic distance
    eckit::geometry::Point3 obsPoint(lons_[jj], lats_[jj], 0.0);
    eckit::geometry::Point2 obsPoint2(lons_[jj], lats_[jj]);

    // Convert chord back to geodesic if KDTree was used
    double horDist;
    if (options_.searchMethod == SearchMethod::KDTREE) {
      eckit::geometry::Point3 refCart, obsCart;
      atlas::util::Earth::convertSphericalToCartesian(refPoint2D, refCart);
      atlas::util::Earth::convertSphericalToCartesian(obsPoint2, obsCart);
      double chord = refCart.distance(obsCart);
      double sinHalfAlpha = std::min(chord / (2.0 * options_.radius_earth), 1.0);
      horDist = 2.0 * options_.radius_earth * std::asin(sinHalfAlpha);
    } else {
      // Brute force: use great circle distance directly
      horDist = eckit::geometry::Sphere::distance(options_.radius_earth,
                                                  refPoint2D, obsPoint2);
    }

    // Vertical distance
    double vertDist = std::abs(vCoordRef - static_cast<double>(vCoord_[jj]));

    // Combined normalized distance
    double L_H = horDist  / horLS;
    double L_V = vertDist / vertLS;
    double combined = std::sqrt(L_H * L_H + L_V * L_V);

    // Only include obs inside the unit ellipsoid (combined < 1)
    if (combined < 1.0) {
      localObsPairs.push_back({jj, combined});
    }
  }

  // Apply maxnobs: keep closest (smallest combined distance) obs
  const boost::optional<int> & maxnobs = options_.maxnobs;
  if ((maxnobs != boost::none) &&
      (static_cast<int>(localObsPairs.size()) > *maxnobs)) {
    std::sort(localObsPairs.begin(), localObsPairs.end(),
              [](const std::pair<size_t, double> & a,
                 const std::pair<size_t, double> & b) {
                return a.second < b.second;
              });
    localObsPairs.resize(*maxnobs);
  }

  // --- Step 3: apply Gaspari-Cohn to combined distance ---
  for (const auto & pair : localObsPairs) {
    const size_t idx       = pair.first;
    const double combined  = pair.second;
    const double locFactor = oops::gc99(combined);
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      locvector[jvar + idx * nvars] = locFactor;
    }
  }
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

template<typename ITERATOR>
double Obs3DLocalization<ITERATOR>::computeLocalization(
      const eckit::geometry::Point3 & p1,
      const eckit::geometry::Point3 & p2) const {
  oops::Log::trace() << "Obs3DLocalization::computeLocalization(Point3)" << std::endl;

  const double horLS  = horLengthscale_;
  const double vertLS = vertLengthscale_;

  // Horizontal geodesic distance (meters). Consistent with ObsHorLocalization.
  eckit::geometry::Point2 p1_2(p1[0], p1[1]);
  eckit::geometry::Point2 p2_2(p2[0], p2[1]);
  const double horDist = atlas::util::Earth::distance(p1_2, p2_2);

  // Vertical distance from the z-component, with optional log transform.
  // Consistent with ObsVertLocalization.
  double vCoord1 = p1[2];
  double vCoord2 = p2[2];
  if (options_.logTransform.value()) {
    if (vCoord1 == 0) { vCoord1 = FLT_EPSILON; }
    if (vCoord2 == 0) { vCoord2 = FLT_EPSILON; }
    vCoord1 = std::log(vCoord1);
    vCoord2 = std::log(vCoord2);
  }
  const double vertDist = std::abs(vCoord1 - vCoord2);

  // Combined normalized (ellipsoid) distance.
  const double L_H = horDist  / horLS;
  const double L_V = vertDist / vertLS;
  const double combined = std::sqrt(L_H * L_H + L_V * L_V);

  // Outside the unit ellipsoid -> zero; otherwise Gaspari-Cohn taper.
  if (combined >= 1.0) return 0.0;
  return oops::gc99(combined);
}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
void Obs3DLocalization<ITERATOR>::print(std::ostream & os) const {
  os << "Obs3DLocalization: Gaspari-Cohn"
     << " with horizontal lengthscale=" << options_.horLengthscale
     << " and vertical lengthscale=" << options_.vertLengthscale << std::endl;
}

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBS3DLOCALIZATION_H_
