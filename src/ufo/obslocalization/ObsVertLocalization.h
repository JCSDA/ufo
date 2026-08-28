/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSVERTLOCALIZATION_H_
#define UFO_OBSLOCALIZATION_OBSVERTLOCALIZATION_H_

#include <algorithm>
#include <cfloat>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point3.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/generic/gc99.h"
#include "oops/generic/soar.h"
#include "oops/util/missingValues.h"

#include "ufo/obslocalization/ObsLocalizationBase.h"
#include "ufo/obslocalization/ObsVertLocParameters.h"

namespace ufo {

/// Vertical observation space localization
template<class ITERATOR>
class ObsVertLocalization: public ObsLocalizationBase<ITERATOR> {
 public:
  ObsVertLocalization(const eckit::Configuration &, const ioda::ObsSpace &);

  /// Compute localization and save localization values in \p locvector.
  /// Missing values indicate that observation is outside of localization.
  /// The lengthscale from ObsVertLocParameters is used.
  void computeLocalization(const ITERATOR &,
                           ioda::ObsVector & locvector) const override;

  /// Compute localization between two points using the z-coordinate as the
  /// vertical coordinate. Returns a value between 0.0 and 1.0.
  double computeLocalization(const eckit::geometry::Point3 & p1,
                             const eckit::geometry::Point3 & p2) const override;

 protected:
  struct LocalObs {
    /// The list of indexes for ObsVector pointing to the valid local obs.
    std::vector<int> index;

    /// The vertical distance of each local ob from the search point.
    std::vector<double> distance;

    /// The maximum search distance that was used for this local obs search.
    double lengthscale;
  };

  /// For a given distance, returns the local observations and their distances.
  /// Intended to be called by \c computeLocalization() .
  const LocalObs getLocalObs(const ITERATOR &, double lengthscale) const;

  /// Compute localization using the set of \p localobs and save
  /// localization values in \p locvector. Missing values are set for obs
  /// outside of localization.
  /// Intended to be called by \c computeLocalization() .
  void localizeLocalObs(const ITERATOR &,
                                ioda::ObsVector & locvector,
                                const LocalObs & localobs) const;

  /// Get the lengthscale specified in the parameters.
  double lengthscale() const {return options_.lengthscale;}

 private:
  ObsVertLocParameters options_;
  std::vector<float> vCoord_;

  void print(std::ostream &) const override;

  /// TODO(travis) distribution name is needed for temporary fix, should be removed eventually
  std::string distName_;
};

// -----------------------------------------------------------------------------

template<typename ITERATOR>
ObsVertLocalization<ITERATOR>::ObsVertLocalization(const eckit::Configuration & config,
                                                const ioda::ObsSpace & obsspace)
  : options_()
{
  options_.validateAndDeserialize(config);

  // Reject invalid parameters up front so both compute paths
  // (Point3/Point3 and iterator-based) can rely on valid values.
  if (options_.lengthscale <= 0.0) {
    throw eckit::BadParameter("vertical lengthscale parameter should be > 0.0");
  }
  if (options_.localizationFunction.value().compare("SOAR") == 0 &&
      options_.SOARexpDecayH == util::missingValue<double>()) {
    throw eckit::BadParameter("soar decay parameter is not specified");
  }

  // check that this distribution supports local obs space
  // TODO(travis) this has been moved to computeLocalization as a quick fix for a bug.
  distName_ = obsspace.distribution()->name();

  // Get vertical coordinate of all observations.
  if (options_.assignConstantVcoordToObs.value()) {
    vCoord_.resize(obsspace.nlocs(), options_.constantVcoordValue.value() );
  } else {
    obsspace.get_db(options_.iodaVerticalCoordinateGroup, options_.iodaVerticalCoordinate, vCoord_);
  }
  if (options_.logTransform.value()) {
    for (unsigned int jj = 0; jj < vCoord_.size(); ++jj) {
      if (vCoord_[jj] == 0) { vCoord_[jj] = FLT_EPSILON; }
      vCoord_[jj] = std::log(vCoord_[jj]);
    }
  }
}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
void ObsVertLocalization<ITERATOR>::computeLocalization(const ITERATOR & i,
                                                 ioda::ObsVector & locvector) const {
  oops::Log::trace() << "ObsVertLocalization::computeLocalization" << std::endl;

  // get the set of local observations using the lengthscale given in the
  // config file options.
  const LocalObs & localobs = getLocalObs(i, options_.lengthscale);

  // compute localization of those local obs. Note that since this is
  // a virtual method, it could be overriden by dervied classes
  localizeLocalObs(i, locvector, localobs);
}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
double ObsVertLocalization<ITERATOR>::computeLocalization(
      const eckit::geometry::Point3 & p1,
      const eckit::geometry::Point3 & p2) const {
  oops::Log::trace() << "ObsVertLocalization::computeLocalization(Point3)" << std::endl;

  // Extract vertical coordinates from Point3 z-component
  double vCoord1 = p1[2];
  double vCoord2 = p2[2];

  // Apply log transform if configured
  if (options_.logTransform.value()) {
    if (vCoord1 == 0) { vCoord1 = FLT_EPSILON; }
    if (vCoord2 == 0) { vCoord2 = FLT_EPSILON; }
    vCoord1 = std::log(vCoord1);
    vCoord2 = std::log(vCoord2);
  }

  double distance = options_.distance(vCoord1, vCoord2);
  double ls = options_.lengthscale;

  if (distance >= ls) {
    return 0.0;
  }

  if (options_.localizationFunction.value().compare("Box Car") == 0) {
    return 1.0;
  } else if (options_.localizationFunction.value().compare("Gaspari Cohn") == 0) {
    return oops::gc99(distance / ls);
  } else if (options_.localizationFunction.value().compare("SOAR") == 0) {
    return oops::soar(distance * options_.SOARexpDecayH);
  } else {
    throw eckit::BadParameter("Vertical correlation function not recognized "
                              + options_.localizationFunction.value());
  }
}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
void ObsVertLocalization<ITERATOR>::localizeLocalObs(const ITERATOR & i,
                                              ioda::ObsVector & locvector,
                                              const LocalObs & localobs) const {
  oops::Log::trace() << "ObsVertLocalization::localizeLocalObs" << std::endl;

  const double missing = util::missingValue<double>();
  const size_t nvars = locvector.nvars();
  const size_t nlocal = localobs.index.size();

  ioda::ObsVector locvectorTmp(locvector);
  for (size_t jj = 0; jj < locvector.size(); ++jj) {
    locvector[jj] = missing;
  }

  if (options_.localizationFunction.value().compare("Box Car") == 0) {
      for (size_t jlocal = 0; jlocal < nlocal; ++jlocal) {
        for (size_t jvar = 0; jvar < nvars; ++jvar) {
          locvector[jvar + localobs.index[jlocal] * nvars] =
            locvectorTmp[jvar + localobs.index[jlocal] * nvars];
        }
      }
  } else if (options_.localizationFunction.value().compare("Gaspari Cohn") == 0) {
      for (size_t jlocal = 0; jlocal < nlocal; ++jlocal) {
        double locFactor = oops::gc99(localobs.distance[jlocal] / localobs.lengthscale);
        for (size_t jvar = 0; jvar < nvars; ++jvar) {
          locvector[jvar + localobs.index[jlocal] * nvars] = locFactor*
            locvectorTmp[jvar + localobs.index[jlocal] * nvars];
        }
      }
  } else if (options_.localizationFunction.value().compare("SOAR") == 0) {
      const double SOARexpDecayH = options_.SOARexpDecayH;
      for (size_t jlocal = 0; jlocal < nlocal; ++jlocal) {
        double locFactor = oops::soar(localobs.distance[jlocal]*SOARexpDecayH);
        for (size_t jvar = 0; jvar < nvars; ++jvar) {
          locvector[jvar + localobs.index[jlocal] * nvars] = locFactor*
            locvectorTmp[jvar + localobs.index[jlocal] * nvars];
        }
      }
  } else {
      std::string message = "Vertical correlation function not recognized "
                            +options_.localizationFunction.value();
      throw eckit::BadParameter(message);
  }

  // make sure that locvector has the same missing value as on input
  for (size_t jj = 0; jj < locvector.size(); ++jj) {
    if (locvectorTmp[jj] == missing) {locvector[jj] = missing;}
  }
}

template<typename ITERATOR>
const typename ObsVertLocalization<ITERATOR>::LocalObs
ObsVertLocalization<ITERATOR>::getLocalObs(const ITERATOR & i,
                                    double lengthscale) const {
  oops::Log::trace() << "ObsVertLocalization::getLocalObs" << std::endl;

  // check that this distribution supports local obs space
  // TODO(travis) this should be in the constructor, but currently
  //  breaks LETKF when using a split observer/solver
  if ( distName_ != "Halo" && distName_ != "InefficientDistribution" ) {
    std::string message = "Can not use ObsVertLocalization with distribution=" + distName_;
    throw eckit::BadParameter(message);
  }

  LocalObs localobs;
  localobs.lengthscale = lengthscale;

  eckit::geometry::Point3 refPoint = *i;
  double vCoordAtIterator = refPoint[2];
  if (options_.logTransform.value()) {
    if (vCoordAtIterator == 0) { vCoordAtIterator = FLT_EPSILON; }
    vCoordAtIterator = std::log(vCoordAtIterator);
  }
  size_t nlocs = vCoord_.size();
  for (unsigned int jj = 0; jj < nlocs; ++jj) {
    double localDist = options_.distance(vCoordAtIterator, vCoord_[jj]);
    if ( localDist < lengthscale ) {
      localobs.index.push_back(jj);
      localobs.distance.push_back(localDist);
    }
  }
  // truncate to maxnobs if needed
  const boost::optional<int> & maxnobs = options_.maxnobs;
  if ( (maxnobs != boost::none) && (static_cast<int>(localobs.index.size()) > *maxnobs ) ) {
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

  return localobs;
}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
void ObsVertLocalization<ITERATOR>::print(std::ostream & os) const {
  os << "ObsVertLocalization with " << options_.lengthscale
     << " lengthscale and " << options_.localizationFunction.value()
     << " localization function"<< std::endl;
}

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSVERTLOCALIZATION_H_
