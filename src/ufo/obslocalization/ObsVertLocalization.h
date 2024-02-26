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

#include "oops/base/ObsLocalizationBase.h"
#include "oops/util/missingValues.h"

#include "ufo/obslocalization/ObsVertLocParameters.h"
#include "ufo/ObsTraits.h"

namespace ufo {

/// Vertical observation space localization
template<class MODEL>
class ObsVertLocalization: public oops::ObsLocalizationBase<MODEL, ObsTraits> {
  typedef typename MODEL::GeometryIterator   GeometryIterator_;

 public:
  ObsVertLocalization(const eckit::Configuration &, const ioda::ObsSpace &);

  /// Compute localization and save localization values in \p locvector.
  /// Missing values indicate that observation is outside of localization.
  /// The lengthscale from ObsVertLocParameters is used.
  void computeLocalization(const GeometryIterator_ &,
                           ioda::ObsVector & locvector) const override;

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
  const LocalObs getLocalObs(const GeometryIterator_ &, double lengthscale) const;

  /// Compute localization using the set of \p localobs and save
  /// localization values in \p locvector. Missing values are set for obs
  /// outside of localization.
  /// Intended to be called by \c computeLocalization() .
  void localizeLocalObs(const GeometryIterator_ &,
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

template<typename MODEL>
ObsVertLocalization<MODEL>::ObsVertLocalization(const eckit::Configuration & config,
                                                const ioda::ObsSpace & obsspace)
  : options_()
{
  options_.validateAndDeserialize(config);
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
      vCoord_[jj] = log(vCoord_[jj]);
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsVertLocalization<MODEL>::computeLocalization(const GeometryIterator_ & i,
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

template<typename MODEL>
void ObsVertLocalization<MODEL>::localizeLocalObs(const GeometryIterator_ & i,
                                              ioda::ObsVector & locvector,
                                              const LocalObs & localobs) const {
  oops::Log::trace() << "ObsVertLocalization::computeLocalization(lengthscale)" << std::endl;

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
      if (SOARexpDecayH == util::missingValue<double>()) {
        std::string message = "soar horizontal decay parameter is not specified";
        throw eckit::BadParameter(message);
      }
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

template<typename MODEL>
const typename ObsVertLocalization<MODEL>::LocalObs
ObsVertLocalization<MODEL>::getLocalObs(const GeometryIterator_ & i,
                                    double lengthscale) const {
  oops::Log::trace() << "ObsVertLocalization::getLocalObs" << std::endl;

  if ( lengthscale <= 0.0 ) {
    throw eckit::BadParameter("lengthscale parameter should be >= 0.0");
  }

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
    vCoordAtIterator = log(vCoordAtIterator);
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

  return localobs;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsVertLocalization<MODEL>::print(std::ostream & os) const {
  os << "ObsVertLocalization with " << options_.lengthscale
     << " lengthscale and " << options_.localizationFunction.value()
     << " localization function"<< std::endl;
}

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSVERTLOCALIZATION_H_
