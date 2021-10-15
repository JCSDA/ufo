/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSLOCGC99_H_
#define UFO_OBSLOCALIZATION_OBSLOCGC99_H_

#include <ostream>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/generic/gc99.h"

#include "ufo/obslocalization/ObsLocalization.h"
#include "ufo/obslocalization/ObsLocParameters.h"

namespace ufo {

/// Horizontal Gaspari-Cohn observation space localization
template<class MODEL>
class ObsLocGC99: public ufo::ObsLocalization<MODEL> {
  typedef typename MODEL::GeometryIterator   GeometryIterator_;
  typedef typename ObsLocalization<MODEL>::LocalObs LocalObs_;

 public:
  ObsLocGC99(const eckit::Configuration &, const ioda::ObsSpace &);

 protected:
  /// Compute GC99 localization using the set of \p localobs, the same lengthscale
  /// used to search for those obs, and save localization values in \p locvector.
  /// (missing values indicate that observation is outside of localization)
  void localizeLocalObs(const GeometryIterator_ &,
                        ioda::ObsVector & locvector,
                        const LocalObs_ &) const override;

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

template<typename MODEL>
ObsLocGC99<MODEL>::ObsLocGC99(const eckit::Configuration & config,
                              const ioda::ObsSpace & obsspace):
       ObsLocalization<MODEL>::ObsLocalization(config, obsspace) {
  oops::Log::debug() <<  *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocGC99<MODEL>::localizeLocalObs(const GeometryIterator_ & i,
                                        ioda::ObsVector & locvector,
                                        const LocalObs_ & localobs) const {
  oops::Log::trace() << "ObsLocGC99::computeLocalization" << std::endl;

  // Apply box car localization
  ObsLocalization<MODEL>::localizeLocalObs(i, locvector, localobs);

  // Apply Gaspari-Cohn localization
  const size_t nvars = locvector.nvars();
  const size_t nlocal = localobs.index.size();
  for (size_t jlocal = 0; jlocal < nlocal; ++jlocal) {
    double locFactor = oops::gc99(localobs.distance[jlocal] / localobs.lengthscale);

    // obsdist is calculated at each location; need to update R for each variable
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      locvector[jvar + localobs.index[jlocal] * nvars] *= locFactor;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocGC99<MODEL>::print(std::ostream & os) const {
  const double lengthscale = ObsLocalization<MODEL>::lengthscale();

  os << "Gaspari-Cohn horizontal localization with " << lengthscale
     << " lengthscale" << std::endl;
}

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSLOCGC99_H_
