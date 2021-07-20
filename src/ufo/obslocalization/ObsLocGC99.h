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

 public:
  ObsLocGC99(const eckit::Configuration &, const ioda::ObsSpace &);

  /// compute localization and save localization values in \p locfactor
  /// (missing values indicate that observation is outside of localization)
  void computeLocalization(const GeometryIterator_ &,
                           ioda::ObsVector & locfactor) const override;

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

template<typename MODEL>
ObsLocGC99<MODEL>::ObsLocGC99(const eckit::Configuration & config,
                              const ioda::ObsSpace & obsspace):
       ObsLocalization<MODEL>::ObsLocalization(config, obsspace) {
  const ObsLocParameters & options = ObsLocalization<MODEL>::localizationOptions();
  oops::Log::debug() <<  "Gaspari-Cohn horizontal localization with " << options.lengthscale
     << " lengthscale" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocGC99<MODEL>::computeLocalization(const GeometryIterator_ & i,
                                            ioda::ObsVector & locvector) const {
  oops::Log::trace() << "ObsLocGC99::computeLocalization" << std::endl;

  // do distance search and compute box-car locvector
  ObsLocalization<MODEL>::computeLocalization(i, locvector);

  // return refs to internals of ObsLocalization
  const std::vector<int> & localobs = ObsLocalization<MODEL>::localobs();
  const std::vector<double> & horizontalObsdist = ObsLocalization<MODEL>::horizontalObsdist();
  const ObsLocParameters  & options = ObsLocalization<MODEL>::localizationOptions();

  const size_t nvars = locvector.nvars();
  for (size_t jlocal = 0; jlocal < localobs.size(); ++jlocal) {
    double locFactor = oops::gc99(horizontalObsdist[jlocal] / options.lengthscale);
    // obsdist is calculated at each location; need to update R for each variable
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      locvector[jvar + localobs[jlocal] * nvars] *= locFactor;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocGC99<MODEL>::print(std::ostream & os) const {
  const ObsLocParameters  & options = ObsLocalization<MODEL>::localizationOptions();
  os << "Gaspari-Cohn horizontal localization with " << options.lengthscale
     << " lengthscale" << std::endl;
}

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSLOCGC99_H_
