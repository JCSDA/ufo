/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSLOCSOAR_H_
#define UFO_OBSLOCALIZATION_OBSLOCSOAR_H_

#include <ostream>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/generic/soar.h"

#include "ufo/obslocalization/ObsLocalization.h"
#include "ufo/obslocalization/ObsLocSOARParameters.h"

namespace ufo {

/// Horizontal SOAR observation space localization
template<class MODEL>
class ObsLocSOAR: public ufo::ObsLocalization<MODEL> {
  typedef typename MODEL::GeometryIterator   GeometryIterator_;
  typedef typename ObsLocalization<MODEL>::LocalObs LocalObs_;

 public:
  ObsLocSOAR(const eckit::Configuration &, const ioda::ObsSpace &);

 protected:
  /// Compute SOAR localization using the set of \p localobs and save localization
  /// values in \p locvector.
  /// (missing values indicate that observation is outside of localization)
  void localizeLocalObs(const GeometryIterator_ &,
                        ioda::ObsVector & locvector,
                        const LocalObs_ & localobs) const override;

 private:
  void print(std::ostream &) const override;

  ObsLocSOARParameters options_;
};
// -----------------------------------------------------------------------------

template<typename MODEL>
ObsLocSOAR<MODEL>::ObsLocSOAR(const eckit::Configuration & config,
                              const ioda::ObsSpace & obsspace):
       ObsLocalization<MODEL>::ObsLocalization(config, obsspace) {
  options_.deserialize(config);
  oops::Log::debug()<< "SOAR horizontal localization with " << options_.SOARexpDecayH
     << " soar decay" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocSOAR<MODEL>::localizeLocalObs(const GeometryIterator_ & i,
                                        ioda::ObsVector & locvector,
                                        const LocalObs_ & localobs) const {
  // Apply box car localization
  ObsLocalization<MODEL>::localizeLocalObs(i, locvector, localobs);

  // Apply SOAR localization
  const double SOARexpDecayH = options_.SOARexpDecayH;
  const size_t nvars = locvector.nvars();
  for (size_t jlocal = 0; jlocal < localobs.index.size(); ++jlocal) {
    double locFactor = oops::soar(localobs.distance[jlocal]*SOARexpDecayH);
    // obsdist is calculated at each location; need to update R for each variable
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      locvector[jvar + localobs.index[jlocal] * nvars] = locFactor;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocSOAR<MODEL>::print(std::ostream & os) const {
  os << "SOAR horizontal localization with " << options_.SOARexpDecayH
     << " soar decay" << std::endl;
}

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSLOCSOAR_H_
