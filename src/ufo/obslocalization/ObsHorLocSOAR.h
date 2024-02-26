/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSHORLOCSOAR_H_
#define UFO_OBSLOCALIZATION_OBSHORLOCSOAR_H_

#include <ostream>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/generic/soar.h"

#include "ufo/obslocalization/ObsHorLocalization.h"
#include "ufo/obslocalization/ObsHorLocSOARParameters.h"

namespace ufo {

/// Horizontal SOAR observation space localization
template<class MODEL>
class ObsHorLocSOAR: public ufo::ObsHorLocalization<MODEL> {
  typedef typename MODEL::GeometryIterator   GeometryIterator_;
  typedef typename ObsHorLocalization<MODEL>::LocalObs LocalObs_;

 public:
  ObsHorLocSOAR(const eckit::Configuration &, const ioda::ObsSpace &);

 protected:
  /// Compute SOAR localization using the set of \p localobs and save localization
  /// values in \p locvector.
  /// (missing values indicate that observation is outside of localization)
  void localizeLocalObs(const GeometryIterator_ &,
                        ioda::ObsVector & locvector,
                        const LocalObs_ & localobs) const override;

 private:
  void print(std::ostream &) const override;

  ObsHorLocSOARParameters options_;
};
// -----------------------------------------------------------------------------

template<typename MODEL>
ObsHorLocSOAR<MODEL>::ObsHorLocSOAR(const eckit::Configuration & config,
                                    const ioda::ObsSpace & obsspace)
  : ObsHorLocalization<MODEL>::ObsHorLocalization(config, obsspace), options_()
{
  options_.validateAndDeserialize(config);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsHorLocSOAR<MODEL>::localizeLocalObs(const GeometryIterator_ & i,
                                        ioda::ObsVector & locvector,
                                        const LocalObs_ & localobs) const {
  // Apply box car localization
  ObsHorLocalization<MODEL>::localizeLocalObs(i, locvector, localobs);

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
void ObsHorLocSOAR<MODEL>::print(std::ostream & os) const {
  os << "SOAR horizontal localization with " << options_.SOARexpDecayH
     << " soar decay" << std::endl;
}

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSHORLOCSOAR_H_
