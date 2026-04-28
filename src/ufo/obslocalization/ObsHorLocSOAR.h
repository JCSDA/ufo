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
template<class ITERATOR>
class ObsHorLocSOAR: public ufo::ObsHorLocalization<ITERATOR> {
  typedef typename ObsHorLocalization<ITERATOR>::LocalObs LocalObs_;

 public:
  ObsHorLocSOAR(const eckit::Configuration &, ioda::ObsSpace &);

  double computeLocalization(const eckit::geometry::Point3 &,
                             const eckit::geometry::Point3 &) const override;
 protected:
  /// Compute SOAR localization using the set of \p localobs and save localization
  /// values in \p locvector.
  /// (missing values indicate that observation is outside of localization)
  void localizeLocalObs(const ITERATOR &,
                        ioda::ObsVector & locvector,
                        const LocalObs_ & localobs) const override;

 private:
  void print(std::ostream &) const override;

  ObsHorLocSOARParameters options_;
};
// -----------------------------------------------------------------------------

template<typename ITERATOR>
ObsHorLocSOAR<ITERATOR>::ObsHorLocSOAR(const eckit::Configuration & config,
                                       ioda::ObsSpace & obsspace)
  : ObsHorLocalization<ITERATOR>::ObsHorLocalization(config, obsspace), options_()
{
  options_.validateAndDeserialize(config);
}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
double ObsHorLocSOAR<ITERATOR>::computeLocalization(
      const eckit::geometry::Point3 & p1,
      const eckit::geometry::Point3 & p2) const {
  eckit::geometry::Point2 p1_2(p1[0], p1[1]);
  eckit::geometry::Point2 p2_2(p2[0], p2[1]);
  double distance = atlas::util::Earth::distance(p1_2, p2_2);
  double lengthscale = ObsHorLocalization<ITERATOR>::lengthscale();
  double loc = oops::soar(distance * options_.SOARexpDecayH);
  return (distance >= lengthscale) ? 0.0 : loc;
}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
void ObsHorLocSOAR<ITERATOR>::localizeLocalObs(const ITERATOR & i,
                                        ioda::ObsVector & locvector,
                                        const LocalObs_ & localobs) const {
  // Apply box car localization
  ObsHorLocalization<ITERATOR>::localizeLocalObs(i, locvector, localobs);

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

template<typename ITERATOR>
void ObsHorLocSOAR<ITERATOR>::print(std::ostream & os) const {
  os << "SOAR horizontal localization with " << options_.SOARexpDecayH
     << " soar decay" << std::endl;
}

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSHORLOCSOAR_H_
