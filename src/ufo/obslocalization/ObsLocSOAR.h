/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSLOCSOAR_H_
#define UFO_OBSLOCALIZATION_OBSLOCSOAR_H_

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

#include "oops/generic/soar.h"
#include "oops/util/Printable.h"

#include "ufo/obslocalization/ObsLocalization.h"
#include "ufo/obslocalization/ObsLocParameters.h"

namespace ufo {

/// Horizontal SOAR observation space localization
template<class MODEL>
class ObsLocSOAR: public ufo::ObsLocalization<MODEL> {
  typedef typename MODEL::GeometryIterator   GeometryIterator_;

 public:
  ObsLocSOAR(const eckit::Configuration &, const ioda::ObsSpace &);

  /// compute localization and save localization values in \p obsvector and
  /// localization flags (1: outside of localization; 0: inside localization area)
  /// in \p outside
  void computeLocalization(const GeometryIterator_ &, ioda::ObsDataVector<int> & outside,
                           ioda::ObsVector & obsvector) const override;

 private:
  void print(std::ostream &) const;
};
// -----------------------------------------------------------------------------

template<typename MODEL>
ObsLocSOAR<MODEL>::ObsLocSOAR(const eckit::Configuration & config,
                              const ioda::ObsSpace & obsspace):
       ObsLocalization<MODEL>::ObsLocalization(config, obsspace) {
  const ObsLocParameters & options = ObsLocalization<MODEL>::localizationOptions();
  oops::Log::debug()<< "SOAR horizontal localization with " << options.SOARexpDecayH
     << " soar decay" << std::endl;
}
// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocSOAR<MODEL>::computeLocalization(const GeometryIterator_ & i,
                                            ioda::ObsDataVector<int> & outside,
                                            ioda::ObsVector & locvector) const {
  oops::Log::trace() << "ObsLocSOAR::computeLocalization" << std::endl;
  // do distance search and compute box-car locvector
  ObsLocalization<MODEL>::computeLocalization(i, outside, locvector);

  // return refs to internals of ObsLocalization
  const std::vector<int> & localobs = ObsLocalization<MODEL>::localobs();
  const std::vector<double> & horizontalObsdist = ObsLocalization<MODEL>::horizontalObsdist();
  const ObsLocParameters  & options = ObsLocalization<MODEL>::localizationOptions();

  const double SOARexpDecayH = options.SOARexpDecayH;
  const size_t nvars = locvector.nvars();
  for (size_t jlocal = 0; jlocal < localobs.size(); ++jlocal) {
    double locFactor = oops::soar(horizontalObsdist[jlocal]*SOARexpDecayH);
    // obsdist is calculated at each location; need to update R for each variable
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      locvector[jvar + localobs[jlocal] * nvars] *= locFactor;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocSOAR<MODEL>::print(std::ostream & os) const {
  const ObsLocParameters & options = ObsLocalization<MODEL>::localizationOptions();
  os << "SOAR horizontal localization with " << options.SOARexpDecayH
     << " soar decay" << std::endl;
}

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSLOCSOAR_H_
