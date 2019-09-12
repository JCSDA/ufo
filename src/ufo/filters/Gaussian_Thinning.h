/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_GAUSSIAN_THINNING_H_
#define UFO_FILTERS_GAUSSIAN_THINNING_H_

#include <ostream>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "ioda/ObsDataVector.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Gaussian_Thinning: Thin observations to a gaussian grid

class Gaussian_Thinning : public FilterBase,
                          private util::ObjectCounter<Gaussian_Thinning> {
 public:
  static const std::string classname() {return "ufo::Gaussian_Thinning";}

  Gaussian_Thinning(ioda::ObsSpace &, const eckit::Configuration &,
           boost::shared_ptr<ioda::ObsDataVector<int> >,
           boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~Gaussian_Thinning();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, std::vector<std::vector<bool>> &) const override;

  static int ll_to_idx(float in_lon, float in_lat, int n_lats, const std::vector<int> &n_lons);
  static int pres_to_idx(float in_pres, int n_vmesh, float vertical_mesh, float vertical_max);
  static int dist_to_centroid(float ob_lon, float ob_lat, float c_lon, float c_lat);
};

}  // namespace ufo

#endif  // UFO_FILTERS_GAUSSIAN_THINNING_H_
