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

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// Gaussian_Thinning: Thin observations to a gaussian grid

class Gaussian_Thinning : public util::Printable,
                  private util::ObjectCounter<Gaussian_Thinning> {
 public:
  static const std::string classname() {return "ufo::Gaussian_Thinning";}

  Gaussian_Thinning(ioda::ObsSpace &, const eckit::Configuration &,
           boost::shared_ptr<ioda::ObsDataVector<int> >,
           boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~Gaussian_Thinning();

  void preProcess() const;
  void priorFilter(const GeoVaLs &) const {}
  void postFilter(const ioda::ObsVector &, const ObsDiagnostics &) const {}

  const oops::Variables & requiredGeoVaLs() const {return geovars_;}
  const oops::Variables & requiredHdiagnostics() const {return diagvars_;}

 private:
  static int ll_to_idx(float in_lon, float in_lat, int n_lats, const std::vector<int> &n_lons);
  static int pres_to_idx(float in_pres, int n_vmesh, float vertical_mesh, float vertical_max);
  static int dist_to_centroid(float ob_lon, float ob_lat, float c_lon, float c_lat);
  void print(std::ostream &) const;

  ioda::ObsSpace & obsdb_;
  const eckit::LocalConfiguration config_;
  const oops::Variables geovars_;
  const oops::Variables diagvars_;
  ioda::ObsDataVector<int> & flags_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_GAUSSIAN_THINNING_H_
