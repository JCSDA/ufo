/*
 * CC0 - No Copyright (Public Domain)
 * 
 * The person who associated a work with this deed has dedicated the work to the public domain by waiving all of his or her rights to the work worldwide under copyright law, including all related and neighboring rights, to the extent allowed by law.
 *
 * You can copy, modify, distribute and perform the work, even for commercial purposes, all without asking permission. See Other Information below.
 *
 * Other Information
 * In no way are the patent or trademark rights of any person affected by CC0, nor are the rights that other persons may have in the work or in how the work is used, such as publicity or privacy rights.
 *
 * The person who associated a work with this deed makes no warranties about the work, and disclaims liability for all uses of the work, to the fullest extent permitted by applicable law.
 *
 * When using or citing the work, you should not imply endorsement by the author or the affirmer.
 */

#ifndef UFO_FILTERS_OBSPOLYGONCHECK_H_
#define UFO_FILTERS_OBSPOLYGONCHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

class ObsPolygonCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsPolygonCheckParameters, FilterParametersBase)
  oops::RequiredParameter <std::string> polyFileName
    {"polygon wkt file",
     "Path to a Well Known Text (WKT) file containing the polygon as longitude, latitude pairs. "
     "Observations that lie within this polygon are retained. "
     "File contents should look like: \"POLYGON((180 -3, -125 18, 98.4 0, 180 -3))\""
     this};
  oops::RequiredParameter <float> insideLon
    {"inside point longitude",
     "Longitude of the sample point within the domain, used to find which side of the polygon is inside.",
     this};
  oops::RequiredParameter <float> insideLat
    {"inside point latitude",
     "Latitude of the sample point within the domain, used to find which side of the polygon is inside.",
     this};
};

/// PolygonCheck: find obs within a specified polygon.

class ObsPolygonCheck : public FilterBase,
                       private util::ObjectCounter<ObsPolygonCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ObsPolygonCheckParameters Parameters_;

  static const std::string classname() {return "ufo::ObsPolygonCheck";}

  ObsPolygonCheck(ioda::ObsSpace &, const Parameters_ &,
                 std::shared_ptr<ioda::ObsDataVector<int> >,
                 std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsPolygonCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::domain;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSPOLYGONCHECK_H_
