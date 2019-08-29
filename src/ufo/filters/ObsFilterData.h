/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_OBSFILTERDATA_H_
#define UFO_FILTERS_OBSFILTERDATA_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

/// ObsFilterData: check observation closeness to background

class ObsFilterData : public util::Printable,
                      private util::ObjectCounter<ObsFilterData> {
 public:
  static const std::string classname() {return "ufo::ObsFilterData";}

  explicit ObsFilterData(ioda::ObsSpace &);
  ~ObsFilterData();

  void associate(const GeoVaLs &);
  void associate(const ioda::ObsVector &);

  std::vector<float> get(const std::string &) const;
  bool has(const std::string &) const;

  size_t nlocs() const;

 private:
  void print(std::ostream &) const;

  ioda::ObsSpace & obsdb_;
  const GeoVaLs mutable * gvals_;
  const ioda::ObsVector mutable * hofx_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFILTERDATA_H_
