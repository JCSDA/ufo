/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_GNSSRO_QC_BACKGROUNDCHECKROGSI_H_
#define UFO_GNSSRO_QC_BACKGROUNDCHECKROGSI_H_

#include <ostream>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

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

/// BackgroundCheckROGSI: check observation closeness to background

class BackgroundCheckROGSI : public FilterBase,
                             private util::ObjectCounter<BackgroundCheckROGSI> {
 public:
  static const std::string classname() {return "ufo::BackgroundCheckROGSI";}

  BackgroundCheckROGSI(ioda::ObsSpace &, const eckit::Configuration &,
                       boost::shared_ptr<ioda::ObsDataVector<int> >,
                       boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~BackgroundCheckROGSI();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, std::vector<std::vector<bool>> &) const override;
};

}  // namespace ufo

#endif  // UFO_GNSSRO_QC_BACKGROUNDCHECKROGSI_H_
