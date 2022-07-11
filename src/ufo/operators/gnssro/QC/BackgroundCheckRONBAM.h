/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GNSSRO_QC_BACKGROUNDCHECKRONBAM_H_
#define UFO_OPERATORS_GNSSRO_QC_BACKGROUNDCHECKRONBAM_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// BackgroundCheckRONBAM: check observation closeness to background

class BackgroundCheckRONBAM : public FilterBase,
                             private util::ObjectCounter<BackgroundCheckRONBAM> {
 public:
  static const std::string classname() {return "ufo::BackgroundCheckRONBAM";}

  BackgroundCheckRONBAM(ioda::ObsSpace &, const eckit::Configuration &,
                       std::shared_ptr<ioda::ObsDataVector<int> >,
                       std::shared_ptr<ioda::ObsDataVector<float> >);
  ~BackgroundCheckRONBAM();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::fguess;}
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_QC_BACKGROUNDCHECKRONBAM_H_
