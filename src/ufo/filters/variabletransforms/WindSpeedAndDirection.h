/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_VARIABLETRANSFORMS_WINDSPEEDANDDIRECTION_H_
#define UFO_FILTERS_VARIABLETRANSFORMS_WINDSPEEDANDDIRECTION_H_

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

/// \brief Wind Speed And Direction filter
///
/// The filter performs a variable conversion from the wind components, eastward_wind and
/// northward_wind, to wind_speed and wind_from_direction. The newly calculated variables
/// are included in the same obs space. The filter does not have any configuration options.
///
/// Example:
///
/// \code{.yaml}
/// obs filters:
/// - filter: Wind Speed And Direction
/// \endcode

class WindSpeedAndDirection : public FilterBase,
                              private util::ObjectCounter<WindSpeedAndDirection> {
 public:
  static const std::string classname() {return "ufo::WindSpeedAndDirection";}

  WindSpeedAndDirection(ioda::ObsSpace &, const eckit::Configuration &,
                        std::shared_ptr<ioda::ObsDataVector<int> >,
                        std::shared_ptr<ioda::ObsDataVector<float> >);
  ~WindSpeedAndDirection();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::preQC;}
};

}  // namespace ufo

#endif  // UFO_FILTERS_VARIABLETRANSFORMS_WINDSPEEDANDDIRECTION_H_
