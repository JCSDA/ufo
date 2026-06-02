/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */
#ifndef UFO_FILTERS_SATNAME_H_
#define UFO_FILTERS_SATNAME_H_

#include <algorithm>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "oops/util/missingValues.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}
namespace ufo {
  static const std::string missing_value_string = util::missingValue<std::string>();
  static const int missing_value_int = util::missingValue<int>();
//
// table in yaml file to relate satellite name and wmo number label in yaml "Satellite_id"
//
class SatnameParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SatnameParameters, Parameters)
 public:
  oops::RequiredParameter<std::string> Satname{"Sat name", this};
  oops::RequiredParameter<float> Satnumber{"Sat ID", this};
};
//
// various satillite instrument characteristics label in yaml "Satellite_comp"
//
class FrequencyBandParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FrequencyBandParameters, Parameters)
 public:
  oops::RequiredParameter<float> minFrequency{"min frequency", this};
  oops::RequiredParameter<float> maxFrequency{"max frequency", this};
  oops::OptionalParameter<int> satobchannel{"satobchannel", this};
  oops::RequiredParameter<std::string> windChannel{"wind channel", this};
  oops::Parameter<int> windChannelID{"wind channel id", missing_value_int, this};
};
//
//  range of satellite wmo numbers "min obs type" and "max obs type"
//
class SatIDRangeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SatIDRangeParameters, Parameters)
 public:
  oops::RequiredParameter<int> minSatID{"min WMO Satellite id", this};
  oops::RequiredParameter<int> maxSatID{"max WMO Satellite id", this};
  oops::Parameter<std::vector<FrequencyBandParameters>> Satellite_comp{
      "Satellite_comp", {}, this};
  oops::Parameter<std::vector<SatnameParameters>> Satellite_id{
      "Satellite_id", {}, this};
};
// Parameters controlling the operation of the SatName filter.
class SatNameParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(SatNameParameters, FilterParametersBase)

 public:
  oops::Parameter<std::vector<SatIDRangeParameters>>
    SatNameAssignments{ "SatName assignments", {}, this};};
class SatName : public FilterBase,
 private util::ObjectCounter<SatName> {
 public:
// The type of parameters accepted by the constructor of this filter.
// This typedef is used by the FilterFactory.
  typedef SatNameParameters Parameters_;
  static const std::string classname() {return "ufo::SatName";}
  SatName(ioda::ObsSpace &, const Parameters_ &,
                 std::shared_ptr<ioda::ObsDataVector<int> >,
                 std::shared_ptr<ioda::ObsDataVector<float> >);
  ~SatName();
 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
  std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::pass;}
  Parameters_ parameters_;
};
}  // namespace ufo
#endif  // UFO_FILTERS_SATNAME_H_
