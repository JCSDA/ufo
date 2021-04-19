/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
//-------------------------------------------------------------------------------
//
// This filter is the first step in the AMV processing.
// It creates a  string for classifying   Atmospheric Motion Vectors (AMV)'s into groups.
// This is then used for later AMV processing modules
//
// Inputs:
// data from an AMV  WMO bufr file in netcdf4 format
//
// Required :
//
// "MetaData", "sensor_central_frequency", cfreq
// "MetaData", "originating_subcentre", orsub
// "MetaData", "satellite_identifier", satid
// "MetaData", "wind_computation_method", compm
//
// Outputs:
// "MetaData", "satwind_id", wind_id
//-------------------------------------------------------------------------------
//
// Method:
//
// There are a number of satellite series (both geostationary and polar) that currently provide AMVs
// Each satellite or satellite pairs, provides AMVs using difference channels in the visible
// infrared and water vapour bands.
//
//
// satellite frequency and computational method are used to generate the AMV character string,
// an example is:
// "NOAA16" and "ir108"  are combined "NOAA16ir108" and added to the output NETCDF4 file
// as a variable satwind_id
//
//
//                 satids channel values
//
//                 hrvis vis06 vis08
//                 ir16 ir38 ir87 ir97 ir108 ir120
//                 wv62 wv67 wv73
//                 cswv62 cswv67 cswv73
//                 mixwv62 mixwv67 mixwv73
//                 misstyp           unknown satellite identified
//
//
#include "ufo/filters/SatName.h"
#include <algorithm>
#include <map>
#include <set>
#include <vector>
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/StringUtils.h"
namespace ufo {
std::string Sat_Characteristics(int SatID, float centralFrequency, int satobchannel,
                        const std::vector<SatIDRangeParameters> &SatIDRanges) {
    std::string Missing_Sat_windChan = "misstyp";
    for (const SatIDRangeParameters &SatIDRange : SatIDRanges)
        if (SatIDRange.minSatID <= SatID && SatID <= SatIDRange.maxSatID)
        for (const FrequencyBandParameters &frequencyBand : SatIDRange.Satellite_comp.value())
        if (frequencyBand.minFrequency <= centralFrequency &&
           centralFrequency <= frequencyBand.maxFrequency &&
           (frequencyBand.satobchannel.value() == boost::none ||
            satobchannel == *frequencyBand.satobchannel.value())) {
        return frequencyBand.windChannel;
  }
      return Missing_Sat_windChan;
}
std::string sat_id(int SatID,
                   const std::vector<SatIDRangeParameters> &SatIDRanges) {
    std::string satmn = "misatid";
    for (const SatIDRangeParameters &SatIDRange : SatIDRanges)
      for (const SatnameParameters &SatNames : SatIDRange.Satellite_id.value())
        if (SatID == SatNames.Satnumber.value()) {return SatNames.Satname;}
    return satmn;
}
// -----------------------------------------------------------------------------
SatName::SatName(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                               std::shared_ptr<ioda::ObsDataVector<int> > flags,
                               std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::debug() << "SatName: config (constructor) = " << parameters_ << std::endl;
}
// -----------------------------------------------------------------------------
SatName::~SatName() {}
// -----------------------------------------------------------------------------
void SatName::applyFilter(const std::vector<bool> & apply,
                                 const Variables & filtervars,
                                 std::vector<std::vector<bool>> & flagged) const {
  std::vector<float> cfreq(obsdb_.nlocs());
  std::vector<int> orsub(obsdb_.nlocs());
  std::vector<int> orgin(obsdb_.nlocs());
  std::vector<int> satid(obsdb_.nlocs());
  std::vector<int> compm(obsdb_.nlocs());
  std::vector<std::string> wind_id(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "sensor_central_frequency", cfreq);
  obsdb_.get_db("MetaData", "originating_subcentre", orsub);
  obsdb_.get_db("MetaData", "originating_centre", orgin);
  obsdb_.get_db("MetaData", "satellite_identifier", satid);
  obsdb_.get_db("MetaData", "wind_computation_method", compm);
//
  for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
    wind_id[jobs] = sat_id(satid[jobs], parameters_.SatNameAssignments.value())
                      + Sat_Characteristics(satid[jobs], cfreq[jobs], compm[jobs],
                         parameters_.SatNameAssignments.value());
  }
  obsdb_.put_db("MetaData", "satwind_id", wind_id);
//  only print the first 10 observations while testing
  for (size_t ii = 0; ii < 10 ; ++ii) {
     oops::Log::trace() << " freq " << cfreq[ii] << "orsub " <<  orsub[ii] << "orgin "
       <<  orgin[ii] << " satid " << satid[ii] << " compm " << compm[ii]
       << " wind id  "<< wind_id[ii] << std::endl;}
}
// -----------------------------------------------------------------------------
void SatName::print(std::ostream & os) const {
  os << "SatName: config = " << parameters_ << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace ufo
