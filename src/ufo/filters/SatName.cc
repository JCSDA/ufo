/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/SatName.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include "ioda/distribution/Accumulator.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

// With floats std::to_string may yield unexpected results, so use stringstream instead
template < typename Type > std::string to_str(const Type & t)
{
  std::ostringstream os;
  os << t;
  return os.str ();
}

std::tuple<std::string, int> get_channel_info(int SatID, float centralFrequency, int satobchannel,
                                              const std::vector<SatIDRangeParameters>
                                              &SatIDRanges) {
  for (const SatIDRangeParameters &SatIDRange : SatIDRanges) {
    if (SatIDRange.minSatID <= SatID && SatID <= SatIDRange.maxSatID) {
      for (const FrequencyBandParameters &frequencyBand : SatIDRange.Satellite_comp.value()) {
        if (frequencyBand.minFrequency <= centralFrequency &&
            centralFrequency <= frequencyBand.maxFrequency &&
            (frequencyBand.satobchannel.value() == boost::none ||
            satobchannel == *frequencyBand.satobchannel.value())) {
          return std::make_tuple(frequencyBand.windChannel,
                                 frequencyBand.windChannelID);
        }
      }
    }
  }
  return std::make_tuple(missing_value_string, missing_value_int);
}
std::string get_sat_name(int SatID, const std::vector<SatIDRangeParameters> &SatIDRanges) {
  for (const SatIDRangeParameters &SatIDRange : SatIDRanges) {
    for (const SatnameParameters &SatNames : SatIDRange.Satellite_id.value()) {
      if (SatID == SatNames.Satnumber.value()) {
        return SatNames.Satname;
      }
    }
  }
  return missing_value_string;
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
/*! \brief A filter that creates a string variable that makes it simpler to
 *  identify Atmospheric Motion Vector (AMV) / Satwind observations by
 *  combining satellite and channel information. This is useful for later
 *  AMV processing where we want to apply filters to subsets of observations.
 *
 *  \details To identify the type of motion that has been tracked, AMV BUFR observations
 *  are supplied with a channel central frequency (Hz) and a wind computation method:
 *
 * 002023 SATELLITE DERIVED WIND COMPUTATION METHOD
 * *  0  Reserved
 * *  1  INFRARED     MOTION OBSERVED IN THE INFRARED CHANNEL
 * *  2  VISIBLE      MOTION OBSERVED IN THE VISIBLE CHANNEL
 * *  3  VAPOUR CLOUD MOTION OBSERVED IN THE WATER VAPOUR CHANNEL
 * *  4  COMBINATION  MOTION OBSERVED IN A COMBINATION OF SPECTRAL CHANNELS
 * *  5  VAPOUR CLEAR MOTION OBSERVED IN THE WATER VAPOUR CHANNEL IN CLEAR AIR
 * *  6  OZONE        MOTION OBSERVED IN THE OZONE CHANNEL
 * *  7  VAPOUR       MOTION OBSERVED IN WATER VAPOUR CHANNEL (CLOUD OR CLEAR)
 * *  13 Root-mean-square
 *
 * The most common use of the wind computation method is to distinguish between clear-sky and
 * cloudy water vapour targets.
 *
 * This filter combines this channel information, together with the satellite name, to
 * create a string that defines the satellite/channel combination of each observation.
 * We also output a diagnostic variable which provides information on unidentified
 * satellites or channels, plus an optional channel number integer.
 *
 * Required :
 * *  "MetaData", "sensorCentralFrequency"
 * *  "MetaData", "satelliteIdentifier"
 * *  "MetaData", "windComputationMethod"
 *
 * Outputs:
 * *  "MetaData", "satwindIdentifier"
 * *  "Diag", "satwindIdentifier"
 * *  "MetaData", "sensorChannelNumber"
 *
 * Example:
 * The following yaml will attempt to identify two infrared channels with computation method
 * value of 1 and central frequencies falling betwen the min and max frequency bounds.
 * If there are observations that match these conditions they are labelled with the respective
 * "wind channel" string.
 * If observations are identified from GOES-16 (platform number 270) they are also labelled with
 * the respective "Sat name" string.
 * This will fill MetaData/satwindIdentifier with values "GOES16ir112","GOES16ir38" if these are present
 * in the observations.
 * If either the satellite or channel are not identified, then MetaData/satwindIdentifier is set to
 * "MISSING". To help track down why observations are set to missing we also output a diagnostic
 * string variable, Diag/satwindIdentifier. When the observation is not identified, this has the form:
 *   id<satellite identifier>_comp<cloud motion method>_freq<central frequency>.
 * E.g. if the satellite is identified but the channel is not: "GOES16_comp3_freq0.484317e14",
 *      if the satellite is not identified but the channel is: "id270ir112".
 *
 *  \code{.unparsed}
 *  obs filters:
 *  - filter: satname
 *    SatName assignments:
 *    - min WMO Satellite id: 1
 *      max WMO Satellite id: 999
 *      Satellite_comp:
 *      - satobchannel: 1
 *        min frequency: 2.6e+13
 *        max frequency: 2.7e+13
 *        wind channel: ir112
 *        wind channel id: 10
 *      - satobchannel: 1
 *        min frequency: 7.5e+13
 *        max frequency: 8.2e+13
 *        wind channel: ir38
 *        wind channel id: 11
 *      Satellite_id:
 *      - Sat ID: 270
 *        Sat name: GOES16
 * \endcode
 *
 */

void SatName::applyFilter(const std::vector<bool> & apply,
                          const Variables & filtervars,
                          std::vector<std::vector<bool>> & flagged) const {
  std::vector<float> cfreq(obsdb_.nlocs());
  std::vector<int> satid(obsdb_.nlocs());
  std::vector<int> compm(obsdb_.nlocs());
  // initialise output vectors to missing data string
  std::vector<std::string> satwind_id(obsdb_.nlocs(), missing_value_string);
  std::vector<std::string> diag_id(obsdb_.nlocs(), missing_value_string);
  std::vector<int> channel_id(obsdb_.nlocs(), missing_value_int);
  // get variables from ObsSpace
  obsdb_.get_db("MetaData", "sensorCentralFrequency", cfreq);
  obsdb_.get_db("MetaData", "satelliteIdentifier", satid);
  obsdb_.get_db("MetaData", "windComputationMethod", compm);
  // define counter variables to be summed over all processors at the end of the routine
  std::unique_ptr<ioda::Accumulator<size_t>> countSatAccumulator =
      obsdb_.distribution()->createAccumulator<size_t>();
  std::unique_ptr<ioda::Accumulator<size_t>> countChanAccumulator =
      obsdb_.distribution()->createAccumulator<size_t>();

  for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
    std::string satellite_name;
    satellite_name = get_sat_name(satid[jobs], parameters_.SatNameAssignments.value());
    const auto[channel_name, channel_number] =
      get_channel_info(satid[jobs], cfreq[jobs], compm[jobs],
                       parameters_.SatNameAssignments.value());
    // Fill channel number
    channel_id[jobs] = channel_number;
    // if both satellite and channel name have been identified, then combine and fill satwind_id
    if (satellite_name != missing_value_string &&
        channel_name != missing_value_string) {
      satwind_id[jobs] = satellite_name + channel_name;
    }
    // if the satellite has not been identified then output the satid number to the diagnostic,
    // otherwise output the found satellite name
    std::string satellite_diag;
    if (satellite_name == missing_value_string) {
      satellite_diag = "id" + to_str(satid[jobs]);
      countSatAccumulator->addTerm(jobs, 1);
    } else {
      satellite_diag = satellite_name;
    }
    // if the channel has not been identified then output the computation method and
    // central frequency to the diagnostic, otherwise output the found channel name
    std::string channel_diag;
    if (channel_name == missing_value_string) {
      channel_diag = "_comp" + to_str(compm[jobs]) +
                     "_freq" + to_str(cfreq[jobs]/1.0e14) + "e14";
      countChanAccumulator->addTerm(jobs, 1);
    } else {
      channel_diag = channel_name;
    }
    // combine diagnostic strings and fill diag_id
    diag_id[jobs] = satellite_diag + channel_diag;
  }
  obsdb_.put_db("MetaData", "satwindIdentifier", satwind_id);
  obsdb_.put_db("Diag", "satwindIdentifier", diag_id);
  obsdb_.put_db("MetaData", "sensorChannelNumber", channel_id);
  // sum number of unidentified satellites and channels
  const std::size_t count_missing_sat = countSatAccumulator->computeResult();
  const std::size_t count_missing_chan = countChanAccumulator->computeResult();
  oops::Log::info() << "SatName: " << count_missing_sat
                     << " observations with unidentified satellite id" << std::endl;
  oops::Log::info() << "SatName: " << count_missing_chan
                     << " observations with unidentified channel" << std::endl;
  }
// -----------------------------------------------------------------------------
void SatName::print(std::ostream & os) const {
  os << "SatName: config = " << parameters_ << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace ufo
