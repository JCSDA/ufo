/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_GETSCALARORFILTERDATA_H_
#define UFO_FILTERS_GETSCALARORFILTERDATA_H_


#include <string>
#include <vector>

namespace ufo {

class ObsFilterData;
/// Function to fill in a vector with either a scalar or data from ObsFilterData
//  - if input string contains a float, output vector would be filled in with that number
//    for nlocs size (e.g. if string is "4.0", output vector would contain nlocs x 4.0)
//  - otherwise output vector would contain data from ObsFilterData (e.g. if string is
//    MyFunction@Function, output vector would contatin that. If there's no MyFunction@Function
//    in ObsFilterData, the function will abort.
//  To be used in error inflation, thresholds for BackgroundCheck, etc
std::vector<float> getScalarOrFilterData(const std::string &, const ObsFilterData &);

}  // namespace ufo

#endif  // UFO_FILTERS_GETSCALARORFILTERDATA_H_
