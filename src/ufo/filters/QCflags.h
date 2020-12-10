/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_QCFLAGS_H_
#define UFO_FILTERS_QCFLAGS_H_

namespace ufo {

namespace QCflags {
  constexpr int pass    = 0;   // we like that one!
  constexpr int missing = 1;   // missing values prevent use of observation
  constexpr int preQC   = 2;   // observation rejected by pre-processing
  constexpr int bounds  = 3;   // observation value out of bounds
  constexpr int domain  = 4;   // observation not within domain of use
  constexpr int black   = 5;   // observation black listed
  constexpr int Hfailed = 6;   // H(x) computation failed
  constexpr int thinned = 7;   // observation removed due to thinning
  constexpr int diffref = 8;   // metadata too far from reference
  constexpr int clw     = 9;   // observation removed due to cloud field
  constexpr int fguess  = 10;  // observation too far from guess
  constexpr int seaice  = 11;  // observation based sea ice detection, also flags land points
  constexpr int track   = 12;  // observation removed as inconsistent with the rest of track
  constexpr int buddy   = 13;  // observation rejected by the buddy check
  constexpr int derivative = 14;  // observation removed due to metadata derivative value
  constexpr int profile = 15;  // observation rejected by at least one profile QC check
  constexpr int ratioref = 16;   // ratio of two values outside of range
};  // namespace QCflags

}  // namespace ufo

#endif  // UFO_FILTERS_QCFLAGS_H_
