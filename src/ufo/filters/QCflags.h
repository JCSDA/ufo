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
  constexpr int pass    = 0;  // we like that one!
  constexpr int passive = 1;  // H(x) is computed (for monitoring, BC...) but obs not assimilated
// Single digit values reserved for DA use.
// For now only 0, 1 and >1 are used but keeping space for other potential use cases.

// Actual rejection flags
  constexpr int missing = 10;  // missing values prevent use of observation
  constexpr int preQC   = 11;  // observation rejected by pre-processing
  constexpr int bounds  = 12;  // observation value out of bounds
  constexpr int domain  = 13;  // observation not within domain of use
  constexpr int black   = 14;  // observation black listed
  constexpr int Hfailed = 15;  // H(x) computation failed
  constexpr int thinned = 16;  // observation removed due to thinning
  constexpr int diffref = 17;  // metadata too far from reference
  constexpr int clw     = 18;  // observation removed due to cloud field
  constexpr int fguess  = 19;  // observation too far from guess
  constexpr int seaice  = 20;  // observation based sea ice detection, also flags land points
  constexpr int track   = 21;  // observation removed as inconsistent with the rest of track
  constexpr int buddy   = 22;  // observation rejected by the buddy check
  constexpr int derivative = 23;  // observation removed due to metadata derivative value
  constexpr int profile = 24;  // observation rejected by at least one profile QC check
  constexpr int onedvar = 25;  // observation failed to converge in 1dvar check
  constexpr int bayesianQC = 26;  // observation failed due to Bayesian background check
  constexpr int modelobthresh = 27;  // observation failed modelob threshold check
  constexpr int history = 28;  // observation failed when compared with historical data
  constexpr int processed = 29;  // observation processed but deliberately H(x) not calculated
  constexpr int superrefraction = 30;  // observation rejected by GNSSRO super refraction QC
  constexpr int superob = 31;  // superob value not set at this location
  /// \brief Return true if the QC flag \p qcflag indicates that an observation has been rejected,
  /// false otherwise.
  inline bool isRejected(int qcflag) {
    // TODO(wsmigaj): don't treat passive observations as rejected
    return qcflag != pass;
  }

  /// \brief Return true if the QC flag \p qcflag indicates that an observation is defective.
  ///
  /// This means that the observed value is missing, the observation has been rejected by
  /// pre-processing or its model equivalent cannot be computed.
  inline bool isDefective(int qcflag) {
    return qcflag == missing || qcflag == preQC || qcflag == Hfailed;
  }
};  // namespace QCflags

}  // namespace ufo

#endif  // UFO_FILTERS_QCFLAGS_H_
