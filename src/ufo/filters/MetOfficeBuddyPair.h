/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_METOFFICEBUDDYPAIR_H_
#define UFO_FILTERS_METOFFICEBUDDYPAIR_H_

namespace ufo {

/// \brief Properties of a pair of observations checked against each other during buddy check.
struct MetOfficeBuddyPair {
  MetOfficeBuddyPair(int obsIdA_, int obsIdB_, double distanceInKm_,
            double rotationAInRad_, double rotationBInRad_)
    : obsIdA(obsIdA_), obsIdB(obsIdB_), distanceInKm(distanceInKm_),
      rotationAInRad(rotationAInRad_), rotationBInRad(rotationBInRad_)
  {}

  int obsIdA;  //< Id of the first observation
  int obsIdB;  //< Id of the second observation
  double distanceInKm;   //< Horizontal distance between observations (in km)
  double rotationAInRad;  //< Direction of ob B from ob A (in rad)
  double rotationBInRad;  //< Reciprocal direction of ob A from B (in rad)
};

}  // namespace ufo

#endif  // UFO_FILTERS_METOFFICEBUDDYPAIR_H_
