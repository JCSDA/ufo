/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_QCFLAGS_H_
#define UFO_QCFLAGS_H_

namespace ufo {

namespace QCflags {
  constexpr int pass    = 0;   // we like that one!
  constexpr int missing = 1;   // missing values prevent use of observation
  constexpr int preQC   = 2;   // observation rejected by pre-processing
  constexpr int bounds  = 3;   // observation value out of bounds
  constexpr int domain  = 4;   // observation not within domain of use
  constexpr int black   = 5;   // observation black listed
  constexpr int Hfailed = 6;   // H(x) computation failed
  constexpr int fguess  = 10;  // observation too far from guess
};

}  // namespace ufo

#endif  // UFO_QCFLAGS_H_
