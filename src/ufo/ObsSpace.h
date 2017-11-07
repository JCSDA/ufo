/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSSPACE_H_
#define UFO_OBSSPACE_H_

#include <map>
#include <ostream>
#include <string>

#include "oops/interface/ObsSpaceBase.h"
#include "util/DateTime.h"
#include "util/Logger.h"
#include "util/Printable.h"

#include "Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ufo {
  class ObsVec;

/// Wrapper around ObsHelpQG, mostly to hide the factory
class ObsSpace : public oops::ObsSpaceBase {

 public:
  ObsSpace(const eckit::Configuration &, const util::DateTime &, const util::DateTime &);
  ObsSpace(const ObsSpace &);
  ~ObsSpace();

  void getdb(const std::string & col, int & keyData) const {
//    helper_->getdb(obsname_, col, keyData);
  }
  void putdb(const std::string & col, const int & keyData) const {
//    helper_->putdb(obsname_, col, keyData);
  }

  int locations(const util::DateTime & t1, const util::DateTime & t2) const {
    int key_locs;
//    key_locs = helper_->locations(obsname_, t1, t2);
    return key_locs;
  }

  void generateDistribution(const eckit::Configuration & conf) {
//    helper_->generateDistribution(conf, obsname_, winbgn_, winend_, nobs_);
  }

//  void printJo(const ObsVec &, const ObsVec &);

  int nobs() const {return nobs_;}
  const std::string & obsname() const {return obsname_;}
  const util::DateTime & windowStart() const {return winbgn_;}
  const util::DateTime & windowEnd() const {return winend_;}

  int & toFortran() {return keyOspace_;}
  const int & toFortran() const {return keyOspace_;}

 private:
  void print(std::ostream &) const;
  ObsSpace & operator= (const ObsSpace &);
  std::string obsname_;
  unsigned int nobs_;
  const util::DateTime winbgn_;
  const util::DateTime winend_;
  F90Ospace keyOspace_;

  static std::map < std::string, int > theObsFileCount_;
};

}  // namespace ufo

#endif  // UFO_OBSSPACE_H_
