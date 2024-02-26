/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIASPRECONDITIONER_H_
#define UFO_OBSBIASPRECONDITIONER_H_

#include <string>
#include <vector>
#include <boost/noncopyable.hpp>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace ufo {
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------

class ObsBiasPreconditioner : public util::Printable,
                              private boost::noncopyable,
                              private util::ObjectCounter<ObsBiasPreconditioner> {
 public:
  static const std::string classname() {return "ufo::ObsBiasPreconditioner";}

// Constructor, destructor
  explicit ObsBiasPreconditioner(const std::vector<double> &);
  ~ObsBiasPreconditioner() {}

// Linear algebra operators
  void multiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;

 private:
  void print(std::ostream &) const {}
  const std::vector<double> precond_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASPRECONDITIONER_H_
