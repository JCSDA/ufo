/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "ufo/obslocalization/ObsLocalizationBase.h"

#include "oops/util/Printable.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Observation space localization
template<typename ITERATOR>
class ObsLocalization : public util::Printable {
 public:
  ObsLocalization(const eckit::Configuration &, ioda::ObsSpace &);
  ~ObsLocalization();

  void computeLocalization(const ITERATOR &, ioda::ObsVector &) const;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<ObsLocalizationBase<ITERATOR>> obsloc_;
};
// -----------------------------------------------------------------------------

template<typename ITERATOR>
ObsLocalization<ITERATOR>::ObsLocalization(const eckit::Configuration & config,
                                           ioda::ObsSpace & os)
  : obsloc_(ObsLocalizationFactory<ITERATOR>::create(config, os))
{}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
ObsLocalization<ITERATOR>::~ObsLocalization() {
  obsloc_.reset();
}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
void ObsLocalization<ITERATOR>::computeLocalization(const ITERATOR & iter,
                                                    ioda::ObsVector & loc) const {
  obsloc_->computeLocalization(iter, loc);
}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
void ObsLocalization<ITERATOR>::print(std::ostream & os) const {
  os << *obsloc_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
