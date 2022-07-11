/*
 * (C) Copyright 2021 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_ACCEPTLIST_H_
#define UFO_FILTERS_ACCEPTLIST_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// \brief Parameters controlling the AcceptList filter.
class AcceptListParameters : public FilterParametersBase {
  // We call this macro instead of "plain" OOPS_CONCRETE_PARAMETERS() because we want to customize
  // the constructor definition (passing "accept" as the default action name to the base class
  // constructor).
  OOPS_CONCRETE_PARAMETERS_ENABLE_COPY_AND_MOVE(AcceptListParameters, FilterParametersBase)

 public:
  /// Set the default action to `accept`.
  explicit AcceptListParameters(oops::Parameters* parent = nullptr)
  : FilterParametersBase(parent, "accept")
  {}

 private:
  // This filter doesn't take any extra parameters.
};

/// \brief A filter that, by default, performs the `accept` action on observations selected by the
/// `where` clause, i.e. resets the QC flags of these observations to `pass` unless they are
/// currently set to 'missing', 'preQC' or 'Hfailed'.
class AcceptList : public FilterBase,
                   private util::ObjectCounter<AcceptList> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef AcceptListParameters Parameters_;

  static const std::string classname() {return "ufo::AcceptList";}

  AcceptList(ioda::ObsSpace &, const Parameters_ &,
             std::shared_ptr<ioda::ObsDataVector<int> >,
             std::shared_ptr<ioda::ObsDataVector<float> >);

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  // "black" seems the right flag to use in the unlikely case of this filter's action being set to
  // "reject" (which is currently the only action that looks at the value returned by
  // this function).
  int qcFlag() const override {return QCflags::black;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_ACCEPTLIST_H_
