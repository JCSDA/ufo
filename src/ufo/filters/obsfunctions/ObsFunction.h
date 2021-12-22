/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTION_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTION_H_

#include <memory>

#include <boost/noncopyable.hpp>

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
}

namespace ufo {
  class ObsFilterData;
  template <typename FunctionValue> class ObsFunctionBase;
  class Variable;
  class Variables;

// -----------------------------------------------------------------------------

/// \brief A function of observation data.
///
/// \tparam FunctionValue
///   Type of the values produced by the function. Must be `float`, `int`, `std::string`
///   or `util::DateTime`.
template <typename FunctionValue>
class ObsFunction : private boost::noncopyable {
 public:
/// constructor takes function name (for factory) on input
  explicit ObsFunction(const Variable &);
  ~ObsFunction();

/// compute(metadata, obs values, output)
  void compute(const ObsFilterData &,
               ioda::ObsDataVector<FunctionValue> &) const;
/// required variables
  const ufo::Variables & requiredVariables() const;

 private:
  std::unique_ptr<ObsFunctionBase<FunctionValue>> obsfct_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTION_H_
