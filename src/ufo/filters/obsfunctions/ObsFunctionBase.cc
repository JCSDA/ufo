/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"

#include <map>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "ufo/filters/Variable.h"

namespace ufo {

const char *ObsFunctionTraits<float>::valueTypeName = "float";
const char *ObsFunctionTraits<int>::valueTypeName = "int";
const char *ObsFunctionTraits<std::string>::valueTypeName = "std::string";
const char *ObsFunctionTraits<util::DateTime>::valueTypeName = "util::DateTime";

const char *ObsFunctionTraits<float>::groupName = "ObsFunction";
const char *ObsFunctionTraits<int>::groupName = "IntObsFunction";
const char *ObsFunctionTraits<std::string>::groupName = "StringObsFunction";
const char *ObsFunctionTraits<util::DateTime>::groupName = "DateTimeObsFunction";

// -----------------------------------------------------------------------------

template <typename FunctionValue>
ObsFunctionFactory<FunctionValue>::ObsFunctionFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw eckit::UserError(name + " already registered in ufo::ObsFunctionFactory<" +
                           std::string(ObsFunctionTraits<FunctionValue>::valueTypeName) + ">",
                           Here());
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
ObsFunctionBase<FunctionValue> * ObsFunctionFactory<FunctionValue>::create(const Variable & var) {
  oops::Log::trace() << "ObsFunctionBase::create starting" << std::endl;
  typename std::map<std::string, ObsFunctionFactory*>::iterator jloc =
     getMakers().find(var.variable());
  if (jloc == getMakers().end()) {
    throw eckit::UserError(var.variable() + " does not exist in ufo::ObsFunctionFactory<" +
                           std::string(ObsFunctionTraits<FunctionValue>::valueTypeName) + ">",
                           Here());
  }
  ObsFunctionBase<FunctionValue> * ptr = jloc->second->make(var.options());
  oops::Log::trace() << "ObsFunctionFactory::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
bool ObsFunctionFactory<FunctionValue>::functionExists(const std::string & name) {
  return (getMakers().find(name) != getMakers().end());
}

// -----------------------------------------------------------------------------

// Explicit instantiations for the supported value types
template class ObsFunctionFactory<float>;
template class ObsFunctionFactory<int>;
template class ObsFunctionFactory<std::string>;
template class ObsFunctionFactory<util::DateTime>;

// -----------------------------------------------------------------------------

}  // namespace ufo
