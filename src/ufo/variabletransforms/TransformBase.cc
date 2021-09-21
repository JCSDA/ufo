/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <map>
#include <string>
#include "oops/util/Logger.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

TransformBase::TransformBase(const VariableTransformsParameters& options,
                             const ObsFilterData& data,
                             const std::shared_ptr<ioda::ObsDataVector<int>>& flags)
                             : options_(options), data_(data) , flags_(*flags) {
  method_ = formulas::resolveMethods(options.Method.value());
  formulation_  = formulas::resolveFormulations(options.Formulation.value(),
                                                options.Method.value());
  UseValidDataOnly_ = options.UseValidDataOnly.value();
  AllowSuperSaturation_ = options.AllowSuperSaturation.value();
  obsName_ = data_.obsspace().obsname();
}

TransformFactory::TransformFactory(const std::string& name) {
  if (getMakers().find(name) != getMakers().end())
    throw eckit::BadParameter(
        name + " already registered in ufo::TransformFactory.", Here());

  getMakers()[name] = this;
}

std::unique_ptr<TransformBase> TransformFactory::create(
    const std::string& name, const VariableTransformsParameters& options,
    const ObsFilterData& data,
    const std::shared_ptr<ioda::ObsDataVector<int>>& flags) {

  oops::Log::trace() << "          --> TransformFactory::create" << std::endl;
  oops::Log::trace() << "              --> name: " << name << std::endl;

  typename std::map<std::string, TransformFactory*>::iterator jloc =
      getMakers().find(name);

  if (jloc == getMakers().end()) {
    std::string makerNameList;
    for (const auto& makerDetails : getMakers())
      makerNameList += "\n  " + makerDetails.first;
    std::cout << "       --> makerNameList" << makerNameList << std::endl;
    std::cout << "                        " << name
              << " does not exist in ufo::TransformFactory. "
              << "Possible values:" << makerNameList << std::endl;

    throw eckit::BadParameter(name +
                                  " does not exist in ufo::TransformFactory. "
                                  "Possible values:" +
                                  makerNameList,
                              Here());
  }

  std::unique_ptr<TransformBase> ptr = jloc->second->make(options, data, flags);
  oops::Log::trace() << "TransformBase::create done" << std::endl;

  return ptr;
}
}  // namespace ufo
