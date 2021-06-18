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
                             ioda::ObsSpace& os,
                             const std::shared_ptr<ioda::ObsDataVector<int>>& flags,
                             const std::vector<bool> &apply)
    : options_(options), obsdb_(os) , flags_(*flags), apply_(apply) {
  method_ = formulas::resolveMethods(options.Method.value());
  formulation_  = formulas::resolveFormulations(options.Formulation.value(),
                                                options.Method.value());
  UseValidDataOnly_ = options.UseValidDataOnly.value();
  AllowSuperSaturation_ = options.AllowSuperSaturation.value();
  obsName_ = os.obsname();
}

void TransformBase::filterObservation(const std::string &variableName,
                                       std::vector<float> &obsVector) const {
  if (flags_.has(variableName)) {
    const float missing = missingValueFloat;
    const std::vector<int> *varFlags = &flags_[variableName];

    std::transform(obsVector.begin(), obsVector.end(),  // Input range 1
                   varFlags->begin(),  // First element of input range vector 2 (must be same size)
                   obsVector.begin(),  // First element of output range (must be same size)
                   [missing](float obsvalue, int flag)
                   { return flag == QCflags::missing || flag == QCflags::bounds
                     ? missing : obsvalue; });
  }
}

void TransformBase::getObservation(const std::string &originalTag, const std::string &varName,
                                   std::vector<float> &obsVector, bool require) const {
  const size_t nlocs = obsdb_.nlocs();

  if (obsdb_.has("DerivedValue", varName)) {
    obsVector = std::vector<float>(nlocs);
    obsdb_.get_db("DerivedValue", varName, obsVector);
    // Set obsValue to missingValueFloat if flag is equal to QCflags::missing or QCflags::bounds
    if (UseValidDataOnly()) filterObservation(varName, obsVector);
  } else if (obsdb_.has(originalTag, varName)) {
    obsVector = std::vector<float>(nlocs);
    obsdb_.get_db(originalTag, varName, obsVector);
    // Set obsValue to missingValueFloat if flag is equal to QCflags::missing or QCflags::bounds
    if (UseValidDataOnly()) filterObservation(varName, obsVector);
  }

  if (require && obsVector.empty()) {
    throw eckit::BadValue("The parameter `" + varName + "@" + originalTag +
                          "` does not exist in the ObsSpace ", Here());
  }
}

TransformFactory::TransformFactory(const std::string& name) {
  if (getMakers().find(name) != getMakers().end())
    throw eckit::BadParameter(
        name + " already registered in ufo::TransformFactory.", Here());

  getMakers()[name] = this;
}

std::unique_ptr<TransformBase> TransformFactory::create(
    const std::string& name, const VariableTransformsParameters& options,
    ioda::ObsSpace& os, const std::shared_ptr<ioda::ObsDataVector<int>>& flags,
    const std::vector<bool>& apply) {

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

  std::unique_ptr<TransformBase> ptr = jloc->second->make(options, os, flags, apply);
  oops::Log::trace() << "TransformBase::create done" << std::endl;

  return ptr;
}
}  // namespace ufo
