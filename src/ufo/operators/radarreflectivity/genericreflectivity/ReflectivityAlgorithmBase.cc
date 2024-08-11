/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/radarreflectivity/genericreflectivity/ReflectivityAlgorithmBase.h"

#include <map>
#include <memory>
#include <string>

#include "oops/util/Logger.h"

namespace ufo {

ReflectivityAlgorithmBase::ReflectivityAlgorithmBase
(const ReflectivityAlgorithmParametersBase & params,
 const ioda::ObsSpace & odb,
 const int idxRefl,
 oops::Variables & reqvars,
 oops::Variables & reqvarsTL)
  : obsdb_(odb),
    idxReflectivity_(idxRefl)
{
  // Subclasses of this class should fill `reqvars` and `reqvarsTL` here.
}

void ReflectivityAlgorithmBase::simulateObs(const GeoVaLs & gv,
                                            ioda::ObsVector & ovec,
                                            ObsDiagnostics & diag,
                                            const QCFlags_t & flag) const {
  oops::Log::trace() << "ReflectivityAlgorithmBase::simulateObs starting" << std::endl;
  this->simulateObsImpl(gv, ovec, diag, flag);
  oops::Log::trace() << "ReflectivityAlgorithmBase::simulateObs done" << std::endl;
}

void ReflectivityAlgorithmBase::setTrajectory(const GeoVaLs & gv,
                                              ObsDiagnostics & diag,
                                              const QCFlags_t & flag) {
  oops::Log::trace() << "ReflectivityAlgorithmBase::setTrajectory starting" << std::endl;
  this->setTrajectoryImpl(gv, diag, flag);
  oops::Log::trace() << "ReflectivityAlgorithmBase::setTrajectory done" << std::endl;
}

void ReflectivityAlgorithmBase::simulateObsTL(const GeoVaLs & dx,
                                              ioda::ObsVector & dy,
                                              const QCFlags_t & flag) const {
  oops::Log::trace() << "ReflectivityAlgorithmBase::simulateObsTL starting" << std::endl;
  this->simulateObsTLImpl(dx, dy, flag);
  oops::Log::trace() << "ReflectivityAlgorithmBase::simulateObsTL done" << std::endl;
}

void ReflectivityAlgorithmBase::simulateObsAD(GeoVaLs & dx,
                                              const ioda::ObsVector & dy,
                                              const QCFlags_t & flag) const {
  oops::Log::trace() << "ReflectivityAlgorithmBase::simulateObsAD starting" << std::endl;
  this->simulateObsADImpl(dx, dy, flag);
  oops::Log::trace() << "ReflectivityAlgorithmBase::simulateObsAD done" << std::endl;
}

void ReflectivityAlgorithmBase::print(std::ostream & os) const {
  this->printImpl(os);
}

ReflectivityAlgorithmFactory::ReflectivityAlgorithmFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end())
    throw eckit::BadParameter(name +
                              " already registered in ufo::ReflectivityAlgorithmFactory.", Here());
  getMakers()[name] = this;
}

std::unique_ptr<ReflectivityAlgorithmBase>
ReflectivityAlgorithmFactory::create
(const ReflectivityAlgorithmParametersBase & params,
 const ioda::ObsSpace & obsdb,
 const int idxRefl,
 oops::Variables & reqvars,
 oops::Variables & reqvarsTL) {
  oops::Log::trace() << "ReflectivityAlgorithmBase::create starting" << std::endl;
  const std::string & name = params.reflectivityAlgorithmName;
  typename std::map<std::string, ReflectivityAlgorithmFactory*>::iterator jloc =
    getMakers().find(name);
  if (jloc == getMakers().end()) {
    std::string makerNameList;
    for (const auto & makerDetails : getMakers()) makerNameList += "\n  " + makerDetails.first;
    throw eckit::BadParameter(name + " does not exist in ufo::ReflectivityAlgorithmFactory. "
                              "Possible values:" + makerNameList, Here());
  }
  std::unique_ptr<ReflectivityAlgorithmBase> ptr =
    jloc->second->make(params, obsdb, idxRefl, reqvars, reqvarsTL);
  oops::Log::trace() << "ReflectivityAlgorithmBase::create done" << std::endl;
  return ptr;
}

std::unique_ptr<ReflectivityAlgorithmParametersBase>
ReflectivityAlgorithmFactory::createParameters(const std::string & name) {
  oops::Log::trace() << "ReflectivityAlgorithmBase::createParameters starting" << std::endl;
  typename std::map<std::string, ReflectivityAlgorithmFactory*>::iterator jloc =
    getMakers().find(name);
  if (jloc == getMakers().end()) {
    std::string makerNameList;
    for (const auto & makerDetails : getMakers()) makerNameList += "\n  " + makerDetails.first;
    throw eckit::BadParameter(name + " does not exist in ufo::ReflectivityAlgorithmFactory. "
                              "Possible values:" + makerNameList, Here());
  }
  std::unique_ptr<ReflectivityAlgorithmParametersBase> ptr = jloc->second->makeParameters();
  oops::Log::trace() << "ReflectivityAlgorithmBase::createParameters done" << std::endl;
  return ptr;
}

}  // namespace ufo
