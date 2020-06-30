/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileData.h"

namespace ufo {
  ProfileData::ProfileData(ioda::ObsSpace &obsdb,
                           const ProfileConsistencyCheckParameters &options,
                           const ProfileIndices &profileIndices)
    : ProfileDataBase(obsdb, options, profileIndices)
  {
    retrieveAllData();
  }

  void ProfileData::retrieveAllData()
  {
    // Retrieve variables from obsdb
    retrieveDataVector("air_pressure", "MetaData", pressures_);
    retrieveDataVector("air_temperature", "ObsValue", tObs_);
    retrieveDataVector("air_temperature", "HofX", tBkg_);
    retrieveDataVector("relative_humidity", "ObsValue", RHObs_);
    retrieveDataVector("relative_humidity", "HofX", RHBkg_);
    retrieveDataVector("dew_point_temperature", "ObsValue", tdObs_);
    retrieveDataVector("geopotential_height", "ObsValue", zObs_);
    retrieveDataVector("geopotential_height", "HofX", zBkg_);
    retrieveDataVector("eastward_wind", "ObsValue", uObs_);
    retrieveDataVector("northward_wind", "ObsValue", vObs_);
    retrieveDataVector("PstarBackgr", "MetaData", PstarBackgr_);
    retrieveDataVector("station_id", "MetaData", stationID_);
  }

  void ProfileData::fillProfileValues()
  {
    // Fill variable information for a particular profile
    fillProfileData(pressures_, pressures_prof_);
    fillProfileData(tObs_, tObs_prof_);
    fillProfileData(tBkg_, tBkg_prof_);
    fillProfileData(RHObs_, RHObs_prof_);
    fillProfileData(RHBkg_, RHBkg_prof_);
    fillProfileData(tdObs_, tdObs_prof_);
    fillProfileData(zObs_, zObs_prof_);
    fillProfileData(zBkg_, zBkg_prof_);
    fillProfileData(uObs_, uObs_prof_);
    fillProfileData(vObs_, vObs_prof_);
    fillProfileData(PstarBackgr_, PstarBackgr_prof_);
    fillProfileData(stationID_, stationID_prof_);
  }
}  // namespace ufo


