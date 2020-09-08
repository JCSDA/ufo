/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_OBSFILTERDATA_H_
#define UFO_FILTERS_OBSFILTERDATA_H_

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/Variable.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/*! \brief ObsFilterData provides access to all data related to an ObsFilter
 *
 * \details ObsFilterData can always provide access to all data from ObsSpace
 * and optionally to data from H(x) ObsVector, GeoVaLs and ObsDiagnostics.
 * The latter three can be associated with ObsFilterData by using associate()
 * method.
 *
 */
class ObsFilterData : public util::Printable,
                      private util::ObjectCounter<ObsFilterData> {
 public:
  static const std::string classname() {return "ufo::ObsFilterData";}

  //! Constructs ObsFilterData and associates ObsSpace with it
  explicit ObsFilterData(ioda::ObsSpace &);
  ~ObsFilterData();

  //! Associates GeoVaLs with this ObsFilterData
  void associate(const GeoVaLs &);
  //! Associates H(x) ObsVector with this ObsFilterData
  void associate(const ioda::ObsVector &, const std::string &);
  //! Associates ObsDiagnostics from ObsOperator with this ObsFilterData
  void associate(const ObsDiagnostics &);
  //! Associates ObsDataVector with this ObsFilterData (float)
  void associate(const ioda::ObsDataVector<float> &, const std::string &);
  //! Associates ObsDataVector with this ObsFilterData (int)
  void associate(const ioda::ObsDataVector<int> &, const std::string &);

  //! Gets requested data from ObsFilterData
  void get(const Variable &, std::vector<float> &) const;
  //! Gets requested data at requested level from ObsFilterData
  void get(const Variable &, const int, std::vector<float> &) const;
  //! Gets requested data from ObsFilterData
  void get(const Variable &, std::vector<std::string> &) const;
  //! Gets requested data from ObsFilterData
  void get(const Variable &, std::vector<int> &) const;
  //! Gets requested data from ObsFilterData
  void get(const Variable &, std::vector<util::DateTime> &) const;
  //! Gets requested data from ObsFilterData (ObsDataVector has to be allocated)
  void get(const Variable &, ioda::ObsDataVector<float> &) const;
  //! Gets requested data from ObsFilterData (ObsDataVector has to be allocated)
  void get(const Variable &, ioda::ObsDataVector<int> &) const;
  //! Checks if requested data exists in ObsFilterData
  bool has(const Variable &) const;

  //! Determines dtype of the provided variable
  ioda::ObsDtype dtype(const Variable &) const;

  //! Returns number of locations
  size_t nlocs() const;
  //! Returns number of levels for specified variable if 3D GeoVaLs or ObsDiags
  size_t nlevs(const Variable &) const;
  //! Returns reference to ObsSpace associated with ObsFilterData
  ioda::ObsSpace & obsspace() const {return obsdb_;}
  //! Returns reference to GeoVaLs required by 1DVar
  const GeoVaLs * getGeoVaLs() const {return gvals_;}
 private:
  void print(std::ostream &) const;
  bool hasVector(const std::string &, const std::string &) const;
  bool hasDataVector(const std::string &, const std::string &) const;
  bool hasDataVectorInt(const std::string &, const std::string &) const;

  ioda::ObsSpace & obsdb_;                 //!< ObsSpace associated with this object
  const GeoVaLs mutable * gvals_;          //!< pointer to GeoVaLs associated with this object
  std::map<std::string, const ioda::ObsVector *> ovecs_;  //!< Associated ObsVectors
  const ObsDiagnostics mutable * diags_;   //!< pointer to ObsDiagnostics associated with object
  std::map<std::string, const ioda::ObsDataVector<float> *> dvecsf_;  //!< Associated ObsDataVectors
  std::map<std::string, const ioda::ObsDataVector<int> *> dvecsi_;  //!< Associated ObsDataVectors
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFILTERDATA_H_
