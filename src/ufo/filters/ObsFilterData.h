/*
 * (C) Copyright 2019-2021 UCAR
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

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/filters/DiagnosticFlag.h"

namespace util {
  class DateTime;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  enum class ObsDtype;
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;
  class Variable;

// -----------------------------------------------------------------------------
/*! \brief ObsFilterData provides access to all data related to an ObsFilter
 *
 * \details ObsFilterData can always provide access to all data from ObsSpace, values computed on
 * the fly by ObsFunctions, and optionally to data from the H(x) ObsVector, GeoVaLs and
 * ObsDiagnostics. The latter three can be associated with ObsFilterData by using associate()
 * method.
 */
class ObsFilterData : public util::Printable,
                      private util::ObjectCounter<ObsFilterData> {
 public:
  static const std::string classname() {return "ufo::ObsFilterData";}

  //! Constructs ObsFilterData and associates an ObsSpace with it.
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

  //! \brief Fills a `std::vector` with values of the specified variable.
  //!
  //! \param varname
  //!   The requested variable.
  //! \param[out] values
  //!   Vector to be filled with values of the requested variable.
  //!
  //! An exception is thrown if the requested variable does not exist or is not of the correct type.
  void get(const Variable &varname, std::vector<float> &values,
           bool skipDerived = false) const;
  //! \overload
  void get(const Variable &varname, std::vector<int> &values,
           bool skipDerived = false) const;
  //! \overload
  void get(const Variable &varname, std::vector<std::string> &values,
           bool skipDerived = false) const;
  //! \overload
  void get(const Variable &varname, std::vector<util::DateTime> &values,
           bool skipDerived = false) const;
  //! \overload
  void get(const Variable &varname, std::vector<DiagnosticFlag> &values,
           bool skipDerived = false) const;

  //! \brief Fills a `std::vector` with values of the specified variable at a single level.
  //!
  //! \param varname
  //!   The requested variable, which must belong to one of the following groups: GeoVaLs, ObsDiag
  //!   and ObsBiasTerm.
  //! \param level
  //!   Requested level.
  //! \param[out]
  //!   Vector to be filled with values of the requested variable.
  //!
  //! An exception is thrown if the requested variable does not exist or if it is not in one of
  //! the groups listed above.
  void get(const Variable & varname, const int level,
           std::vector<float> & values) const;
  //! \overload
  void get(const Variable & varname, const int level,
           std::vector<double> & values) const;

  //! brief Fills a `ioda::ObsDataVector` with values of the specified variable.
  //!
  //! \param varname
  //!   The requested variable.
  //! \param[out] values
  //!   An `ObsDataVector` to be filled with values of the requested variable. Must be
  //!   pre-populated, i.e. contain a row corresponding to each requested channel of that variable.
  //!
  //! An exception is thrown if the requested variable does not exist or is not of the correct type.
  void get(const Variable &varname, ioda::ObsDataVector<float> &values,
           bool skipDerived = false) const;
  //! \overload
  void get(const Variable &varname, ioda::ObsDataVector<int> &values,
           bool skipDerived = false) const;
  //! \overload
  void get(const Variable &varname, ioda::ObsDataVector<std::string> &values,
           bool skipDerived = false) const;
  //! \overload
  void get(const Variable &varname, ioda::ObsDataVector<util::DateTime> &values,
           bool skipDerived = false) const;
  //! \overload
  void get(const Variable &varname, ioda::ObsDataVector<DiagnosticFlag> &values,
           bool skipDerived = false) const;

  //! Returns true if variable `varname` is known to ObsFilterData, false otherwise.
  bool has(const Variable &varname) const;

  //! Determines dtype of the provided variable
  ioda::ObsDtype dtype(const Variable &) const;

  //! \brief Returns the number of locations in the associated ObsSpace.
  size_t nlocs() const;
  //! \brief Returns the number of levels in the specified variable.
  //!
  //! This is useful primarily for variables belonging to the GeoVaLs, ObsDiag and ObsBiasTerm
  //! groups. For all other variables this function returns 1.
  size_t nlevs(const Variable &) const;
  //! Returns reference to ObsSpace associated with ObsFilterData
  ioda::ObsSpace & obsspace() const {return obsdb_;}
  //! Returns reference to GeoVaLs required by 1DVar
  const GeoVaLs * getGeoVaLs() const {return gvals_;}
  //! Returns reference to ObsDiagnostics
  const ObsDiagnostics * getObsDiags() const {return diags_;}
 private:
  void print(std::ostream &) const;
  bool hasVector(const std::string &, const std::string &) const;
  bool hasDataVector(const std::string &, const std::string &) const;
  bool hasDataVectorInt(const std::string &, const std::string &) const;

  template <typename T>
  void getVector(const Variable &varname, std::vector<T> &values,
                 bool skipDerived = false) const;
  /// Called by the overloads of get() taking an ioda::ObsDataVector of strings or datetimes.
  template <typename T>
  void getNonNumeric(const Variable &varname, ioda::ObsDataVector<T> &values,
                     bool skipDerived = false) const;

  ioda::ObsSpace & obsdb_;                 //!< ObsSpace associated with this object
  const GeoVaLs mutable * gvals_;          //!< pointer to GeoVaLs associated with this object
  std::map<std::string, const ioda::ObsVector *> ovecs_;  //!< Associated ObsVectors
  const ObsDiagnostics mutable * diags_;   //!< pointer to ObsDiagnostics associated with object
  std::map<std::string, const ioda::ObsDataVector<float> *> dvecsf_;  //!< Associated ObsDataVectors
  std::map<std::string, const ioda::ObsDataVector<int> *> dvecsi_;  //!< Associated ObsDataVectors
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFILTERDATA_H_
