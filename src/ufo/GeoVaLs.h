/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_GEOVALS_H_
#define UFO_GEOVALS_H_

#include <Eigen/Core>

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/missingValues.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

#include "ufo/Fortran.h"

namespace oops {
  template <typename OBS> class Locations;
}

namespace util {
  template <typename T> struct Range;
}

namespace ioda {
  class ObsSpace;
  class Distribution;
}

namespace ufo {
struct ObsTraits;

/// \brief Parameters controlling GeoVaLs read/write
class GeoVaLsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeoVaLsParameters, Parameters)

 public:
  /// Filename for I/O.
  /// Note: this parameter is optional because filename does not need
  /// to be specified for GeoVaLs ctor when Variables are empty.
  oops::OptionalParameter<std::string> filename{"filename", this};
  /// A multiplier for how many locations in the geovals file per
  /// a single location in the obs file. There needs to be at least
  /// loc_multiplier * obs_all_nlocs locations in the geovals file.
  oops::Parameter<int> loc_multiplier{"loc_multiplier", 1, this};
  /// Flip GeoVals levels after reading the geoval file when levels_are_top_down is false.
  oops::Parameter<bool> levels_are_top_down{"levels_are_top_down", true, this};
};

// -----------------------------------------------------------------------------

/// GeoVaLs: geophysical values at locations
///
/// Instances of this class store values of (model) variables interpolated along collections of
/// paths sampling the observation locations. Different variables may be interpolated along
/// different sets of paths.
class GeoVaLs : public util::Printable,
                private util::ObjectCounter<GeoVaLs> {
  /// \brief A reference to a read-only vector-valued expression.
  ///
  /// For example, an Eigen::Vector or an Eigen::Map (the latter can be used as a view onto
  /// a chunk of memory stored in another container, such as a std::vector).
  template <typename T>
  using ConstVectorRef = Eigen::Ref<const Eigen::Vector<T, Eigen::Dynamic>>;

  /// \brief A reference to a read-only matrix-valued expression.
  ///
  /// For example, an Eigen::Matrix or an Eigen::Map (the latter can be used as a view onto
  /// a chunk of memory stored in another container, such as a std::vector).
  template <typename T>
  using ConstMatrixRef = Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;

  /// \brief A reference to a writable matrix-valued expression.
  ///
  /// For example, an Eigen::Matrix or an Eigen::Map (the latter can be used as a view onto
  /// a chunk of memory stored in another container, such as a std::vector).
  template <typename T>
  using MatrixRef = Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;

 public:
  typedef GeoVaLsParameters Parameters_;
  typedef oops::Locations<ObsTraits> Locations_;

  static const std::string classname() {return "ufo::GeoVaLs";}

  GeoVaLs(const Locations_ &, const oops::Variables &,
          const std::vector<size_t> &);

// Deprecated constructor - Please do not use this constructor in new code.
  GeoVaLs(std::shared_ptr<const ioda::Distribution>, const oops::Variables &);
// Deprecated constructor - Please do not use this constructor in new code.
  GeoVaLs(const Locations_ &, const oops::Variables &);
// Constructor for tests - Please do not use this constructor in new code.
  GeoVaLs(const Parameters_ &, const ioda::ObsSpace &, const oops::Variables &);

  GeoVaLs(const GeoVaLs &, const int &);
  GeoVaLs(const GeoVaLs &);

  ~GeoVaLs();

  GeoVaLs & operator = (const GeoVaLs &);
  GeoVaLs & operator *= (const double);
  GeoVaLs & operator *= (const std::vector<float> &);
  GeoVaLs & operator += (const GeoVaLs &);
  GeoVaLs & operator -= (const GeoVaLs &);
  GeoVaLs & operator *= (const GeoVaLs &);
  double dot_product_with(const GeoVaLs &) const;

  /// \brief Deprecated method. Allocates GeoVaLs for \p vars variables with
  /// \p nlev number of levels
  /// \details Please do not use in any new code. This method is currently
  /// only used for ObsDiagnostics and will be removed soon. Rely on
  ///
  ///     GeoVaLs(const Locations &, const oops::Variables &,
  ///             const std::vector<size_t> &)
  ///
  /// to allocate GeoVaLs correctly.
  /// Fails if at least one of the \p vars doesn't exist in GeoVaLs.
  /// Only allocates variables that haven't been allocated before.
  /// Fails if one of \p vars is already allocated with number of levels
  /// different than \p nlev; doesn't reallocate variables that are already
  /// allocated with \p nlev.
  void allocate(const int & nlev, const oops::Variables & vars);

  void zero();
  void reorderzdir(const std::string &, const std::string &);
  void random();
  double rms() const;
  double normalizedrms(const GeoVaLs &) const;

  bool has(const std::string & var) const {return vars_.has(var);}
  const oops::Variables & getVars() const {return vars_;}

  size_t nlevs(const std::string & var) const;
  /// Return the number of profiles in the variable \p var (i.e. the number of paths along which
  /// it has been interpolated).
  size_t nprofiles(const std::string & var) const;
  /// Get 2D GeoVaLs for variable \p var (fails for 3D GeoVaLs)
  void get(std::vector<double> &, const std::string & var) const;
  /// Get 2D GeoVaLs for variable \p var (fails for 3D GeoVaLs), and convert to float
  void get(std::vector<float> &, const std::string & var) const;
  /// Get 2D GeoVaLs for variable \p var (fails for 3D GeoVaLs), and convert to int
  void get(std::vector<int> &, const std::string & var) const;

  /// Get GeoVaLs at a specified level
  void getAtLevel(std::vector<double> &, const std::string &, const int) const;
  /// Get GeoVaLs at a specified level and convert to float
  void getAtLevel(std::vector<float> &, const std::string &, const int) const;
  /// Get GeoVaLs at a specified level and convert to int
  void getAtLevel(std::vector<int> &, const std::string &, const int) const;

  /// Get a specified profile of the variable `var` (i.e. the result of its interpolation along a
  /// specified path)
  void getProfile(std::vector<double> &vals, const std::string &var, const int profileIndex) const;
  /// Get a specified profile of the variable `var` and convert to float
  void getProfile(std::vector<float> &, const std::string &, const int profileIndex) const;
  /// Get a specified profile of the variable `var` and convert to int
  void getProfile(std::vector<int> &, const std::string &, const int profileIndex) const;

  /// Get GeoVaLs at a specified location
  ///
  /// This function works only if the variable `var` has been interpolated along one path per
  /// observation location; otherwise it throws an exception. Use the getProfile()
  /// function to handle the general case.
  void getAtLocation(std::vector<double> &vals, const std::string &var, const int loc) const;
  /// Get GeoVaLs at a specified location and convert to float
  ///
  /// This function works only if the variable `var` has been interpolated along one path per
  /// observation location; otherwise it throws an exception. Use the getProfile()
  /// function to handle the general case.
  void getAtLocation(std::vector<float> &, const std::string &, const int) const;
  /// Get GeoVaLs at a specified location and convert to int
  ///
  /// This function works only if the variable `var` has been interpolated along one path per
  /// observation location; otherwise it throws an exception. Use the getProfile()
  /// function to handle the general case.
  void getAtLocation(std::vector<int> &, const std::string &, const int) const;

  /// Put GeoVaLs for double variable \p var at level \p lev.
  void putAtLevel(const std::vector<double> & vals, const std::string & var, const int lev) const;
  /// Put GeoVaLs for float variable \p var at level \p lev.
  void putAtLevel(const std::vector<float> & vals, const std::string & var, const int lev) const;
  /// Put GeoVaLs for int variable \p var at level \p lev.
  void putAtLevel(const std::vector<int> & vals, const std::string & var, const int lev) const;

  /// Store the specified profile of variable \p var (i.e. the values of this variable along the
  /// specified interpolation path)
  void putProfile(const std::vector<double> & vals, const std::string & var,
                  const int profileIndex) const;
  /// Store the specified profile of variable \p var (i.e. the values of this variable along the
  /// specified interpolation path)
  void putProfile(const std::vector<float> & vals, const std::string & var,
                  const int profileIndex) const;
  /// Store the specified profile of variable \p var (i.e. the values of this variable along the
  /// specified interpolation path)
  void putProfile(const std::vector<int> & vals, const std::string & var,
                  const int profileIndex) const;

  /// Put GeoVaLs for double variable \p var at location \p loc.
  ///
  /// This function works only if the variable `var` has been interpolated along one path per
  /// observation location; otherwise it throws an exception. Use the putProfile()
  /// function to handle the general case.
  void putAtLocation(const std::vector<double> & vals, const std::string & var,
                     const int loc) const;
  /// Put GeoVaLs for float variable \p var at location \p loc.
  ///
  /// This function works only if the variable `var` has been interpolated along one path per
  /// observation location; otherwise it throws an exception. Use the putProfile()
  /// function to handle the general case.
  void putAtLocation(const std::vector<float> & vals, const std::string & var, const int loc) const;
  /// Put GeoVaLs for int variable \p var at location \p loc.
  ///
  /// This function works only if the variable `var` has been interpolated along one path per
  /// observation location; otherwise it throws an exception. Use the putProfile()
  /// function to handle the general case.
  void putAtLocation(const std::vector<int> & vals, const std::string & var, const int loc) const;

  /// \brief Retrieve a vector mapping the index of each observation location to the range of
  /// indices of the profiles obtained by interpolating the variable `var` along the paths sampling
  /// that location.
  ///
  /// For instance, the profiles obtained by sampling the `i`th location have indices from
  /// `profileIndicesGroupedByLocation[i].begin` up to but not including
  /// `profileIndicesGroupedByLocation[i].end`.
  void getProfileIndicesGroupedByLocation(
      const std::string &var,
      std::vector<util::Range<size_t>> &profileIndicesGroupedByLocation) const;

  void read(const Parameters_ &, const ioda::ObsSpace &);
  void write(const Parameters_ &) const;

  /// \brief Return the number of observation locations.
  ///
  /// Note that each GeoVaL may contain multiple profiles obtained by interpolating the
  /// corresponding model variable along paths sampling the same location, and the set of
  /// interpolation paths may vary from one variable to another. Call nprofiles() to retrieve the
  /// number of profiles in a specific GeoVaL.
  size_t nlocs() const;

  void fill(const std::string &name, const ConstVectorRef<size_t> &indx,
            const ConstMatrixRef<double> &vals, const bool levelsTopDown);
  void fillAD(const std::string &name, const ConstVectorRef<size_t> &indx,
              MatrixRef<double> vals, const bool levelsTopDown) const;

  int & toFortran() {return keyGVL_;}
  const int & toFortran() const {return keyGVL_;}

 private:
  void print(std::ostream &) const;

  /// Convert data stored in a Locations_ object into the form required by
  /// the Fortran GeoVaLs setup subroutines.
  void fillSetupInputs(const Locations_ & locations,
                       size_t & nlocs, std::vector<size_t> & numPathsByMethod,
                       std::vector<size_t> & samplingMethodByVar) const;

  /// Finish setting up the Fortran GeoVaLs object by letting it know which interpolation paths
  /// sample which observation locations.
  void setupSamplingMethods(const Locations_ & locations);

  // -----------------------------------------------------------------------------
  /*! \brief Take the input vector and recast to type<T> whilst respecting
             missing values */
  template <typename T>
  void cast(const std::vector<double> & doubleVals,
            std::vector<T> & vals) const {
    const T missing = util::missingValue(T());
    const double missingDouble = util::missingValue(double());
    std::transform(doubleVals.begin(),
                   doubleVals.end(),
                   vals.begin(),
                   [missing, missingDouble](double val){
                     return val == missingDouble ? missing : static_cast<T>(val);
                   });
  }

  F90goms keyGVL_;
  oops::Variables vars_;
  std::shared_ptr<const ioda::Distribution> dist_;   /// observations MPI distribution
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_GEOVALS_H_
