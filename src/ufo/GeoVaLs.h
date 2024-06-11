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

namespace eckit {
  class Configuration;
}

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

// -----------------------------------------------------------------------------

/// \brief GeoVaL format.
///
/// \see GeoVaLs for more information about the formats.
///
/// \note The enumerator values must be kept in sync with the constants defined in
/// `ufo_geovals_mod.F90`.
enum class GeoVaLFormat : int {
  /// \brief Not a format in itself, but a special value signaling that the code should defer to
  /// the format returned by GeoVaLs::defaultFormat().
  ///
  /// The default format can be adjusted by calling GeoVaLs::setDefaultFormat(), either directly or
  /// using a ScopedDefaultGeoVaLFormatChange object.
  DEFAULT = 0,

  /// \brief The sampled format (with one GeoVaL profile per interpolation path).
  SAMPLED = 1,

  /// \brief The reduced format (with one GeoVaL profile per observation location).
  REDUCED = 2
};

std::ostream & operator<<(std::ostream &os, GeoVaLFormat format);

// -----------------------------------------------------------------------------

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

/// \brief Stores values of geophysical variables at observation locations.
///
/// The values of each variable can be stored in up to two formats:
///
/// (a) *Sampled*, i.e. as profiles obtained by direct interpolation of a model field along a
///     set of paths sampling the observation locations. The mapping between locations and paths
///     does not need to be 1-to-1; in particular, each location may be sampled by multiple paths.
///     Different variables may be interpolated along different sets of paths. This format of
///     variable values is most commonly used by observation operators.
///
/// (b) *Reduced*, i.e. computed by replacing each set of profiles obtained by interpolation along
///     the paths sampling a single location by a single profile. This format of variable values
///     is most commonly used by observation filters and bias operators.
///
/// In the simplest and most common case, in which each location is sampled by just a single path,
/// the two formats are identical to each other. The UFO implementation of the GeoVaLs interface
/// automatically detects whether this is the case for given variable and if so, stores only one
/// copy of its values. When there are multiple paths per location, values in the sampled format
/// are typically converted to the reduced format by an appropriate form of averaging.
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

  GeoVaLs(const Locations_ & locations,
          const oops::Variables & vars, const std::vector<size_t> & nlevs);

// Deprecated constructor - Please do not use this constructor in new code.
  GeoVaLs(std::shared_ptr<const ioda::Distribution> dist, const oops::Variables & vars);
// Deprecated constructor - Please do not use this constructor in new code.
  GeoVaLs(const Locations_ & locations, const oops::Variables & vars);
// Constructor for tests - Please do not use this constructor in new code.
  GeoVaLs(const eckit::Configuration &, const ioda::ObsSpace &, const oops::Variables &);

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

  /// \brief Instructs a GeoVaLs object to store some GeoVaLs in the reduced format.
  ///
  /// Space will be allocated for the reduced representation of these variables unless they are
  /// already stored in the sampled format and their reduced and sampled representations are
  /// identical.
  ///
  /// \param vars
  ///   Names of the variables to be stored in the reduced format. They must not already be stored
  ///   in this format.
  /// \param nlevs
  ///   A vector of length `vars.size()` whose ith element indicates how many values per
  ///   location will be stored in the GeoVaL corresponding to the ith variable in `vars`.
  void addReducedVars(const oops::Variables & vars, const std::vector<size_t> & nlevs);

  void zero();
  void reorderzdir(const std::string &, const std::string &);
  void random();
  double rms() const;
  double normalizedrms(const GeoVaLs &) const;

  /// \brief Return true if this GeoVaLs object contains variable `var` stored in format `format`.
  bool has(const oops::Variable & var, GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// \brief Return the list of variables stored in the sampled format.
  const oops::Variables & getVars() const {return vars_;}
  /// \brief Return the list of variables stored in the reduced format.
  const oops::Variables & getReducedVars() const {return reducedVars_;}

  /// \brief Return true if the reduced and sampled formats of variable `var` are stored in the
  /// same block of memory (and therefore are identical).
  bool areReducedAndSampledFormatsAliased(const oops::Variable & var) const;

  /// \brief Return the format of variables accessed by functions such as get...() and put...()
  /// if their `format` parameter is set to GeoVaLFormat::DEFAULT.
  GeoVaLFormat defaultFormat() const;
  /// \brief Set the format of variables accessed by functions such as get...() and put...()
  /// if their `format` parameter is set to GeoVaLFormat::DEFAULT.
  void setDefaultFormat(GeoVaLFormat format);

  /// Return the number of levels in the variable \p var stored in format \p format.
  size_t nlevs(const oops::Variable & var, GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Return the number of profiles in the variable \p var (i.e. the number of paths along which
  /// it has been interpolated) stored in format \p format.
  size_t nprofiles(const oops::Variable & var, GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Get 2D GeoVaLs for variable \p var stored in format \p format (fails for 3D GeoVaLs)
  void get(std::vector<double> &, const oops::Variable & var,
           GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Get 2D GeoVaLs for variable \p var stored in format \p format (fails for 3D GeoVaLs),
  /// and convert to float
  void get(std::vector<float> &, const oops::Variable & var,
           GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Get 2D GeoVaLs for variable \p var stored in format \p format (fails for 3D GeoVaLs),
  /// and convert to int
  void get(std::vector<int> &, const oops::Variable & var,
           GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;

  /// Get GeoVaLs at a specified level
  void getAtLevel(std::vector<double> &, const oops::Variable &, const int,
                  GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Get GeoVaLs at a specified level and convert to float
  void getAtLevel(std::vector<float> &, const oops::Variable &, const int,
                  GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Get GeoVaLs at a specified level and convert to int
  void getAtLevel(std::vector<int> &, const oops::Variable &, const int,
                  GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;

  /// Get a specified profile of the variable `var` stored in format \p format
  ///
  /// For variables stored in the sampled format, each profile contains the results of variable
  /// interpolation along a specific path. For variables stored in the reduced format, each profile
  /// contains values computed by reducing all profiles sampling a specific observation location to
  /// a single profile.
  void getProfile(std::vector<double> &vals, const oops::Variable &var, const int profileIndex,
                  GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Get a specified profile of the variable `var` and convert to float
  void getProfile(std::vector<float> &, const oops::Variable &, const int profileIndex,
                  GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Get a specified profile of the variable `var` and convert to int
  void getProfile(std::vector<int> &, const oops::Variable &, const int profileIndex,
                  GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;

  /// Get GeoVaLs at a specified location
  ///
  /// For variables stored in the sampled format, this function works only if the variable `var`
  /// has been interpolated along one path per observation location; otherwise it throws an
  /// exception. Use the getProfile() function to handle the general case.
  void getAtLocation(std::vector<double> &vals, const oops::Variable &var, const int loc,
                     GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Get GeoVaLs at a specified location and convert to float
  ///
  /// For variables stored in the sampled format, this function works only if the variable `var`
  /// has been interpolated along one path per observation location; otherwise it throws an
  /// exception. Use the getProfile() function to handle the general case.
  void getAtLocation(std::vector<float> &, const oops::Variable &, const int,
                     GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Get GeoVaLs at a specified location and convert to int
  ///
  /// For variables stored in the sampled format, this function works only if the variable `var`
  /// has been interpolated along one path per observation location; otherwise it throws an
  /// exception. Use the getProfile() function to handle the general case.
  void getAtLocation(std::vector<int> &, const oops::Variable &, const int,
                     GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;

  /// Put GeoVaLs at level \p lev for variable \p var of type \c double stored in format \p format.
  void putAtLevel(const std::vector<double> & vals, const oops::Variable & var, const int lev,
                  GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Put GeoVaLs at level \p lev for variable \p var of type \c float stored in format \p format.
  void putAtLevel(const std::vector<float> & vals, const oops::Variable & var, const int lev,
                  GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Put GeoVaLs at level \p lev for variable \p var of type \c int stored in format \p format.
  void putAtLevel(const std::vector<int> & vals, const oops::Variable & var, const int lev,
                  GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;

  /// Store the specified profile of variable \p var stored in format \p format.
  ///
  /// For variables stored in the sampled format, each profile contains the results of variable
  /// interpolation along a specific path. For variables stored in the reduced format, each profile
  /// contains values computed by reducing all profiles sampling a specific observation location to
  /// a single profile.
  void putProfile(const std::vector<double> & vals, const oops::Variable & var,
                  const int profileIndex, GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Store the specified profile of variable \p var stored in format \p format.
  void putProfile(const std::vector<float> & vals, const oops::Variable & var,
                  const int profileIndex, GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Store the specified profile of variable \p var stored in format \p format.
  void putProfile(const std::vector<int> & vals, const oops::Variable & var,
                  const int profileIndex, GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;

  /// Put GeoVaLs at location \p loc for variable \p var of type \c double
  /// stored in format \p format.
  ///
  /// This function works only if the variable `var` has been interpolated along one path per
  /// observation location; otherwise it throws an exception. Use the putProfile()
  /// function to handle the general case.
  void putAtLocation(const std::vector<double> & vals, const oops::Variable & var, const int loc,
                     GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Put GeoVaLs at location \p loc for variable \p var of type \c float
  /// stored in format \p format.
  ///
  /// This function works only if the variable `var` has been interpolated along one path per
  /// observation location; otherwise it throws an exception. Use the putProfile()
  /// function to handle the general case.
  void putAtLocation(const std::vector<float> & vals, const oops::Variable & var, const int loc,
                     GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;
  /// Put GeoVaLs at location \p loc for variable \p var of type \c int
  /// stored in format \p format.
  ///
  /// This function works only if the variable `var` has been interpolated along one path per
  /// observation location; otherwise it throws an exception. Use the putProfile()
  /// function to handle the general case.
  void putAtLocation(const std::vector<int> & vals, const oops::Variable & var, const int loc,
                     GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;

  /// \brief Retrieve a vector mapping the index of each observation location to the range of
  /// indices of the profiles obtained by interpolating the variable `var` along the paths sampling
  /// that location.
  ///
  /// For instance, the profiles obtained by sampling the `i`th location have indices from
  /// `profileIndicesGroupedByLocation[i].begin` up to but not including
  /// `profileIndicesGroupedByLocation[i].end`.
  void getProfileIndicesGroupedByLocation(
      const oops::Variable &var,
      std::vector<util::Range<size_t>> &profileIndicesGroupedByLocation,
      GeoVaLFormat format = GeoVaLFormat::DEFAULT) const;

  void read(const eckit::Configuration &, const ioda::ObsSpace &);
  void write(const eckit::Configuration &) const;

  /// \brief Return the number of observation locations.
  ///
  /// Note that each GeoVaL stored in the sampled format may contain multiple profiles obtained by
  /// interpolating the corresponding model variable along paths sampling the same location, and
  /// the set of interpolation paths may vary from one variable to another. Call nprofiles() to
  /// retrieve the number of profiles in a specific GeoVaL.
  size_t nlocs() const;

  // This always acts on sampled GeoVaLs
  void fill(const oops::Variable &name, const ConstVectorRef<size_t> &indx,
            const ConstMatrixRef<double> &vals, const bool levelsTopDown);
  void fillAD(const oops::Variable &name, const ConstVectorRef<size_t> &indx,
              MatrixRef<double> vals, const bool levelsTopDown) const;

  int & toFortran() {return keyGVL_;}
  const int & toFortran() const {return keyGVL_;}

 private:
  void print(std::ostream &) const;

  /// Convert data stored in a Locations_ object into the form required by
  /// the Fortran GeoVaLs setup subroutines.
  ///
  /// \note `isSamplingMethodTrivial` is declared as a `unique_ptr` to an array rather than as a
  /// vector because we need to pass a `bool*` pointer to Fortran, but the internal representation
  /// of the data held by a `std::vector<bool>` may be different from an array of `bool`s.
  void fillSetupInputs(const Locations_ & locations,
                       size_t & nlocs, std::vector<size_t> & numPathsByMethod,
                       std::vector<size_t> & samplingMethodByVar,
                       std::unique_ptr<bool[]> & isSamplingMethodTrivial) const;

  /// Finish setting up the Fortran GeoVaLs object by letting it know which interpolation paths
  /// sample which observation locations.
  void setupSamplingMethods(const Locations_ & locations);

  // -----------------------------------------------------------------------------
  /*! \brief Take the input vector and recast to type<T> whilst respecting
             missing values */
  template <typename T>
  void cast(const std::vector<double> & doubleVals,
            std::vector<T> & vals) const {
    const T missing = util::missingValue<T>();
    const double missingDouble = util::missingValue<double>();
    std::transform(doubleVals.begin(),
                   doubleVals.end(),
                   vals.begin(),
                   [missing, missingDouble](double val){
                     return val == missingDouble ? missing : static_cast<T>(val);
                   });
  }

  /// \brief If `format` is set to DEFAULT, convert it to ORIGINAL or REDUCED depending on the
  /// value returned by defaultFormat().
  GeoVaLFormat explicitFormat(GeoVaLFormat format) const;

  /// \brief Like explicitFormat(), except that the result is converted to `int` to allow it to be
  /// passed to Fortran.
  int explicitFormatAsInt(GeoVaLFormat format) const;

  F90goms keyGVL_;
  oops::Variables vars_;
  oops::Variables reducedVars_;
  std::shared_ptr<const ioda::Distribution> dist_;   /// observations MPI distribution
};

}  // namespace ufo

#endif  // UFO_GEOVALS_H_
