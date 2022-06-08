# UFO Directory Structure

## The source directory tree (what you get with `git clone ufo`)

- `.github`: Used for issue templates
- `CI`: Used for CI scripts
- `cmake`: CMake-specific files go here. This directory is used to group instructions to simplify the top-level CMakeLists.txt file.
- `docs`: Doxygen documentation. Top-level organization and long documentation pages are found here.
- `share`
  - `ufo`: ufo-specific resource files. This includes scripts consumed by downstream projects and files needed for ufo components to run. These files are not test files.
- `src`: Source files (C++, C, Fortran, Python) and private headers (i.e. not exposed outside of ufo)
  - `ufo`: ufo library sources
    - Base / top-level ufo classes and definitions go here
    - `basis`: Functions for simobs
    - `errors`: Definitions for observation error covariance matrices
    - `filters`
        - General filter source files and headers go here.
        - `actions`: Filter actions go here.
        - `obsfunctions`: ObsFunctions are all found here.
    - `fov`: Field of view calculations for satellite footprints.
    - `obslocalization`: Model interface classes for horizontal and vertical localization of observations.
    - `operators`: Contains all forward, TL, and AD operators (aerosols, crtm, gnssro, identity, scatwind, etc.).
    - `predictors`: For bias correction codes.
    - `profile`: Contains routines applicable to handling vertical profiles of data, such as with radiosondes.
    - `utils`: Utility functions and classes. This folder contains constant definitions, string utilities, distance calculators, routines to determine indices used for variable iterations, and classes for data file manipulation and interpolation.
    - `variabletransforms`: Routines to transform variables. Ex: height from pressure. Ex 2: potential temperature from pressure and temperature.
  - `mains`: Executable program sources
- `test`: Library tests and example applications
  - `mains`: Test program sources
  - `ufo`: Test program headers
  - `testinput`: Test YAMLs
    - `instrument_tests`: Instrument tests that integrate several components of UFO. These tests show how to configure instruments to match operational configurations.
    - `unit_tests`: Tests of individual components of ufo. The layout here roughly matches the layout in the `src` directory (see above).
- `tools`: Repository tools (cpplint.py script, test wrapper scripts, etc.)
- `.clang-format`, `.clang-tidy`, `CPPLINT.cfg`: Linters and code formatting configuration
- `.gitattributes` and `.gitignore`: git-specific configuration files
- `CMakeLists.txt`: Top-level CMake project configuration
- `COPYING` and `LICENSE.md`: Apache license information
- `README.md`: Top-level project description


## The install tree (what you get with `make install`)

The installation directory structure roughly follows the structure specified in the [Filesystem Hierarchy Standard](https://en.wikipedia.org/wiki/Filesystem_Hierarchy_Standard) and the [GNU Coding Standards](https://www.gnu.org/prep/standards/html_node/Directory-Variables.html).

The install tree is organized as follows:

- Root directory (```CMAKE_INSTALL_PREFIX```)
  - `bin`: ufo binaries
  - `include`
    - `ufo`: ufo's C++ and Fortran headers. This corresponds to the headers in the `src` directory.
  - `lib`: ufo's libraries
    - *libufo.[so,dylib]*: The ufo C++/Fortran library
    - `cmake`
      - `ufo`: CMake project definitions used by downstream projects to import ufo using the ```find_package``` command.
    - `python[3.8,3.9,*]`: Python modules for interfacing with ufo (**To be implemented**) 
  - `module`
    - `ufo`: Fortran module files. These are subdivided according to compiler family and compiler version because Fortran modules are compiler-dependent (e.g. `module/ufo/GNU/11.2.0`).
  - `share`
    - `doc`
      - `ufo`
        - `html`: HTML manual pages
    - `ufo`: ufo-specific resource files
    - `man`: GNU man-formatted ufo documentation

