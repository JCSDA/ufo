# PPRO: Parameterized Polarimetric Radar Operator

## Overview

PPRO is the Parameterized Polarimetric Radar Operator integrated into JEDI-UFO for computing dual-polarization radar variables from atmospheric model state variables.

> **Note**: This is the built-in UFO version. For standalone library usage, see the external [ppro-lib](https://github.com/JCSDA-internal/ppro-lib) repository.

## Scientific Background

The P-PRO operator provides two polarimetric forward operators:

**Zhang21 Method** (Zhang et al. 2021): Assumes double-moment microphysics schemes with **exponential** particle size distribution (PSD). Dual-pol radar quantities are parameterized as polynomial functions of mean diameter (Dm) and mass content of hydrometeors (Wx), by fitting to integral of **T-matrix**-based backscattering amplitude over the PSD. Two melting models are implemented in this operator.

**TCWA2 Method** (Tsai et al. 2026): Similar to Zhang21, but assumes a **gamma** PSD. Dual-pol radar quantities are parameterized as more complex functions of gamma-PSD parameters (λ, α). Environmental state (e.g., T), density and shape of particles, and model-provided melting fractions (smlf/gmlf) are also used in the simulation. The single scattering calculation is based on the **Rayleigh** approximation, following Ryzhkov et al., 2011.

## Operator Comparison: Zhang21 vs TCWA2

| Feature               | Zhang21                                            | TCWA2                                                    |
|-----------------------|----------------------------------------------------|----------------------------------------------------------|
| **PSD assumption**    | Exponential                                        | Gamma                                                    |
| **Scattering method** | T-matrix                                           | Rayleigh approximation                                   |
| **Fitted functions**  | Polynomials of Dm and Wx (mass content)            | Complex functions of more parameters                     |
| **Microphysics**      | Thompson, WSM6, NSSL, TCWA2                        | TCWA2 (Taiwan CWA)                                       |
| **Melting treatment** | Melting models (Zhang21 or Liu et al. 2024)        | Melted fraction from TCWA2 MP scheme                     |
| **Required inputs**   | qr, qs, qg [+ qh, nr, ns, ng, nh, vg, vh for NSSL] | qr, qs, qg, qc, qi, nr, ns, ng, ni, smlf, gmlf           |
| **Ice crystal**       | Not used                                           | Explicit (qi, ni)                                        |
| **Cloud water**       | Not used                                           | Explicit (qc)                                            |

**Key Algorithmic Differences**:
- **Zhang21**: Exponential PSD + T-matrix scattering → Dm, Wx → Polynomial evaluation → Zh, ZDR, KDP, ρhv
- **TCWA2**: Gamma PSD + Rayleigh scattering → λ, α, Dx, smlf/gmlf → Fitted analytic functions → Zh, ZDR, KDP

**Recommendation**: 
- Use **Zhang21** for multiple microphysics schemes (e.g., Thompson, WSM6, NSSL, TCWA2)
- Use **TCWA2** for Taiwan CWA double-moment microphysics scheme TCWA2

The operator computes dual-polarization radar variables including:
- **Zhh**: Horizontal reflectivity (dBZ)
- **ZDR**: Differential reflectivity (dB)
- **KDP**: Specific differential phase (deg/km)
- **ρhv**: Copolar correlation coefficient

## Features

- **Multiple Polarimetric Operators:**
  - **Zhang21**: Exponential PSD, T-matrix scattering, polynomial functions (Zhang et al. 2021)
  - **TCWA2**: Gamma PSD, Rayleigh scattering, complex functions (Tsai et al. 2026)
  
- **Microphysics Schemes:**
  - **Thompson**: Double-moment for rain (with Zhang21 operator)
  - **WSM6**: Single-moment (with Zhang21 operator)
  - **NSSL**: Double-moment with hail (with Zhang21 operator)
  - **TCWA2**: Taiwan CWA double-moment microphysics scheme (Tsai et al. 2026, works with both TCWA2 and Zhang21 operator)
  
- **Radar Frequencies:**
  - S-band (11 cm wavelength) and C-band (5.3 cm wavelength) support
  
- **Advanced Physics:**
  - Melting layer treatment (Zhang21: melting model, Liu et al. 2024; TCWA2: melting fraction from TCWA2 MP scheme)
  - Pure and melting hydrometeor categories (rain, snow, graupel, hail, ice)
  - Gamma PSD (Particle Size Distribution) parameterization (TCWA2)

## Directory Structure

```
ppro/
├── CMakeLists.txt              # CMake build configuration
├── README.md                   # This file
├── PPRO_CONFIG_GUIDE.md        # User configuration guide for YAML setup
├── ObsPPRO.h                   # C++ UFO operator interface (header)
├── ObsPPRO.cc                  # C++ UFO operator implementation
├── Zhang21Forward.h            # Zhang21 forward operator (header)
├── Zhang21Forward.cc           # Zhang21 forward operator (implementation)
├── TCWA2Forward.h              # TCWA2 forward operator (header)
├── TCWA2Forward.cc             # TCWA2 forward operator (implementation)
└── PPROCoefsData.h             # Zhang21 polynomial coefficients (embedded)
```

## Usage in JEDI-UFO

See `PPRO_CONFIG_GUIDE.md` for detailed YAML configuration examples.

### Basic Configuration

```yaml
obs operator:
  name: PPRO
  polarimetric operator: <Zhang21|TCWA2>
  microphysics option: <Thompson|WSM6|NSSL|TCWA2>
  melting scheme: liu24                     # For Zhang21 (recommended)
  VertCoord: height_above_mean_sea_level
```

## Module Structure

### Main Operator (`ObsPPRO.h` / `ObsPPRO.cc`)

C++ UFO operator providing the standard interface:
- constructor: parses YAML config, builds the per-instance tuning and the
  required GeoVaLs list, and selects the Zhang21 or TCWA2 operator
- `simulateObs()`: main computation (dispatches to the selected operator per point)
- `print()`: logs operator type, microphysics option, and band

### Zhang21 Operator (`Zhang21Forward.h` / `.cc`)

Exponential PSD + T-matrix polynomial operator:
- `computePoint()`: compute using polynomial functions of Dm and Wx
- `meltingModelLiu24()` / `meltingModelJung08()`: melting layer treatment
- per-hydrometeor PSD/Z helpers (`dmZwsm6()`, `dmZ2moment()`, ...)

### TCWA2 Operator (`TCWA2Forward.h` / `.cc`)

Gamma PSD + Rayleigh scattering fitted formulations:
- `dualpolOpRainTcwa2()`, `dualpolOpIceTcwa2()`, `dualpolOpSnowTcwa2()`,
  `dualpolOpGraupTcwa2()`: per-hydrometeor analytic formulas
- `tcwa2ComputePoint()`: aggregates the per-hydrometeor contributions

## Coefficient Files

Both **Zhang21** and **TCWA2** operators have all coefficients embedded directly in the source code (`PPROCoefsData.h`) - no external coefficient files are required at runtime. This ensures portability and eliminates path configuration issues.

## Version History

- **2026-06**: Full C++ rewrite of the operator from Fortran (Zhang21 and TCWA2), removing the Fortran sources and C-Fortran interface (Rong Kong)
- **2026-04**: Integrated back into UFO as built-in module (removed external library dependency)
- **2025-04**: Multi-operator architecture refactoring, modularization with external ppro-lib library
- **2025-01**: Added TCWA2 support
- **2024**: Initial PPRO implementation with Zhang21, Thompson, WSM6, NSSL support

## Authors

**Lead Developer:**
- **Rong Kong** (NCAR/MMM) - C++ rewrite of the PPRO operator from Fortran (both Zhang21 and TCWA2); library architecture, modularization, multi-operator framework, operator tuning, testing and maintenance

**Operator Implementation (original Fortran):**
- **Zhiquan (Jake) Liu** (NCAR/MMM) - Zhang21 operator implementation
- **Tzu-Chin Tsai** (CWA) - TCWA2 operator implementation

**Extensions and Integration:**
- **Hejun Xie** (NCAR/MMM) - Integration of Zhang21 into the JEDI-UFO framework
- **Tao Sun** (NCAR/MMM) - Hail categories and extended microphysics support

## License

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

## References

1. **Zhang21 Operator**:  
   Zhang, G., J. Gao, and M. Du, 2021: Parameterized forward operators for simulation and assimilation of polarimetric radar data with numerical weather predictions. *Adv. Atmos. Sci.*, **38**(5), 737−754.

2. **Zhang21 Melting Scheme**:  
   Liu, et al., 2024: A New Melting Model and Its Implementation in Parameterized Forward Operators for Polarimetric Radar Data Simulation With Double Moment Microphysics Schemes. *JGR Atmospheres*, **129**.

3. **Application of Zhang21 Operator in MPAS-JEDI**:  
   Kong et al., 2026: Assimilation of Radar Radial Velocity and Polarimetric Observations Using LETKF within MPAS-JEDI: Development of Assimilation Capabilities and Test with an Afternoon Thunderstorm in Taiwan. *Mon. Wea. Rev.*, submitted.

4. **TCWA2 Operator**:  
   Tsai, T.-C., J.-P. Chen, Z. Liu, S.-Y. Jiang, R. Kong, Y.-J. Wu, J. Ban, L.-F. Hsiao, Y.-S. Tang, P.-L. Chang, and J.-S. Hong, 2026: Development of the TCWA2 Bulk Cloud Microphysics Scheme and Its Integration with a Dual-Polarization Radar Operator for Forecasting Applications. *J. Adv. Model. Earth Syst.*, submitted.

## Contact

For questions, issues, or contributions:
- **Technical Lead**: Rong Kong (rkong@ucar.edu)
