# PPRO Operator Configuration Guide

This document describes how to configure the PPRO (Parameterized Polarimetric Radar Operator) in JEDI-UFO for different microphysics schemes and polarimetric operators.

## Overview

The PPRO operator supports:
- **Polarimetric Operators**: 
  - Zhang21 (polynomial functions with coefficients from T-matrix simulations)
  - TCWA2 (analytical gamma distribution formulation)
- **Microphysics Schemes**: Thompson, WSM6, NSSL, TCWA2
- **Radar Bands**: S-band (iband=1), C-band (iband=2)
- **Output Variables**: Zhh, ZDR, KDP, ρhv

## Basic Configuration Structure

```yaml
obs operator:
  name: PPRO
  polarimetric operator: <Zhang21|TCWA2>    # Choose operator type
  microphysics option: <Thompson|WSM6|NSSL|TCWA2>
  melting scheme: liu24                     # Liu et al. 2024 melting model (for Zhang21)
  VertCoord: height_above_mean_sea_level
  # Variable mappings (depends on microphysics scheme)
```

---

## Configuration Examples by Microphysics Scheme

### 1. WSM6 (Single-Moment) with Zhang21 Operator

**Required GeoVaLs: 7 variables**
- Basic: air_density, air_temperature, air_pressure, specific_humidity
- Hydrometeors: qr, qs, qg

```yaml
obs operator:
  name: PPRO
  polarimetric operator: Zhang21        # Polynomial-based from T-matrix (default)
  microphysics option: WSM6
  melting scheme: liu24                 # Liu et al. 2024 melting model (recommended)
  VertCoord: height_above_mean_sea_level
  var_rain_mixing_ratio: qr              # rain water mixing ratio (kg/kg)
  var_snow_mixing_ratio: qs              # snow mixing ratio (kg/kg)
  var_graupel_mixing_ratio: qg           # graupel mixing ratio (kg/kg)
```

**Notes:**
- Simplest configuration with 7 GeoVaLs
- Uses fixed intercept parameters N₀
- Zhang21 operator uses polynomial functions of Dm and fx fitted from T-matrix simulations
- `melting scheme: liu24` enables the Liu et al. (2024) melting layer parameterization
- `melting scheme: jung08` selects the Jung et al. (2008) melting model: the
  rain–ice mixture fraction follows Eq. (2) of Jung et al. (2008),
  F = Fmax·[min(qx/qr, qr/qx)]^0.3, so the mixture is maximal when the rain and
  ice mixing ratios are equal and vanishes when either goes to zero

---

### 2. Thompson (Double-Moment for Rain) with Zhang21 Operator

**Required GeoVaLs: 8 variables**
- WSM6 variables (7) + nr

```yaml
obs operator:
  name: PPRO
  polarimetric operator: Zhang21         # Polynomial-based from T-matrix (default)
  microphysics option: Thompson
  melting scheme: liu24                  # Liu et al. 2024 melting model (recommended)
  VertCoord: height_above_mean_sea_level
  var_rain_mixing_ratio: qr
  var_snow_mixing_ratio: qs
  var_graupel_mixing_ratio: qg
  var_rain_number_concentration: nr      # rain number concentration (#/kg)
```

**Notes:**
- Adds rain number concentration for improved rain representation
- 8 total GeoVaLs required
- Zhang21 operator handles both single-moment (WSM6) and double-moment (Thompson) schemes
- `melting scheme: liu24` enables the Liu et al. (2024) melting layer parameterization

---

### 3. NSSL (Double-Moment with Hail) with Zhang21 Operator

**Required GeoVaLs: 14 variables**
- Basic (4) + mass mixing ratios (4) + number concentrations (4) + volume mixing ratios (2)

```yaml
obs operator:
  name: PPRO
  polarimetric operator: Zhang21         # Polynomial-based from T-matrix (default)
  microphysics option: NSSL
  melting scheme: liu24                  # Liu et al. 2024 melting model (recommended)
  VertCoord: height_above_mean_sea_level
  # Mass mixing ratios
  var_rain_mixing_ratio: qr
  var_snow_mixing_ratio: qs
  var_graupel_mixing_ratio: qg
  var_hail_mixing_ratio: qh              # hail mixing ratio (kg/kg)
  # Number concentrations
  var_rain_number_concentration: nr       # #/kg
  var_snow_number_concentration: ns       # #/kg
  var_graupel_number_concentration: ng    # #/kg
  var_hail_number_concentration: nh       # #/kg
  # Volume mixing ratios
  var_graupel_vol_mixing_ratio: vg       # m³/kg
  var_hail_vol_mixing_ratio: vh          # m³/kg
```

**Notes:**
- Most comprehensive scheme with 14 GeoVaLs
- Includes hail category for severe weather
- Requires volume mixing ratios for graupel and hail
- `melting scheme: liu24` enables the Liu et al. (2024) melting layer parameterization

---

### 4. TCWA2 (Taiwan CWA Analytical Scheme)

**Required GeoVaLs: 15 variables**
- Basic (4) + mass mixing ratios (5) + number concentrations (4) + melted fractions (2)

```yaml
obs operator:
  name: PPRO
  polarimetric operator: TCWA2           # MUST use TCWA2 operator!
  microphysics option: TCWA2             # MUST use TCWA2 microphysics!
  VertCoord: height_above_mean_sea_level
  # Mass mixing ratios
  var_rain_mixing_ratio: qr              # kg/kg
  var_snow_mixing_ratio: qs              # kg/kg
  var_graupel_mixing_ratio: qg           # kg/kg
  var_cloud_mixing_ratio: qc             # cloud water (kg/kg)
  var_ice_mixing_ratio: qi               # ice crystal (kg/kg)
  # Number concentrations
  var_rain_number_concentration: nr       # #/kg
  var_snow_number_concentration: ns       # #/kg
  var_graupel_number_concentration: ng    # #/kg
  var_ice_number_concentration: ni        # #/kg
  # Melted fractions (unique to TCWA2)
  var_snow_melted_fraction: smlf         # snow melting fraction [0-1]
  var_graupel_melted_fraction: gmlf      # graupel melting fraction [0-1]
```

**Notes:**
- **CRITICAL**: Must use `polarimetric operator: TCWA2` AND `microphysics option: TCWA2`
- Uses analytical gamma distribution formulation (no polynomial coefficients needed)
- Requires melted fraction variables for melting layer treatment
- Total 15 GeoVaLs required

**TCWA2 Features:**
- No coefficient files needed (analytical formulation)
- Requires melting layer information (smlf, gmlf)
- Designed for Taiwan CWA cloud microphysics scheme
- Complete ice crystal treatment (qi, ni)

---

## GeoVaLs Summary Table

| Scheme    | # GeoVaLs | Basic | qr,qs,qg | qh | qc,qi | nr,ns,ng | nh | ni | vg,vh | smlf,gmlf |
|-----------|-----------|-------|----------|----|----|----------|----|----|-------|-----------|
| **WSM6**     | 7         | 4     | ✓        |    |    |          |    |    |       |           |
| **Thompson** | 8         | 4     | ✓        |    |    | nr only  |    |    |       |           |
| **NSSL**     | 14        | 4     | ✓        | ✓  |    | ✓        | ✓  |    | ✓     |           |
| **TCWA2**    | 15        | 4     | ✓        |    | ✓  | ✓        |    | ✓  |       | ✓         |

**Legend:**
- Basic (4): air_density, air_temperature, air_pressure, specific_humidity
- qr/qs/qg/qh/qc/qi: mass mixing ratios (kg/kg)
- nr/ns/ng/nh/ni: number concentrations (#/kg)
- vg/vh: volume mixing ratios (m³/kg)
- smlf/gmlf: melted fractions [0-1]

---

## Operator Type Selection

### Zhang21 (Polynomial-Based from T-Matrix Simulations)
- **Method**: Polynomial functions of Dm (mean diameter) and fx (water fraction) fitted from T-matrix simulations
- **Use with**: WSM6, Thompson, NSSL
- **Melting scheme**: `melting scheme: liu24` (Liu et al. 2024, recommended)
- **Pros**: Validated, includes melting layer treatment, efficient polynomial evaluation
- **Coefficients**: Embedded in source (`PPROCoefsData.h`); no external coefficient files needed at runtime
- **Config**: `polarimetric operator: Zhang21` (or omit, it's default)

### TCWA2 (Analytical Gamma Distribution)
- **Method**: Analytical formulation based on gamma distribution PSD parameters
- **Use with**: TCWA2 microphysics **ONLY**
- **Pros**: No polynomial coefficients needed, direct gamma PSD calculation
- **Cons**: Requires melted fraction fields explicitly
- **Config**: `polarimetric operator: TCWA2` (REQUIRED!)
- **Note**: TCWA2 does not compute ρhv from scattering; if `copolarCorrelationCoefficient`
  is requested it returns a placeholder constant (0.95). Do not assimilate ρhv with TCWA2.

---

## Common Issues and Solutions

### 1. Missing GeoVaLs
**Error**: "Variable xxx not found in GeoVaLs"
**Solution**: Ensure your model output includes all required variables for your chosen scheme (see table above)

### 2. Wrong Operator/Scheme Combination
**Error**: "TCWA2 operator requires qi, ni, qc, smlf, gmlf, nr, ns, ng"
**Solution**: 
- If using TCWA2 microphysics → Must set `polarimetric operator: TCWA2`
- If using Thompson/WSM6/NSSL → Use `polarimetric operator: Zhang21` (default)

### 3. Missing Melted Fraction Fields
**Error**: Variable not found for smlf/gmlf
**Solution**: 
- These are required for TCWA2 only
- Other schemes (Thompson/WSM6/NSSL) don't need these fields

---

## Complete Example: TCWA2 with All Settings

```yaml
observations:
  observers:
  - obs space:
      name: Radar
      obsdatain:
        engine:
          type: H5File
          obsfile: radar_obs.nc4
      simulated variables: [equivalentReflectivityFactor, differentialReflectivity, 
                            specificDifferentialPhase, copolarCorrelationCoefficient]
    obs operator:
      name: PPRO
      polarimetric operator: TCWA2
      microphysics option: TCWA2
      VertCoord: height_above_mean_sea_level
      # Variable mappings
      var_rain_mixing_ratio: qr
      var_snow_mixing_ratio: qs
      var_graupel_mixing_ratio: qg
      var_cloud_mixing_ratio: qc
      var_ice_mixing_ratio: qi
      var_rain_number_concentration: nr
      var_snow_number_concentration: ns
      var_graupel_number_concentration: ng
      var_ice_number_concentration: ni
      var_snow_melted_fraction: smlf
      var_graupel_melted_fraction: gmlf
    obs error:
      covariance model: diagonal
```

---

## References

1. **Zhang21 Operator**:  
   Zhang, G., J. Gao, and M. Du, 2021: Parameterized forward operators for simulation and assimilation of polarimetric radar data with numerical weather predictions. *Adv. Atmos. Sci.*, **38**(5), 737−754.

2. **Zhang21 Melting Scheme**:  
   Liu, et al., 2024: A New Melting Model and Its Implementation in Parameterized Forward Operators for Polarimetric Radar Data Simulation With Double Moment Microphysics Schemes. *JGR Atmospheres*, **129**.

3. **Application of Zhang21 Operator in MPAS-JEDI**:  
   Kong et al., 2026: Assimilation of Radar Radial Velocity and Polarimetric Observations Using LETKF within MPAS-JEDI: Development of Assimilation Capabilities and Test with an Afternoon Thunderstorm in Taiwan. *Mon. Wea. Rev.*, submitted.

4. **TCWA2 Operator**:  
   Tsai, T.-C., J.-P. Chen, Z. Liu, S.-Y. Jiang, R. Kong, Y.-J. Wu, J. Ban, L.-F. Hsiao, Y.-S. Tang, P.-L. Chang, and J.-S. Hong, 2026: Development of the TCWA2 Bulk Cloud Microphysics Scheme and Its Integration with a Dual-Polarization Radar Operator for Forecasting Applications. *J. Adv. Model. Earth Syst.*, submitted.

5. **PPRO Library Documentation**:  
   See `README.md` in this directory for library-level documentation.

---

## Version History

- **2026-06**: Full C++ rewrite of the operator from Fortran (Zhang21 and TCWA2); coefficients embedded in source, no external library/files
- **2026-04**: Integrated back into UFO as a built-in module (removed external library dependency)
- **2025-01**: Added TCWA2 support
- **2024**: Initial PPRO implementation with Zhang21, Thompson, WSM6, NSSL support

