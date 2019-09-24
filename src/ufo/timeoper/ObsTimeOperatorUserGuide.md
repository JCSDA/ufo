# Obs Time Operator User Guide

### Yaml usage

For now the yaml should look like (taken from the standalone ufo test file):

```yaml
window_begin: 2018-04-14T21:00:00Z
window_end: 2018-04-15T03:00:00Z

Observations:
  ObsTypes:
  - ObsOperator:
      name: TimeOpInterp
      lin_name: VertInterp    
      windowSub: PT3H    
      ObsOperator:
        name: VertInterp
    ObsSpace:
      name: Radiosonde
      ObsDataIn:
        obsfile: Data/sondes_obs_2018041500_vs.nc4
      simulate:
        variables: [air_temperature]
    GeoVaLs:
      filename: Data/sondes_geoval_2018041500_vs.nc4
      loc_multiplier: 2
    vecequiv: GsiHofX
    tolerance: 1.0e-06
```

* **lin_name** - is needed to find the appropriate factory name for the LinearObsOperator. This is not need when the 
tlad code is added (and  then can be removed)
* **name** - is the key to the time operator factory to instantiate the code in the TimeOperator directory.
* **windowSub** - this is needed to populate the time-weights and calculate stateTime.  By doing it this way we can use
 a similar approach in TLAD to populate the time weights.
* **operator_type** - to allow different options to the timeoperator (such as accumulation of observations).
### Technical overview

We have following files

* **ObsTimeOper.h** and **ObsTimeOper.cc** This does the actual wrapping of one ObsOperator into another.  It is using
public inheritance from **ObsOperatorBase**.
* The geovals are created to be *timeStencil_* times its normal size and the data is ordered with each observation
having *timeStencil_* positions that are continguous.
* I have hardwired *timeStencil_* to 2 in the C++ constructor of ObsTimeOper.
* The **time_weight** is currently calculated in the Fortran locations module in **ufo_timeoper_set_timeweight** in 
**ufo_timeoepr_mod.F90**.
* The ObsOperator name is stored in the Fortran module **time_oper_mod**. This is done so that in the future other time 
weighting functions can be created and used dependent on the ObsOperator name.
* The **simulateObs** part is generic for all possible time operators. Ideally I would like to resize the geovals
object back to its normal size at this point and avoid the missing data indicator issue.

### How to extend
* Add an additional factory name to **ObsTimeOper.cc**
* Go to **ufo_timeoper_F90** and add the new time weighting function. 

### Testing
* I have removed print statements for this PR.
* Note that the values of the states are the same at different times in this contrived test.
* I have added some temporary code into **oops/test/interface/ObsOperator.h** to test the locations
 object generated in time oper.




