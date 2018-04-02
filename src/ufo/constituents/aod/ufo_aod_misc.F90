MODULE ufo_aod_misc

  USE netcdf

  PRIVATE
  PUBLIC :: max_name_length,max_string_length,max_vars
  PUBLIC :: small_value

  REAL, PARAMETER :: small_value=1.e-16

  INTEGER, PARAMETER :: &
       &max_string_length=NF90_MAX_NAME,&
       &max_name_length=56,&
       &max_vars=50
 
END MODULE ufo_aod_misc
