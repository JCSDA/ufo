!
! Observation error namelist for GPSRO observations
! This entry is for satellite 3 (known as Metop-B (DMI)), taken from processing centre 94
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 4 (known as Metop-A (DMI)), taken from processing centre 94
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 5 (known as Metop-C (DMI)), taken from processing centre 94
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 3 (known as Metop-B (EUMETSAT)), taken from processing centre 254
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 4 (known as Metop-A (EUMETSAT)), taken from processing centre 254
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 5 (known as Metop-C (EUMETSAT)), taken from processing centre 254
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 42 (known as TerraSAR-X (GFZ)), taken from processing centre 78
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 43 (known as TanDEM-X (GFZ)), taken from processing centre 78
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 522 (known as FY-3C (CMA)), taken from processing centre 38
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 523 (known as FY-3D (CMA)), taken from processing centre 38
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 740 (known as COSMIC-1 FM1 (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 741 (known as COSMIC-1 FM2(UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 742 (known as COSMIC-1 FM3 (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 743 (known as COSMIC-1 FM4 (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 744 (known as COSMIC-1 FM5 (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 745 (known as COSMIC-1 FM6 (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 825 (known as KOMPSAT-5 (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 269 (known as Spire constellation), taken from processing centre 178
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 750 (known as COSMIC-2 E1 (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 751 (known as COSMIC-2 E2 (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 752 (known as COSMIC-2 E3 (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 753 (known as COSMIC-2 E4 (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 754 (known as COSMIC-2 E5 (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 755 (known as COSMIC-2 E6 (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
!
! Observation error namelist for GPSRO observations
! This entry is for satellite 44 (known as PAZ (UCAR)), taken from processing centre 60
!
! Each entry is for a given range of average troposphere temperatures, which will be interpolated between.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background bending angle).
!
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3e-06
  av_temp = 225
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.048, 0.020, 0.008, 0.016, 0.300
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3e-06
  av_temp = 230
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.058, 0.014, 0.010, 0.037, 0.280
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3e-06
  av_temp = 235
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.064, 0.029, 0.012, 0.015, 0.220
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3e-06
  av_temp = 240
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.075, 0.033, 0.014, 0.014, 0.175
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3e-06
  av_temp = 245
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.081, 0.036, 0.012, 0.017, 0.146
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3e-06
  av_temp = 250
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.111, 0.048, 0.012, 0.020, 0.160
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3e-06
  av_temp = 255
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.132, 0.061, 0.012, 0.023, 0.152
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3e-06
  av_temp = 260
  max_height = 20000
  heights = 3000.0, 10000.0, 20000.0, 40000.0, 59900.0
  obs_errors = 0.128, 0.072, 0.010, 0.019, 0.179
/
