MODULE MO_BAND_WAVENUMBERS
  USE mo_rrtmgp_kind

  IMPLICIT NONE

  ! INTEGER/REAL PARAMETER DECLARATIONS FOR BUILDING
  INTEGER, PARAMETER :: NBAND_LW = 16, NBAND_SW = 14

  ! BAND PARAMETERS
  REAL(WP), DIMENSION(2, NBAND_LW), PARAMETER :: BND_LIMITS_WAVENUMBER_LW = &
       RESHAPE((/ &
          10._wp,  250.0_wp, &
        250.0_wp,  500.0_wp, &
        500.0_wp,  630.0_wp, &
        630.0_wp,  700.0_wp, &
        700.0_wp,  820.0_wp, &
        820.0_wp,  980.0_wp, &
        980.0_wp, 1080.0_wp, &
       1080.0_wp, 1180.0_wp, &
       1180.0_wp, 1390.0_wp, &
       1390.0_wp, 1480.0_wp, &
       1480.0_wp, 1800.0_wp, &
       1800.0_wp, 2080.0_wp, &
       2080.0_wp, 2250.0_wp, &
       2250.0_wp, 2390.0_wp, &
       2390.0_wp, 2680.0_wp, &
       2680.0_wp, 3250._wp/), &
       (/2,NBAND_LW/))

  REAL(WP), DIMENSION(2, NBAND_SW), PARAMETER :: BND_LIMITS_WAVENUMBER_SW = &
       RESHAPE((/ &
        820.0_wp, 2600.0_wp, &
       2600.0_wp, 3250.0_wp, &
       3250.0_wp, 4000.0_wp, &
       4000.0_wp, 4650.0_wp, &
       4650.0_wp, 5150.0_wp, &
       5150.0_wp, 6150.0_wp, &
       6150.0_wp, 7700.0_wp, &
       7700.0_wp, 8050.0_wp, &
       8050.0_wp,12850.0_wp, &
      12850.0_wp,16000.0_wp, &
      16000.0_wp,22650.0_wp, &
      22650.0_wp,29000.0_wp, &
      29000.0_wp,38000.0_wp, &
      38000.0_wp,50000.0_wp/), &
       (/2,NBAND_SW/))

end MODULE MO_BAND_WAVENUMBERS
