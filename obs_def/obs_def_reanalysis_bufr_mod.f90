! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! BEGIN DART PREPROCESS KIND LIST
!RADIOSONDE_U_WIND_COMPONENT,  KIND_U_WIND_COMPONENT
!RADIOSONDE_V_WIND_COMPONENT,  KIND_V_WIND_COMPONENT
!RADIOSONDE_SURFACE_PRESSURE,  KIND_SURFACE_PRESSURE
!RADIOSONDE_TEMPERATURE,       KIND_TEMPERATURE
!RADIOSONDE_SPECIFIC_HUMIDITY, KIND_SPECIFIC_HUMIDITY
!AIRCRAFT_U_WIND_COMPONENT,    KIND_U_WIND_COMPONENT
!AIRCRAFT_V_WIND_COMPONENT,    KIND_V_WIND_COMPONENT
!AIRCRAFT_TEMPERATURE,         KIND_TEMPERATURE
!AIRCRAFT_SPECIFIC_HUMIDITY,   KIND_SPECIFIC_HUMIDITY
!ACARS_U_WIND_COMPONENT,       KIND_U_WIND_COMPONENT
!ACARS_V_WIND_COMPONENT,       KIND_V_WIND_COMPONENT
!ACARS_TEMPERATURE,            KIND_TEMPERATURE
!ACARS_SPECIFIC_HUMIDITY,      KIND_SPECIFIC_HUMIDITY
!MARINE_SFC_U_WIND_COMPONENT,  KIND_U_WIND_COMPONENT
!MARINE_SFC_V_WIND_COMPONENT,  KIND_V_WIND_COMPONENT
!MARINE_SFC_TEMPERATURE,       KIND_TEMPERATURE
!MARINE_SFC_SPECIFIC_HUMIDITY, KIND_SPECIFIC_HUMIDITY
!LAND_SFC_U_WIND_COMPONENT,    KIND_U_WIND_COMPONENT
!LAND_SFC_V_WIND_COMPONENT,    KIND_V_WIND_COMPONENT
!LAND_SFC_TEMPERATURE,         KIND_TEMPERATURE
!LAND_SFC_SPECIFIC_HUMIDITY,   KIND_SPECIFIC_HUMIDITY
!SAT_U_WIND_COMPONENT,         KIND_U_WIND_COMPONENT
!SAT_V_WIND_COMPONENT,         KIND_V_WIND_COMPONENT
!ATOV_TEMPERATURE,             KIND_TEMPERATURE
!AIRS_TEMPERATURE,             KIND_TEMPERATURE
!AIRS_SPECIFIC_HUMIDITY,       KIND_SPECIFIC_HUMIDITY
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!! No use statements are required for the reanalysis bufr obs_def module
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(RADIOSONDE_U_WIND_COMPONENT, AIRCRAFT_U_WIND_COMPONENT, &
!              ACARS_U_WIND_COMPONENT, MARINE_SFC_U_WIND_COMPONENT, &
!              SAT_U_WIND_COMPONENT, LAND_SFC_U_WIND_COMPONENT)
!            call interpolate(state, location, KIND_U_WIND_COMPONENT, obs_val, istatus)         
!         case(RADIOSONDE_V_WIND_COMPONENT, AIRCRAFT_V_WIND_COMPONENT, &
!              ACARS_V_WIND_COMPONENT, MARINE_SFC_V_WIND_COMPONENT, &
!              SAT_V_WIND_COMPONENT, LAND_SFC_V_WIND_COMPONENT)
!            call interpolate(state, location, KIND_V_WIND_COMPONENT, obs_val, istatus)         
!         case(RADIOSONDE_TEMPERATURE, AIRCRAFT_TEMPERATURE, ACARS_TEMPERATURE, &
!              MARINE_SFC_TEMPERATURE, LAND_SFC_TEMPERATURE, ATOV_TEMPERATURE, &
!              AIRS_TEMPERATURE )
!            call interpolate(state, location, KIND_TEMPERATURE, obs_val, istatus)
!         case(RADIOSONDE_SPECIFIC_HUMIDITY, AIRCRAFT_SPECIFIC_HUMIDITY, & 
!              ACARS_SPECIFIC_HUMIDITY, MARINE_SFC_SPECIFIC_HUMIDITY, &
!              LAND_SFC_SPECIFIC_HUMIDITY, AIRS_SPECIFIC_HUMIDITY)
!            call interpolate(state, location, KIND_SPECIFIC_HUMIDITY, obs_val, istatus)
!            !!! UNITS in original BUFR are g/kg; This is converted to kg/kg by
!            !!! the BUFR to obs_sequence conversion programs making the line below
!            !!! unnecessary. PLEASE pay attention to units for specific humidity in models.
!            !!!obs_val = obs_val * 1000.0_r8
!         case(RADIOSONDE_SURFACE_PRESSURE)
!            call interpolate(state, location, KIND_SURFACE_PRESSURE, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!case(RADIOSONDE_U_WIND_COMPONENT, RADIOSONDE_V_WIND_COMPONENT, RADIOSONDE_SURFACE_PRESSURE, &
!   RADIOSONDE_TEMPERATURE, RADIOSONDE_SPECIFIC_HUMIDITY, &
!   AIRCRAFT_U_WIND_COMPONENT, AIRCRAFT_V_WIND_COMPONENT, &
!   AIRCRAFT_TEMPERATURE, AIRCRAFT_SPECIFIC_HUMIDITY, ACARS_U_WIND_COMPONENT, & 
!   ACARS_V_WIND_COMPONENT, ACARS_TEMPERATURE, ACARS_SPECIFIC_HUMIDITY, &
!   MARINE_SFC_U_WIND_COMPONENT, MARINE_SFC_V_WIND_COMPONENT, MARINE_SFC_TEMPERATURE, & 
!   MARINE_SFC_SPECIFIC_HUMIDITY,  LAND_SFC_U_WIND_COMPONENT, LAND_SFC_V_WIND_COMPONENT, & 
!   LAND_SFC_TEMPERATURE, LAND_SFC_SPECIFIC_HUMIDITY, SAT_U_WIND_COMPONENT, &
!   SAT_V_WIND_COMPONENT, ATOV_TEMPERATURE, AIRS_TEMPERATURE, AIRS_SPECIFIC_HUMIDITY )
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!case(RADIOSONDE_U_WIND_COMPONENT, RADIOSONDE_V_WIND_COMPONENT, RADIOSONDE_SURFACE_PRESSURE, &
!   RADIOSONDE_TEMPERATURE, RADIOSONDE_SPECIFIC_HUMIDITY, &
!   AIRCRAFT_U_WIND_COMPONENT, AIRCRAFT_V_WIND_COMPONENT, &
!   AIRCRAFT_TEMPERATURE, AIRCRAFT_SPECIFIC_HUMIDITY, ACARS_U_WIND_COMPONENT, &
!   ACARS_V_WIND_COMPONENT, ACARS_TEMPERATURE, ACARS_SPECIFIC_HUMIDITY, &
!   MARINE_SFC_U_WIND_COMPONENT, MARINE_SFC_V_WIND_COMPONENT, MARINE_SFC_TEMPERATURE, &
!   MARINE_SFC_SPECIFIC_HUMIDITY,  LAND_SFC_U_WIND_COMPONENT, LAND_SFC_V_WIND_COMPONENT, &
!   LAND_SFC_TEMPERATURE, LAND_SFC_SPECIFIC_HUMIDITY, SAT_U_WIND_COMPONENT, &
!   SAT_V_WIND_COMPONENT, ATOV_TEMPERATURE, AIRS_TEMPERATURE, AIRS_SPECIFIC_HUMIDITY )
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!case(RADIOSONDE_U_WIND_COMPONENT, RADIOSONDE_V_WIND_COMPONENT, RADIOSONDE_SURFACE_PRESSURE, &
!   RADIOSONDE_TEMPERATURE, RADIOSONDE_SPECIFIC_HUMIDITY, &
!   AIRCRAFT_U_WIND_COMPONENT, AIRCRAFT_V_WIND_COMPONENT, &
!   AIRCRAFT_TEMPERATURE, AIRCRAFT_SPECIFIC_HUMIDITY, ACARS_U_WIND_COMPONENT, &
!   ACARS_V_WIND_COMPONENT, ACARS_TEMPERATURE, ACARS_SPECIFIC_HUMIDITY, &
!   MARINE_SFC_U_WIND_COMPONENT, MARINE_SFC_V_WIND_COMPONENT, MARINE_SFC_TEMPERATURE, &
!   MARINE_SFC_SPECIFIC_HUMIDITY,  LAND_SFC_U_WIND_COMPONENT, LAND_SFC_V_WIND_COMPONENT, &
!   LAND_SFC_TEMPERATURE, LAND_SFC_SPECIFIC_HUMIDITY, SAT_U_WIND_COMPONENT, &
!   SAT_V_WIND_COMPONENT, ATOV_TEMPERATURE, AIRS_TEMPERATURE, AIRS_SPECIFIC_HUMIDITY )
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

