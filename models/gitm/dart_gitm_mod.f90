! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module dart_gitm_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! This is the interface between the GITM modules and DART.
! To reduce the possibility of scoping issues, all the
! unrestricted GITM modules are confined to this module.

use ModConstants
use ModKind 
use ModTime 
use ModSizeGitm
use ModPlanet 

use typesizes
use netcdf

use    utilities_mod, only : error_handler, E_ERR, E_WARN, E_MSG

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_gitm_nLons,    &
          get_gitm_nLats,    &
          get_gitm_nAlts,    &
          get_nSpecies,      &
          get_nSpeciesTotal, &
          get_nIons,         &
          get_nSpeciesAll,   &
          decode_gitm_indices

! version controlled file description for error handling, do not edit

character(len=128), parameter :: &
   source   = '$URL$', &
   revision = '$Revision$', &
   revdate  = '$Date$'

character(len=256) :: string1, string2

contains

!===================================================================
! All the public interfaces ... nothing more.
!===================================================================

integer function get_gitm_nLons()
   get_gitm_nLons = nLons
end function get_gitm_nLons

integer function get_gitm_nLats()
   get_gitm_nLats = nLats
end function get_gitm_nLats

integer function get_gitm_nAlts()
   get_gitm_nAlts = nAlts
end function get_gitm_nAlts

integer function get_nSpecies()
   get_nSpecies = nSpecies   ! From ModPlanet, hopefully
end function get_nSpecies

integer function get_nSpeciesTotal()
   get_nSpeciesTotal = nSpeciesTotal   ! From ModPlanet, hopefully
end function get_nSpeciesTotal

integer function get_nIons()
   get_nIons = nIons   ! From ModPlanet, hopefully
end function get_nIons

integer function get_nSpeciesAll()
   get_nSpeciesAll = nSpeciesAll   ! From ModPlanet, hopefully
end function get_nSpeciesAll


subroutine decode_gitm_indices( state_variables, gitm_varname, gitm_dim, gitm_index, &
               long_name, units, varname, kind_string, dart_kind)
! The rosetta stone relating the user input 'strings' to integer indices. 
!
! progvar%varname      = varname
! progvar%long_name    = long_name
! progvar%units        = units
! progvar%gitm_varname = gitm_varname
! progvar%gitm_dim     = gitm_dim
! progvar%gitm_index   = gitm_index

character(len=*), dimension(:),   intent(in)  :: state_variables
character(len=*),                 intent(out) :: gitm_varname
integer,                          intent(out) :: gitm_dim, gitm_index
character(len=NF90_MAX_NAME),     intent(out) :: long_name
character(len=NF90_MAX_NAME),     intent(out) :: units
character(len=NF90_MAX_NAME),     intent(out) :: varname, kind_string
integer,                          intent(out) :: dart_kind

integer :: ngood, nrows, i, varid


ngood = 0
MyLoop : do i = 1, nrows

   varname = trim(state_variables(2*i -1))
   kind_string = trim(state_variables(2*i   ))

   if ( varname == ' ' .and. kind_string == ' ' ) exit MyLoop ! Found end of list.

   if ( varname == ' ' .or. kind_string == ' ' ) then
      string1 = 'gitm_vars_nml:gitm state_variables not fully specified'
      call error_handler(E_ERR,'decode_gitm_indices',string1,source,revision,revdate)
   endif

   long_name    = 'something real'
   units        = 'furlongs/fortnight'

   select case (trim(varname))

   ! The first hunk of these all come from the NDensityS variable, defined to be:
   ! do iSpecies=1,nSpeciesTotal
   !    write(iRestartUnit_) NDensityS(:,:,:,iSpecies,iBlock)
   ! endd

   case ('iO_3P_NDensityS') 
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index   = iO_3P_
      long_name    = 'something surreal'
      units        = 'furlongs/fortnight'

   case ('iO2_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iO2_

   case ('iN2_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iN2_

   case ('iN_4S_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iN_4S_

   case ('iNO_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iNO_

   case ('iN_2D_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iN_2D_

   case ('iN_2P_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iN_2P_

   case ('iH_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iH_

   case ('iHe_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iHe_

   case ('iAr_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iAr_

   case ('iO_1D_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iO_1D_

   ! The next hunk of these all pertain to the IDensityS variable:
   ! do iSpecies=1,nIons
   !    write(iRestartUnit_) IDensityS(:,:,:,iSpecies,iBlock)
   ! enddo

   case ('iO_4SP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iO_4SP_

   case ('iO2P_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iO2P_

   case ('iN2P_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iN2P_

   case ('iNP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iNP_

   case ('iNOP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iNOP_

   case ('iO_2DP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iO_2DP_

   case ('iO_2PP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iO_2PP_

   case ('iHP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iHP_

   case ('iHeP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iHeP_

   case ('ie_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = ie_

   case ('Temperature') ! write(iRestartUnit_)  Temperature(:,:,:,iBlock)
      gitm_varname = 'Temperature'
      gitm_dim     = -1
      gitm_index   = -1

   case ('ITemperature') ! write(iRestartUnit_) ITemperature(:,:,:,iBlock)
      gitm_varname = 'ITemperature'
      gitm_dim     = -1
      gitm_index   = -1

   case ('eTemperature') ! write(iRestartUnit_) eTemperature(:,:,:,iBlock)
      gitm_varname = 'eTemperature'
      gitm_dim     = -1
      gitm_index   = -1

   case ('U_Velocity_component') ! write(iRestartUnit_) Velocity(:,:,:,iBlock)
      gitm_varname = 'Velocity'
      gitm_dim     = 4
      gitm_index   = 1

   case ('V_Velocity_component') ! write(iRestartUnit_) Velocity(:,:,:,iBlock)
      gitm_varname = 'Velocity'
      gitm_dim     = 4
      gitm_index   = 2

   case ('W_Velocity_component') ! write(iRestartUnit_) Velocity(:,:,:,iBlock)
      gitm_varname = 'Velocity'
      gitm_dim     = 4
      gitm_index   = 3

   case ('U_IVelocity_component') ! write(iRestartUnit_) Velocity(:,:,:,iBlock)
      gitm_varname = 'IVelocity'
      gitm_dim     = 4
      gitm_index   = 1

   case ('V_IVelocity_component') ! write(iRestartUnit_) Velocity(:,:,:,iBlock)
      gitm_varname = 'IVelocity'
      gitm_dim     = 4
      gitm_index   = 2

   case ('W_IVelocity_component') ! write(iRestartUnit_) IVelocity(:,:,:,iBlock)
      gitm_varname = 'IVelocity'
      gitm_dim     = 4
      gitm_index   = 3

   case ('iO_3P_VerticalVelocity')
      gitm_varname = 'VerticalVelocity'
      gitm_dim     = 4
      gitm_index   = iO_3P_

   case ('iO2_VerticalVelocity')
      gitm_varname = 'VerticalVelocity'
      gitm_dim     = 4
      gitm_index   = iO2_

   case ('iN2_VerticalVelocity')
      gitm_varname = 'VerticalVelocity'
      gitm_dim     = 4
      gitm_index   = iN2_

   case ('iN_4S_VerticalVelocity')
      gitm_varname = 'VerticalVelocity'
      gitm_dim     = 4
      gitm_index   = iN_4S_

   case ('iNO_VerticalVelocity')
      gitm_varname = 'VerticalVelocity'
      gitm_dim     = 4
      gitm_index   = iNO_

   case default

      write(string1,*)'unknown GITM variable '//trim(varname)
      call error_handler(E_ERR,'define_var_dims',string1,source,revision,revdate)

   end select

   ! Make sure DART kind is valid

!  dart_kind = get_raw_obs_kind_index(kind_string)
!  if( dart_kind < 0 ) then
!     write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(kind_string)
!     call error_handler(E_ERR,'decode_gitm_indices',string1,source,revision,revdate)
!  endif

   ! Record the contents of the DART state vector 

!  if ( debug > 0 ) then
 !    write(logfileunit,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
 !    write(     *     ,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
!  endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'decode_gitm_indices',string1,source,revision,revdate,text2=string2)
endif

end subroutine decode_gitm_indices




!===================================================================
! End of dart_gitm_mod
!===================================================================
end module dart_gitm_mod
