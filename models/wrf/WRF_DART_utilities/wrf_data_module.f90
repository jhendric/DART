MODULE wrf_data_module

!  $Source$
!  $Revision$
!  $Date$

TYPE wrf_data

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

   integer :: ncid  ! netcdf id for file
   integer :: bt_id, bt, sn_id, sn, we_id, we
   integer :: u_id, v_id, w_id, ph_id, phb_id, t_id,   &
                        mu_id, mub_id,                           &
                        qv_id, qc_id, qr_id, qi_id, qs_id, qg_id 
   integer :: ptop_id
   logical :: ice_micro


!---
!  arrays for data

   real, pointer :: u(:,:,:)
   real, pointer :: v(:,:,:)
   real, pointer :: w(:,:,:)
   real, pointer :: ph(:,:,:)
   real, pointer :: phb(:,:,:)
   real, pointer :: t(:,:,:)
   real, pointer :: qv(:,:,:)
   real, pointer :: qc(:,:,:)
   real, pointer :: qr(:,:,:)
   real, pointer :: qi(:,:,:)
   real, pointer :: qs(:,:,:)
   real, pointer :: qg(:,:,:)
   real, pointer :: mu(:,:)
   real, pointer :: mub(:,:)

end type
END MODULE wrf_data_module

