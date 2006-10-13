! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module assim_tools_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 
!
! A variety of operations required by assimilation.

use      types_mod, only : r8
use  utilities_mod,       only : file_exist, get_unit, check_namelist_read, do_output, &
                                 find_namelist_in_file, register_module, error_handler, &
                                 E_ERR, E_MSG, logfileunit
use       sort_mod,       only : index_sort 
use random_seq_mod,       only : random_seq_type, random_gaussian, init_random_seq, &
                                 random_uniform

use obs_sequence_mod,     only : obs_sequence_type, obs_type, get_num_copies, get_num_qc, &
                                 init_obs, get_obs_from_key, get_obs_def, get_obs_values, &
                                 destroy_obs
   
use          obs_def_mod, only : obs_def_type, get_obs_def_location, get_obs_def_time, &
                                 get_obs_def_error_variance, get_obs_kind

use       cov_cutoff_mod, only : comp_cov_factor

use       reg_factor_mod, only : comp_reg_factor

use         location_mod, only : location_type, get_close_maxdist_init, &
                                 get_close_obs_init, get_close_obs, get_close_type

use ensemble_manager_mod, only : ensemble_type, get_my_num_vars, get_my_vars, &
                                 compute_copy_mean_var, get_var_owner_index

use mpi_utilities_mod,    only : my_task_id, broadcast_send, broadcast_recv, &
                                 sum_across_tasks

use adaptive_inflate_mod, only : do_obs_inflate,  do_single_ss_inflate,           &
                                 do_varying_ss_inflate, get_inflate, set_inflate, &
                                 get_sd, set_sd, update_inflation,                &
                                 inflate_ens, adaptive_inflate_type, deterministic_inflate

use time_manager_mod,     only : time_type

use assim_model_mod,      only : get_state_meta_data

implicit none
private

public :: filter_assim

! Indicates if module initialization subroutine has been called yet
logical                :: module_initialized = .false.

type (random_seq_type) :: inc_ran_seq
logical                :: first_inc_ran_call = .true.
character(len = 129)   :: errstring

! Need to read in table for off-line based sampling correction and store it
logical                :: first_get_correction = .true.
real(r8)               :: exp_true_correl(200), alpha(200)                                                                      
! CVS Generated file description for error handling, do not edit
character(len=128), parameter :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

!============================================================================

!---- namelist with default values

! Filter kind selects type of observation space filter
!      1 = EAKF filter
!      2 = ENKF
!      3 = Kernel filter
!      4 = particle filter
!      5 = random draw from posterior
!      6 = deterministic draw from posterior with fixed kurtosis
integer  :: filter_kind                     = 1
real(r8) :: cutoff                          = 0.2_r8
logical  :: sort_obs_inc                    = .false.
logical  :: spread_restoration              = .false.
logical  :: sampling_error_correction       = .false.
integer  :: adaptive_localization_threshold = -1

namelist / assim_tools_nml / filter_kind, cutoff, sort_obs_inc, &
   spread_restoration, sampling_error_correction, adaptive_localization_threshold

!============================================================================

contains

!-------------------------------------------------------------

subroutine assim_tools_init()

integer :: iunit, io

call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "assim_tools_nml", iunit)
read(iunit, nml = assim_tools_nml, iostat = io)
call check_namelist_read(iunit, io, "assim_tools_nml")

! Write the namelist values to the log file
call error_handler(E_MSG,'assim_tools_init','assim_tools namelist values',' ',' ',' ')
if (do_output()) write(logfileunit, nml=assim_tools_nml)
if (do_output()) write(     *     , nml=assim_tools_nml)

! FOR NOW, can only do spread restoration with filter option 1 (need to extend this)
if(spread_restoration .and. .not. filter_kind == 1) then
   write(errstring, *) 'cant combine spread_restoration and filter_kind ', filter_kind
   call error_handler(E_ERR,'assim_tools_init', errstring, source, revision, revdate)
endif

end subroutine assim_tools_init

!-------------------------------------------------------------

subroutine filter_assim(ens_handle, obs_ens_handle, obs_seq, keys,           &
   ens_size, num_groups, obs_val_index, inflate, ENS_MEAN_COPY, ENS_SD_COPY, &
   ENS_INF_COPY, ENS_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY,          &
   OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START,            &
   OBS_PRIOR_VAR_END, inflate_only)

type(ensemble_type),         intent(inout) :: ens_handle, obs_ens_handle
type(obs_sequence_type),     intent(in)    :: obs_seq
integer,                     intent(in)    :: keys(:)
integer,                     intent(in)    :: ens_size, num_groups, obs_val_index
type(adaptive_inflate_type), intent(inout) :: inflate
integer,                     intent(in)    :: ENS_MEAN_COPY, ENS_SD_COPY, ENS_INF_COPY
integer,                     intent(in)    :: ENS_INF_SD_COPY
integer,                     intent(in)    :: OBS_KEY_COPY, OBS_GLOBAL_QC_COPY
integer,                     intent(in)    :: OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END
integer,                     intent(in)    :: OBS_PRIOR_VAR_START, OBS_PRIOR_VAR_END
logical,                     intent(in)    :: inflate_only

real(r8) :: obs_prior(ens_size), obs_inc(ens_size), increment(ens_size)
real(r8) :: reg_factor
real(r8) :: net_a(num_groups), reg_coef(num_groups), correl(num_groups)
real(r8) :: cov_factor, obs(1), obs_err_var, my_inflate, my_inflate_sd
real(r8) :: varying_ss_inflate, varying_ss_inflate_sd
real(r8) :: ss_inflate_base, obs_qc, cutoff_rev
real(r8) :: gamma, ens_obs_mean, ens_obs_var, ens_var_deflate
real(r8) :: r_mean, r_var
real(r8) :: orig_obs_prior_mean(num_groups), orig_obs_prior_var(num_groups)
real(r8) :: obs_prior_mean(num_groups), obs_prior_var(num_groups)
real(r8) :: close_obs_dist(obs_ens_handle%my_num_vars)
real(r8) :: close_state_dist(ens_handle%my_num_vars)

integer  :: my_num_obs, i, j, owner, owners_index, my_num_state
integer  :: my_obs(obs_ens_handle%my_num_vars), my_state(ens_handle%my_num_vars)
integer  :: this_obs_key, obs_mean_index, obs_var_index
integer  :: grp_beg(num_groups), grp_end(num_groups), grp_size, grp_bot, grp_top, group
integer  :: close_obs_ind(obs_ens_handle%my_num_vars)
integer  :: close_state_ind(ens_handle%my_num_vars)
integer  :: num_close_obs, obs_index, num_close_states, state_index
integer  :: total_num_close_obs
integer  :: base_obs_kind, my_obs_kind(obs_ens_handle%my_num_vars)
integer  :: my_state_kind(ens_handle%my_num_vars)

type(location_type)  :: my_obs_loc(obs_ens_handle%my_num_vars), base_obs_loc
type(location_type)  :: my_state_loc(ens_handle%my_num_vars)
type(get_close_type) :: gc_obs, gc_state
type(obs_type)       :: observation
type(obs_def_type)   :: obs_def
type(time_type)      :: obs_time

! for performance, local copies 
logical :: local_single_ss_inflate
logical :: local_varying_ss_inflate
logical :: local_obs_inflate

!PAR: THIS SHOULD COME FROM SOMEWHERE ELSE AND BE NAMELIST CONTOLLED
real(r8) :: qc_threshold = 10.0_r8

! Initialize assim_tools_module if needed
if(.not. module_initialized) then
   call assim_tools_init
   module_initialized = .true.
endif

! For performance, make local copies of these settings which
! are really in the inflate derived type.
local_single_ss_inflate = do_single_ss_inflate(inflate)
local_varying_ss_inflate = do_varying_ss_inflate(inflate)
local_obs_inflate = do_obs_inflate(inflate)

! Divide ensemble into num_groups groups
grp_size = ens_size / num_groups
do group = 1, num_groups
   grp_beg(group) = (group - 1) * grp_size + 1
   grp_end(group) = grp_beg(group) + grp_size - 1
enddo

! Put initial value of state space inflation in copy normally used for SD
ens_handle%copies(ENS_SD_COPY, :) = ens_handle%copies(ENS_INF_COPY, :)

! For single state or obs space inflation, the inflation is like a token
! Gets passed from the processor with a given obs on to the next
if(local_single_ss_inflate) then
   my_inflate = ens_handle%copies(ENS_INF_COPY, 1)
   my_inflate_sd = ens_handle%copies(ENS_INF_SD_COPY, 1)
end if

! For obs space inflation, single value comes from storage in adaptive_inflate_mod
if(local_obs_inflate) then
   my_inflate = get_inflate(inflate)
   my_inflate_sd = get_sd(inflate)
endif

! Get info on my number and indices for obs
my_num_obs = get_my_num_vars(obs_ens_handle)
call get_my_vars(obs_ens_handle, my_obs)

! Construct an observation temporary
call init_obs(observation, get_num_copies(obs_seq), get_num_qc(obs_seq))
! Get the locations for all of my observations 
Locations: do i = 1, obs_ens_handle%my_num_vars
   this_obs_key = obs_ens_handle%copies(OBS_KEY_COPY, i)
   call get_obs_from_key(obs_seq, this_obs_key, observation)
   call get_obs_def(observation, obs_def)
   my_obs_loc(i) = get_obs_def_location(obs_def)
   my_obs_kind(i) = get_obs_kind(obs_def)
   ! Need the time for regression diagnostics potentially
   if(i == 1) obs_time = get_obs_def_time(obs_def)
end do Locations

! Get info on my number and indices for state
my_num_state = get_my_num_vars(ens_handle)
call get_my_vars(ens_handle, my_state)
! Get the location of all my state variables
do i = 1, ens_handle%my_num_vars
   call get_state_meta_data(my_state(i), my_state_loc(i), my_state_kind(i))
end do

! PAR: MIGHT BE BETTER TO HAVE ONE PE DEDICATED TO COMPUTING 
! INCREMENTS. OWNING PE WOULD SHIP IT'S PRIOR TO THIS ONE
! BEFORE EACH INCREMENT.

! Get mean and variance of each group's observation priors for adaptive inflation
! These are transmitted with the prior obs dist and increments
if(local_varying_ss_inflate .or. local_single_ss_inflate) then
   do group = 1, num_groups
      grp_bot = grp_beg(group)
      grp_top = grp_end(group)
      obs_mean_index = OBS_PRIOR_MEAN_START + group - 1
      obs_var_index  = OBS_PRIOR_VAR_START  + group - 1
         call compute_copy_mean_var(obs_ens_handle, grp_bot, grp_top, &
           obs_mean_index, obs_var_index) 
   end do
endif

! NOTE THESE COULD ONLY BE DONE ONCE PER RUN!!! FIGURE THIS OUT.
! The computations in the two get_close_maxdist_init are redundant
! Initialize the method for getting state variables close to a given ob on my process
call get_close_maxdist_init(gc_state, 2.0_r8*cutoff)
call get_close_obs_init(gc_state, my_num_state, my_state_loc)

! Initialize the method for getting obs close to a given ob on my process
call get_close_maxdist_init(gc_obs, 2.0_r8*cutoff)
call get_close_obs_init(gc_obs, my_num_obs, my_obs_loc)

! Loop through all the (global) observations sequentially
SEQUENTIAL_OBS: do i = 1, obs_ens_handle%num_vars

   ! Every pe has information about obs sequence
   call get_obs_from_key(obs_seq, keys(i), observation)
   call get_obs_def(observation, obs_def)
   base_obs_loc = get_obs_def_location(obs_def)
   obs_err_var = get_obs_def_error_variance(obs_def)
   base_obs_kind = get_obs_kind(obs_def)
   ! Get the value of the observation
   call get_obs_values(observation, obs, obs_val_index)

   ! Find out who has this observation and where it is
   call get_var_owner_index(i, owner, owners_index)
   ! Following block is done only by the owner of this observation
   !-----------------------------------------------------------------------
   if(my_task_id() == owner) then
      obs_qc = obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, owners_index)
      ! With new DART qc, only value of 0 for this field should be assimilated
      IF_QC_IS_OKAY: if(nint(obs_qc) ==0) then
         obs_prior = obs_ens_handle%copies(1:ens_size, owners_index)
         orig_obs_prior_mean = obs_ens_handle%copies(OBS_PRIOR_MEAN_START: &
            OBS_PRIOR_MEAN_END, owners_index)
         orig_obs_prior_var = obs_ens_handle%copies(OBS_PRIOR_VAR_START: &
            OBS_PRIOR_VAR_END, owners_index)

         ! Compute observation space increments for each group
         do group = 1, num_groups
            grp_bot = grp_beg(group)
            grp_top = grp_end(group)
            call obs_increment(obs_prior(grp_bot:grp_top), &
               grp_size, obs(1), obs_err_var, obs_inc(grp_bot:grp_top), &
               inflate, my_inflate, my_inflate_sd, net_a(group))
         end do

         ! Compute updated values for single state space inflation
         SINGLE_SS_INFLATE: if(local_single_ss_inflate) then
            ss_inflate_base = ens_handle%copies(ENS_SD_COPY, 1)
            do group = 1, num_groups
               if(my_inflate > 0.0_r8 .and. my_inflate_sd > 0.0_r8) then
                  ! For case with single spatial inflation, use gamma = 1.0_r8
                  gamma = 1.0_r8
                  ! Deflate the inflated variance
                  ens_obs_mean = orig_obs_prior_mean(group)
                  ens_obs_var = orig_obs_prior_var(group)
                  ens_var_deflate = ens_obs_var / &
                     (1.0_r8 + gamma*(sqrt(ss_inflate_base) - 1.0_r8))**2
                  
                  ! If this is inflate only (i.e. posterior) remove impact of this obs.
                  if(inflate_only) then
                     r_var = 1.0_r8 / (1.0_r8 / ens_var_deflate - 1.0_r8 / obs_err_var)
                     r_mean = r_var *(ens_obs_mean / ens_var_deflate - obs(1) / obs_err_var)
                  else
                     r_mean = ens_obs_mean
                     r_var = ens_var_deflate
                  endif

                  call update_inflation(inflate, my_inflate, my_inflate_sd, &
                     r_mean, r_var, obs(1), obs_err_var, &
                     gamma)
               endif
            end do
         endif SINGLE_SS_INFLATE
      endif IF_QC_IS_OKAY

      !Broadcast the info from this obs
      ! What gets broadcast depends on what kind of inflation is being done
      if(local_varying_ss_inflate) then
         call my_broadcast_send1(owner, obs_prior, obs_inc, orig_obs_prior_mean, &
            orig_obs_prior_var, obs_qc)
      else if(local_single_ss_inflate .or. local_obs_inflate) then
         call my_broadcast_send2(owner, obs_prior, obs_inc, my_inflate, &
            my_inflate_sd, obs_qc)
      else
         call my_broadcast_send3(owner, obs_prior, obs_inc, obs_qc)
      endif

   ! Next block is done by processes that do NOT own this observation
   !-----------------------------------------------------------------------
   else
      ! What gets broadcast depends on what kind of inflation is being done
      ! I don't store this obs; receive the obs prior and increment from broadcast
      if(local_varying_ss_inflate) then
         call my_broadcast_recv1(owner, obs_prior, obs_inc, orig_obs_prior_mean, &
            orig_obs_prior_var, obs_qc)
      else if(local_single_ss_inflate .or. local_obs_inflate) then
         call my_broadcast_recv2(owner, obs_prior, obs_inc, my_inflate, &
            my_inflate_sd, obs_qc)
      else
         call my_broadcast_recv3(owner, obs_prior, obs_inc, obs_qc)
      endif
   endif
   !-----------------------------------------------------------------------

   ! Everybody is doing this section, cycle if qc is bad
   if(nint(obs_qc) /= 0) cycle SEQUENTIAL_OBS

   ! Can compute prior mean and variance of obs for each group just once here
   do group = 1, num_groups
      grp_bot = grp_beg(group)
      grp_top = grp_end(group)
      obs_prior_mean(group) = sum(obs_prior(grp_bot:grp_top)) / grp_size
      obs_prior_var(group)  = sum(obs_prior(grp_bot:grp_top) * obs_prior(grp_bot:grp_top)) - &
         grp_size * obs_prior_mean(group)**2
   end do

   ! Now everybody updates their close states

   ! Getting close states for each scalar observation for now
   ! Need model state for some distance computations in sigma, have
   ! to pick some state, only current ones are ensemble, just pass 1

   ! PAR ALSO NEED A CALL TO SOMEPLACE THAT LETS PEOPLE MODIFY THE DISTANCES SOMEWHERE

   !------------------------------------------------------
   
   ! Find observations on my process that are close to one being assimilated
   ! Need to get obs density first in case of adaptive localization
   call get_close_obs(gc_obs, base_obs_loc, base_obs_kind, my_obs_loc, my_obs_kind, &
      num_close_obs, close_obs_ind, close_obs_dist)

   ! For adaptive localization, need number of other obs close to the chosen observation
   cutoff_rev = cutoff
   if(adaptive_localization_threshold > 0) then
      call sum_across_tasks(num_close_obs, total_num_close_obs)
      ! Want expected number of close observations to be reduced to some threshold
      if(total_num_close_obs > adaptive_localization_threshold) then
         ! Change the cutoff radius to get the appropriate number in the circle
         ! This is specific to models on sphere's
         ! Need to get thinning out of assim_tools and into something about locations
         ! Compute a new radius if the total_num_close is greater than the desired as
         ! 2*cutoff_rev = sqrt(2*cutoff * adaptive_localization_threshold / total_num_close_obs)
         cutoff_rev =  sqrt((2.0_r8*cutoff)**2 *adaptive_localization_threshold / &
            total_num_close_obs) / 2.0_r8
      endif
   endif

   call get_close_obs(gc_state, base_obs_loc, base_obs_kind, my_state_loc, my_state_kind, &
      num_close_states, close_state_ind, close_state_dist)

   STATE_UPDATE: do j = 1, num_close_states
      state_index = close_state_ind(j)

      ! Get the initial values of inflation for this variable if state varying inflation
      if(local_varying_ss_inflate) then
         varying_ss_inflate = ens_handle%copies(ENS_INF_COPY, state_index)
         varying_ss_inflate_sd = ens_handle%copies(ENS_INF_SD_COPY, state_index)
      else
         varying_ss_inflate    = 0.0_r8
         varying_ss_inflate_sd = 0.0_r8
      endif
     
      ! Compute the distance and covariance factor 
!PAR URGENT: MAKE SURE THIS INDEXING IS CORRECT; SAME FOR OTHER COMP_COV_FACTOR
      cov_factor = comp_cov_factor(close_state_dist(j), cutoff_rev, &
         base_obs_loc, base_obs_kind, my_state_loc(state_index), my_state_kind(state_index))
      if(cov_factor <= 0.0_r8) cycle STATE_UPDATE

      do group = 1, num_groups
         grp_bot = grp_beg(group)
         grp_top = grp_end(group)
         ! Do update of state, correl only needed for varying ss inflate
         if(local_varying_ss_inflate .and. varying_ss_inflate > 0.0_r8 .and. &
            varying_ss_inflate_sd > 0.0_r8) then
            call update_from_obs_inc(obs_prior(grp_bot:grp_top), obs_prior_mean(group), &
               obs_prior_var(group), obs_inc(grp_bot:grp_top), &
               ens_handle%copies(grp_bot:grp_top, state_index), grp_size, &
               increment(grp_bot:grp_top), reg_coef(group), net_a(group), correl(group))
         else
            call update_from_obs_inc(obs_prior(grp_bot:grp_top), obs_prior_mean(group), &
               obs_prior_var(group), obs_inc(grp_bot:grp_top), &
               ens_handle%copies(grp_bot:grp_top, state_index), grp_size, &
               increment(grp_bot:grp_top), reg_coef(group), net_a(group))
         endif
      end do

      ! Compute an information factor for impact of this observation on this state
      if(num_groups == 1) then
          reg_factor = 1.0_r8
      else
         ! Pass the time along with the index for possible diagnostic output
         ! Compute regression factor for this obs-state pair
         reg_factor = comp_reg_factor(num_groups, reg_coef, obs_time, i, my_state(state_index))
      endif

      reg_factor = min(reg_factor, cov_factor)

!PAR NEED TO TURN STUFF OFF MORE EFFICEINTLY
      if(.not. inflate_only) then
         ens_handle%copies(1:ens_size, state_index) = &
            ens_handle%copies(1:ens_size, state_index) + reg_factor * increment
      endif

      ! Compute spatially-variing state space inflation
      if(local_varying_ss_inflate) then
         ! base is the initial inflate value for this state variable
         ss_inflate_base = ens_handle%copies(ENS_SD_COPY, state_index)
         GroupInflate: do group = 1, num_groups
            if(varying_ss_inflate > 0.0_r8 .and. varying_ss_inflate_sd > 0.0_r8) then
               gamma = reg_factor * abs(correl(group))
               ! Deflate the inflated variance using the INITIAL state inflate
               ! value (before these obs started gumming it up).
               ens_obs_mean = orig_obs_prior_mean(group)
               ens_obs_var =  orig_obs_prior_var(group)

               ens_var_deflate = ens_obs_var / &
                  (1.0_r8 + gamma*(sqrt(ss_inflate_base) - 1.0_r8))**2
                  
                  ! If this is inflate only (i.e. posterior) remove impact of this obs.
                  if(inflate_only) then
                     r_var = 1.0_r8 / (1.0_r8 / ens_var_deflate - 1.0_r8 / obs_err_var)
                     r_mean = r_var *(ens_obs_mean / ens_var_deflate - obs(1) / obs_err_var)
                  else
                     r_mean = ens_obs_mean
                     r_var = ens_var_deflate
                  endif

               ! IS A TABLE LOOKUP POSSIBLE TO ACCELERATE THIS?
               call update_inflation(inflate, varying_ss_inflate, varying_ss_inflate_sd, &
                  r_mean, r_var, obs(1), obs_err_var, gamma)
            endif
            ! Copy updated values into ensemble obs storage
            ens_handle%copies(ENS_INF_COPY, state_index) = varying_ss_inflate
            ens_handle%copies(ENS_INF_SD_COPY, state_index) = varying_ss_inflate_sd
         end do GroupInflate
      endif

   end do STATE_UPDATE
   !------------------------------------------------------

   ! Now everybody updates their obs priors (only ones after this one)
   OBS_UPDATE: do j = 1, num_close_obs
      obs_index = close_obs_ind(j)
      ! Only have to update obs that have not yet been used
      if(my_obs(obs_index) > i) then

         ! Compute the distance and the covar_factor
         cov_factor = comp_cov_factor(close_obs_dist(j), cutoff_rev, &
            base_obs_loc, base_obs_kind, my_obs_loc(obs_index), my_obs_kind(obs_index))
         if(cov_factor <= 0.0_r8) cycle OBS_UPDATE

         do group = 1, num_groups
            grp_bot = grp_beg(group)
            grp_top = grp_end(group)
            call update_from_obs_inc(obs_prior(grp_bot:grp_top), obs_prior_mean(group), &
               obs_prior_var(group), obs_inc(grp_bot:grp_top), &
                obs_ens_handle%copies(grp_bot:grp_top, obs_index), grp_size, &
                increment(grp_bot:grp_top), reg_coef(group), net_a(group))
         end do

         ! Compute an information factor for impact of this observation on this state
         if(num_groups == 1) then
             reg_factor = 1.0_r8
         else
            ! Pass the time along with the index for possible diagnostic output
            ! Compute regression factor for this obs-state pair
            ! Negative indicates that this is an observation index
            reg_factor = comp_reg_factor(num_groups, reg_coef, obs_time, i, -1*my_obs(obs_index))
         endif

         reg_factor = min(reg_factor, cov_factor)

         if(.not. inflate_only) then
            obs_ens_handle%copies(1:ens_size, obs_index) = &
              obs_ens_handle%copies(1:ens_size, obs_index) + reg_factor * increment
         endif
      endif
   end do OBS_UPDATE
end do SEQUENTIAL_OBS

! Every pe needs to get the current my_inflate and my_inflate_sd back
if(local_single_ss_inflate) then
   ens_handle%copies(ENS_INF_COPY, :) = my_inflate
   ens_handle%copies(ENS_INF_SD_COPY, :) = my_inflate_sd
end if

! Everybody needs to store the latest value for obs_inflate
if(local_obs_inflate) then
   call set_inflate(inflate, my_inflate)
   call set_sd(inflate, my_inflate_sd)
endif

! Free up the storage
call destroy_obs(observation)

end subroutine filter_assim

!-------------------------------------------------------------

subroutine obs_increment(ens_in, ens_size, obs, obs_var, obs_inc, &
   inflate, my_cov_inflate, my_cov_inflate_sd, net_a)

! Given the ensemble prior for an observation, the observation, and
! the observation error variance, computes increments and adjusts
! observation space inflation values

integer,                     intent(in)    :: ens_size
real(r8),                    intent(in)    :: ens_in(ens_size), obs, obs_var
real(r8),                    intent(out)   :: obs_inc(ens_size)
type(adaptive_inflate_type), intent(inout) :: inflate
real(r8),                    intent(inout) :: my_cov_inflate, my_cov_inflate_sd
real(r8),                    intent(out)   :: net_a

real(r8) :: ens(ens_size), inflate_inc(ens_size)
real(r8) :: prior_mean, prior_var, new_val(ens_size)
integer  :: i, ens_index(ens_size), new_index(ens_size)

! Copy the input ensemble to something that can be modified
ens = ens_in

! Null value of net spread change factor is 1.0
net_a = 0.0_r8

! Compute prior variance and mean from sample
prior_mean = sum(ens) / ens_size
prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)

! If observation space inflation is being done, compute the initial 
! increments and update the inflation factor and its standard deviation
! as needed. my_cov_inflate < 0 means don't do any of this.
if(do_obs_inflate(inflate)) then
   ! If my_cov_inflate_sd is <= 0, just retain current my_cov_inflate setting
   if(my_cov_inflate_sd > 0.0_r8) & 
      ! Gamma set to 1.0 because no distance for observation space
      call update_inflation(inflate, my_cov_inflate, my_cov_inflate_sd, prior_mean, &
         prior_var, obs, obs_var, gamma = 1.0_r8)

   ! Now inflate the ensemble and compute a preliminary inflation increment
   call inflate_ens(inflate, ens, prior_mean, my_cov_inflate, prior_var)
   ! Keep the increment due to inflation alone 
   inflate_inc = ens - ens_in

   ! Need to recompute variance if non-deterministic inflation (mean is unchanged)
   if(.not. deterministic_inflate(inflate)) &
      prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)
endif

! If both obs_var and prior_var are 0 don't know what to do
if(obs_var == 0.0_r8 .and. prior_var == 0.0_r8) call error_handler(E_ERR,&
   'obs_increment', 'Both obs_var and prior_var are zero. This is inconsistent', &
           source, revision, revdate)

! Call the appropriate filter option to compute increments for ensemble
if(filter_kind == 1) then
   call obs_increment_eakf(ens, ens_size, prior_mean, prior_var, &
      obs, obs_var, obs_inc, net_a)
else if(filter_kind == 2) then
   call obs_increment_enkf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
else if(filter_kind == 3) then
   call obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
else if(filter_kind == 4) then
   call obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc)
else if(filter_kind == 5) then
   call obs_increment_ran_kf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
else if(filter_kind == 6) then
   call obs_increment_det_kf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
else 
   call error_handler(E_ERR,'obs_increment', &
              'Illegal value of filter_kind in assim_tools namelist [1-6 OK]', &
              source, revision, revdate)
endif

! Add in the extra increments if doing observation space covariance inflation
if(do_obs_inflate(inflate)) obs_inc = obs_inc + inflate_inc

! To minimize regression errors, may want to sort to minimize increments
! This makes sense for any of the non-deterministic algorithms
! By doing it here, can take care of both standard non-deterministic updates
! plus non-deterministic obs space covariance inflation. This is expensive, so
! don't use it if it's not needed.
if(sort_obs_inc) then
   new_val = ens_in + obs_inc
   ! Sorting to make increments as small as possible
   call index_sort(ens_in, ens_index, ens_size)
   call index_sort(new_val, new_index, ens_size)
   do i = 1, ens_size
      obs_inc(ens_index(i)) = new_val(new_index(i)) - ens_in(ens_index(i))
      ! The following line can provide improved results; understand why
      !!!obs_inc(ens_index(i)) = new_val(new_index(i)) - ens(ens_index(i))
   end do
end if

! Get the net change in spread if obs space inflation was used
if(do_obs_inflate(inflate)) net_a = net_a * sqrt(my_cov_inflate)

end subroutine obs_increment



subroutine obs_increment_eakf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc, a)
!========================================================================
!
! EAKF version of obs increment

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), prior_mean, prior_var, obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
real(r8), intent(out) :: a

real(r8) :: new_mean, var_ratio

! Compute the new mean
if (obs_var /= 0.0_r8) then
   var_ratio = obs_var / (prior_var + obs_var)
   new_mean  = var_ratio * (prior_mean  + prior_var*obs / obs_var)
! If obs is a delta function, it becomes new value
else
   var_ratio = 0.0_r8
   new_mean  = obs
endif

! Compute sd ratio and shift ensemble
a = sqrt(var_ratio)
obs_inc = a * (ens - prior_mean) + new_mean - ens

end subroutine obs_increment_eakf


subroutine obs_increment_ran_kf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
!========================================================================
!
! Forms a random sample of the Gaussian from the update equations.
! This is very close to what a true 'ENSEMBLE' Kalman Filter would 
! look like. Note that outliers, multimodality, etc., get tossed.

integer,   intent(in)  :: ens_size
real(r8),  intent(in)  :: prior_mean, prior_var
real(r8),  intent(in)  :: ens(ens_size), obs, obs_var
real(r8),  intent(out) :: obs_inc(ens_size)

real(r8) :: new_mean, var_ratio
real(r8) :: temp_mean, temp_var, new_ens(ens_size), new_var
integer  :: i

if (obs_var /= 0.0_r8) then
   var_ratio = obs_var / (prior_var + obs_var)
   new_var = var_ratio * prior_var
   new_mean  = var_ratio * (prior_mean  + prior_var*obs / obs_var)
else
   var_ratio = 0.0_r8
   new_var = var_ratio * prior_var
   new_mean  = obs
endif

! Now, just from a random sample from the updated distribution
! Then adjust the mean (what about adjusting the variance?)!
! Definitely need to sort with this; sort is done in main obs_increment
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq)
   first_inc_ran_call = .false.
endif

do i = 1, ens_size
   new_ens(i) = random_gaussian(inc_ran_seq, new_mean, sqrt(prior_var*var_ratio))
end do

! Adjust the mean of the new ensemble
temp_mean = sum(new_ens) / ens_size
new_ens(:) = new_ens(:) - temp_mean + new_mean

! Compute prior variance and mean from sample
temp_var  = sum((new_ens - new_mean)**2) / (ens_size - 1)
! Adjust the variance, also
new_ens = (new_ens - new_mean) * sqrt(new_var / temp_var) + new_mean

! Get the increments
obs_inc = new_ens - ens

end subroutine obs_increment_ran_kf



subroutine obs_increment_det_kf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
!========================================================================
!
! Does a deterministic ensemble layout for the updated Gaussian.
! Note that all outliers, multimodal behavior, etc. get tossed.

integer,  intent(in)   :: ens_size
real(r8), intent(in)  :: prior_mean, prior_var
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: new_mean, var_ratio, temp_var, new_ens(ens_size), new_var
integer :: i

if (obs_var /= 0.0_r8) then
   var_ratio = obs_var / (prior_var + obs_var)
   new_var = var_ratio * prior_var
   new_mean  = var_ratio * (prior_mean  + prior_var*obs / obs_var)
else
   if (prior_var /= 0.0_r8) then
      var_ratio = 0.0_r8
      new_var = var_ratio * prior_var
      new_mean  = obs
   else
      call error_handler(E_ERR,'obs_increment_det_kf', &
           'Both obs_var and prior_var are zero. This is inconsistent', &
           source, revision, revdate)
   endif
endif

! Want a symmetric distribution with kurtosis 3 and variance new_var and mean new_mean
if(ens_size /= 20) then
   write(*, *) 'EXPERIMENTAL version obs_increment_det_kf only works for ens_size 20 now'
   stop
endif

! This has kurtosis of 3.0, verify again from initial uniform
!new_ens(1) = -2.146750_r8
!new_ens(2) = -1.601447_r8
!new_ens(3) = -1.151582_r8
!new_ens(4) = -0.7898650_r8
!new_ens(5) = -0.5086292_r8
!new_ens(6) = -0.2997678_r8
!new_ens(7) = -0.1546035_r8
!new_ens(8) = -6.371084E-02_r8
!new_ens(9) = -1.658448E-02_r8
!new_ens(10) = -9.175255E-04_r8

! This has kurtosis of 3.0, verify again from initial inverse gaussian
!new_ens(1) = -2.188401_r8
!new_ens(2) = -1.502174_r8
!new_ens(3) = -1.094422_r8
!new_ens(4) = -0.8052422_r8
!new_ens(5) = -0.5840152_r8
!new_ens(6) = -0.4084518_r8
!new_ens(7) = -0.2672727_r8
!new_ens(8) = -0.1547534_r8
!new_ens(9) = -6.894587E-02_r8
!new_ens(10) = -1.243549E-02_r8

! This has kurtosis of 2.0, verify again 
new_ens(1) = -1.789296_r8
new_ens(2) = -1.523611_r8
new_ens(3) = -1.271505_r8
new_ens(4) = -1.033960_r8
new_ens(5) = -0.8121864_r8
new_ens(6) = -0.6077276_r8
new_ens(7) = -0.4226459_r8
new_ens(8) = -0.2598947_r8
new_ens(9) = -0.1242189_r8
new_ens(10) = -2.539018E-02_r8

! This has kurtosis of 1.7, verify again 
!new_ens(1) = -1.648638_r8
!new_ens(2) = -1.459415_r8
!new_ens(3) = -1.272322_r8
!new_ens(4) = -1.087619_r8
!new_ens(5) = -0.9056374_r8
!new_ens(6) = -0.7268229_r8
!new_ens(7) = -0.5518176_r8
!new_ens(8) = -0.3816142_r8
!new_ens(9) = -0.2179997_r8
!new_ens(10) = -6.538583E-02_r8
do i = 11, 20
   new_ens(i) = -1.0_r8 * new_ens(20 + 1 - i)
end do

! Right now, this ensemble has mean 0 and some variance
! Compute prior variance and mean from sample
temp_var  = sum((new_ens)**2) / (ens_size - 1)

! Adjust the variance of this ensemble to match requirements and add in the mean
new_ens = new_ens * sqrt(new_var / temp_var) + new_mean

! Get the increments
obs_inc = new_ens - ens

end subroutine obs_increment_det_kf




subroutine obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc)
!------------------------------------------------------------------------
!
! A observation space only particle filter implementation for a
! two step sequential update filter. Second version, 2 October, 2003.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: a, weight(ens_size), rel_weight(ens_size), cum_weight(0:ens_size)
real(r8) :: base, frac, new_val(ens_size), weight_sum
integer  :: i, j, indx(ens_size)

! The factor a is not defined for particle filters
a = -1.0_r8

! Begin by computing a weight for each of the prior ensemble members
do i = 1, ens_size
   weight(i) = exp(-1.0_r8 * (ens(i) - obs)**2 / (2.0_r8 * obs_var))
end do

! Compute relative weight for each ensemble member
weight_sum = sum(weight)
do i = 1, ens_size
   rel_weight(i) = weight(i) / weight_sum
end do

! Compute cumulative weights at boundaries
cum_weight(0) = 0.0_r8
do i = 1, ens_size
   cum_weight(i) = cum_weight(i - 1) + rel_weight(i)
!   write(*,'(1x,i3,3(e10.4,1x))') i, weight(i), rel_weight(i), cum_weight(i)
end do
! Fix up for round-off error if any
cum_weight(ens_size) = 1.0_r8

! Do a deterministic implementation: just divide interval into ens_size parts and see
! which interval this is in (careful to offset; not start at 0)
base = 1.0_r8 / (ens_size * 2.0_r8)

do i = 1, ens_size

   frac = base + (i - 1.0_r8) / ens_size

   ! Now search in the cumulative range to see where this frac falls
   ! Can make this search more efficient by limiting base
   do j = 1, ens_size
      if(cum_weight(j - 1) < frac .and. frac < cum_weight(j)) then
         indx(i) = j
!         write(*, *) i, frac, 'gets index ', j
         goto 111
      end if
   end do

111 continue

end do

! Set the new values for the ensemble members
do i = 1, ens_size
   new_val(i) = ens(indx(i))
!   write(*, *) 'new_val ', i, new_val(i)
end do

! Generate increments
obs_inc = new_val - ens

end subroutine obs_increment_particle



subroutine obs_increment_enkf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment_enkf(ens, ens_size, obs, obs_var, obs_inc)
!

! ENKF version of obs increment

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), prior_mean, prior_var, obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: a, obs_var_inv, prior_var_inv, new_var, new_mean(ens_size)
! real(r8) :: sx, s_x2
real(r8) :: temp_mean, temp_obs(ens_size)
integer  :: i

! The factor a is not defined for kernel filters
a = -1.0_r8

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0_r8 / obs_var

prior_var_inv = 1.0_r8 / prior_var
new_var       = 1.0_r8 / (prior_var_inv + obs_var_inv)

! If this is first time through, need to initialize the random sequence
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq)
   first_inc_ran_call = .false.
endif

! Generate perturbed obs
do i = 1, ens_size
    temp_obs(i) = random_gaussian(inc_ran_seq, obs, sqrt(obs_var))
end do

! Move this so that it has original obs mean
temp_mean = sum(temp_obs) / ens_size
temp_obs(:) = temp_obs(:) - temp_mean + obs

! Loop through pairs of priors and obs and compute new mean
do i = 1, ens_size
   new_mean(i) = new_var * (prior_var_inv * ens(i) + temp_obs(i) / obs_var)
   obs_inc(i)  = new_mean(i) - ens(i)
end do

! Can also adjust mean (and) variance of final sample; works fine
!sx         = sum(new_mean)
!s_x2       = sum(new_mean * new_mean)
!temp_mean = sx / ens_size
!temp_var  = (s_x2 - sx**2 / ens_size) / (ens_size - 1)
!new_mean = (new_mean - temp_mean) * sqrt(new_var / temp_var) + updated_mean
!obs_inc = new_mean - ens


end subroutine obs_increment_enkf



subroutine obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
!

! Kernel version of obs increment

integer, intent(in)             :: ens_size
real(r8), intent(in)            :: ens(ens_size), obs, obs_var
real(r8), intent(out)           :: obs_inc(ens_size)

real(r8) :: obs_var_inv
real(r8) :: prior_mean, prior_cov_inv, new_cov, prior_cov
real(r8) :: sx, s_x2
real(r8) :: weight(ens_size), new_mean(ens_size)
real(r8) :: cum_weight, total_weight, cum_frac(ens_size)
real(r8) :: unif, norm, new_member(ens_size)

integer :: i, j, kernel

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0_r8 / obs_var

! Compute prior mean and covariance
sx         = sum(ens)
s_x2       = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov  = sum((ens - prior_mean)**2) / (ens_size - 1)

prior_cov     = prior_cov / 10.0_r8     ! For kernels, scale the prior covariance
prior_cov_inv = 1.0_r8 / prior_cov

! Compute new covariance once for these kernels
new_cov = 1.0_r8 / (prior_cov_inv + obs_var_inv)

! New mean is computed ens_size times as is weight
do i = 1, ens_size
   new_mean(i) = new_cov*(prior_cov_inv * ens(i) + obs / obs_var)
   weight(i) =  2.71828_r8 ** (-0.5_r8 * (ens(i)**2 * prior_cov_inv + &
      obs**2 * obs_var_inv - new_mean(i)**2 / new_cov))
end do

! Compute total weight
total_weight = sum(weight)
cum_weight   = 0.0_r8
do i = 1, ens_size
   cum_weight  = cum_weight + weight(i)
   cum_frac(i) = cum_weight / total_weight
end do

! If this is first time through, need to initialize the random sequence
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq)
   first_inc_ran_call = .false.
endif

! Generate a uniform random number and a Gaussian for each new member
do i = 1, ens_size
   unif = random_uniform(inc_ran_seq)
   ! Figure out which kernel it's in
   do j = 1, ens_size
      if(unif < cum_frac(j)) then
         kernel = j
         goto 10
      end if
   end do
10 continue

   ! Next calculate a unit normal in this kernel
   norm = random_gaussian(inc_ran_seq, 0.0_r8, sqrt(new_cov))
   ! Now generate the new ensemble member
   new_member(i) = new_mean(kernel) + norm
end do

! Generate the increments
obs_inc = new_member - ens

end subroutine obs_increment_kernel



subroutine update_from_obs_inc(obs, obs_prior_mean, obs_prior_var, obs_inc, &
               state, ens_size, state_inc, reg_coef, net_a, correl_out)
!========================================================================

! Does linear regression of a state variable onto an observation and
! computes state variable increments from observation increments

integer,            intent(in)    :: ens_size
real(r8),           intent(in)    :: obs(ens_size), obs_inc(ens_size)
real(r8),           intent(in)    :: obs_prior_mean, obs_prior_var
real(r8),           intent(in)    :: state(ens_size)
real(r8),           intent(out)   :: state_inc(ens_size), reg_coef
real(r8),           intent(inout) :: net_a
real(r8), optional, intent(inout) :: correl_out

real(r8) :: t(ens_size), obs_state_cov
real(r8) :: restoration_inc(ens_size), state_mean, state_var, correl
real(r8) :: factor, exp_true_correl, mean_factor

! For efficiency, just compute regression coefficient here unless correl is needed
t = obs - obs_prior_mean
obs_state_cov = sum(t * state)

if (obs_prior_var /= 0.0_r8) then
   reg_coef = obs_state_cov/obs_prior_var
else
   reg_coef = 0.0_r8
endif

! If correl_out is present, need correl for adaptive inflation
! Also needed for file correction below
if(present(correl_out) .or. sampling_error_correction) then
   state_var = sum(state * state) - sum(state)**2 / ens_size
   if(obs_prior_var * state_var <= 0.0_r8) then
      correl = 0.0_r8
   else
      correl = obs_state_cov / sqrt(obs_prior_var * state_var)
   endif
   if(correl >  1.0_r8) correl =  1.0_r8
   if(correl < -1.0_r8) correl = -1.0_r8
endif
if(present(correl_out)) correl_out = correl


! BEGIN TEST OF CORRECTION FROM FILE +++++++++++++++++++++++++++++++++++++++++++++++++
! Get the expected actual correlation and the regression weight reduction factor
if(sampling_error_correction) then
   call get_correction_from_file(ens_size, correl, mean_factor, exp_true_correl)
   !write(*, *) correl, exp_true_correl, mean_factor
   ! Watch out for division by zero; if correl is really small regression is safely 0
   if(abs(correl) > 0.001) then
      reg_coef = reg_coef * (exp_true_correl / correl) * mean_factor
   else
      reg_coef = 0.0_r8
   endif
   correl = exp_true_correl
endif

! END TEST OF CORRECTION FROM FILE +++++++++++++++++++++++++++++++++++++++++++++++++



! Then compute the increment as product of reg_coef and observation space increment
state_inc = reg_coef * obs_inc


! Spread restoration algorithm option
if(spread_restoration) then
   ! Don't use this to reduce spread at present (should revisit this line)
   if(net_a > 1.0_r8) net_a = 1.0_r8

   ! Default restoration increment is 0.0
   restoration_inc = 0.0_r8

   ! Compute the factor by which to inflate
   ! These come from correl_error.f90 in system_simulation and the files ens??_pairs and
   ! ens_pairs_0.5 in work under system_simulation. Assume a linear reduction from 1
   ! as a function of the net_a. Assume that the slope of this reduction is a function of
   ! the reciprocal of the ensemble_size (slope = 0.80 / ens_size). These are empirical
   ! for now.
   factor = 1.0_r8 / (1.0_r8 - (1.0_r8 - net_a) * (0.8_r8 / ens_size)) - 1.0_r8

   ! Variance restoration
   state_mean = sum(state) / ens_size
   restoration_inc = factor * (state - state_mean)
   state_inc = state_inc + restoration_inc
endif

end subroutine update_from_obs_inc


!------------------------------------------------------------------------

subroutine get_correction_from_file(ens_size, scorrel, mean_factor, expected_true_correl)

integer,   intent(in) :: ens_size
real(r8),  intent(in) :: scorrel
real(r8), intent(out) :: mean_factor, expected_true_correl

! Reads in a regression error file for a give ensemble_size and uses interpolation
! to get correction factor into the file

integer  :: iunit, i, low_indx, high_indx
real(r8) :: temp, temp2, correl, fract, low_correl, low_exp_correl, low_alpha
real(r8) :: high_correl, high_exp_correl, high_alpha
character(len = 20) :: correction_file_name

if(first_get_correction) then
   ! Compute the file name for this ensemble size
   if(ens_size < 10) then
      write(correction_file_name, 11) 'final_full.', ens_size
   else if(ens_size < 100) then
      write(correction_file_name, 21) 'final_full.', ens_size
   else if(ens_size < 1000) then
      write(correction_file_name, 31) 'final_full.', ens_size
   else if(ens_size < 10000) then
      write(correction_file_name, 41) 'final_full.', ens_size
   else
      write(errstring,*)'Trying to use ',ens_size,' model states -- too many.'
      call error_handler(E_MSG,'get_correction_from_file',errstring,source,revision,revdate)
      call error_handler(E_ERR,'get_correction_from_file','Use less than 10000 ensemble.',source,revision,revdate)

    11   format(a11, i1)
    21   format(a11, i2)
    31   format(a11, i3)
    41   format(a11, i4)
   endif
 
   ! Make sure that the correction file exists, else an error
   if(.not. file_exist(correction_file_name)) then
      write(errstring,*) 'Correction file ', correction_file_name, ' does not exist'
      call error_handler(E_ERR,'get_correction_from_file',errstring,source,revision,revdate)
   endif

   ! Read in file to get the expected value of the true correlation given the sample
   iunit = get_unit()
   open(unit = iunit, file = correction_file_name)
   do i = 1, 200
      read(iunit, *) temp, temp2, exp_true_correl(i), alpha(i)
   end do
   close(iunit)

   first_get_correction = .false.
endif


! First quick modification of correl to expected true correl for test (should interp)
if(scorrel < -1.0_r8) then
   correl = -1.0_r8
   mean_factor = 1.0_r8
else if(scorrel > 1.0_r8) then
   correl = 1.0_r8
   mean_factor = 1.0_r8
else if(scorrel <= -0.995_r8) then
   fract = (scorrel + 1.0_r8) / 0.05_r8
   correl = (exp_true_correl(1) + 1.0_r8) * fract - 1.0_r8
   mean_factor = (alpha(1) - 1.0_r8) * fract + 1.0_r8
else if(scorrel >= 0.995_r8) then
   fract = (scorrel - 0.995_r8) / 0.05_r8
   correl = (1.0_r8 - exp_true_correl(200)) * fract + exp_true_correl(200)
   mean_factor = (1.0_r8 - alpha(200)) * fract + alpha(200)
else
   low_indx = (scorrel + 0.995_r8) / 0.01_r8 + 1.0_r8
   low_correl = -0.995_r8 + (low_indx - 1) * 0.01
   low_exp_correl = exp_true_correl(low_indx)
   low_alpha = alpha(low_indx)
   high_indx = low_indx + 1
   high_correl = low_correl + 0.01
   high_exp_correl = exp_true_correl(high_indx)
   high_alpha = alpha(low_indx)
   fract = (scorrel - low_correl) / (high_correl - low_correl)
   correl = (high_exp_correl - low_exp_correl) * fract + low_exp_correl
   mean_factor = (high_alpha - low_alpha) * fract + low_alpha
endif

expected_true_correl = correl

end subroutine get_correction_from_file

!------------------------------------------------------------------------

! Temporary testing stuff

subroutine my_broadcast_send1(owner, obs_prior, obs_inc, orig_obs_prior_mean, &
         orig_obs_prior_var, qc)

integer, intent(in) :: owner
real(r8), intent(inout) :: obs_prior(:), obs_inc(:), orig_obs_prior_mean(:)
real(r8), intent(inout) :: orig_obs_prior_var(:), qc

real(r8) :: temp(size(obs_inc) + size(orig_obs_prior_mean) + size(orig_obs_prior_var) + 1)

! Pack stuff together into two arrays and use broadcast_send
temp(1:size(obs_inc)) = obs_inc
temp(size(obs_inc) + 1: size(obs_inc) + size(orig_obs_prior_mean)) = orig_obs_prior_mean
temp(size(obs_inc) + size(orig_obs_prior_mean) + 1:size(temp) - 1) = orig_obs_prior_var
temp(size(temp)) = qc
call broadcast_send(owner, obs_prior, temp)

end subroutine my_broadcast_send1

!------------------------------------------------------------------------


subroutine my_broadcast_recv1(owner, obs_prior, obs_inc, orig_obs_prior_mean, &
         orig_obs_prior_var, qc)

integer, intent(in) :: owner
real(r8), intent(inout) :: obs_prior(:), obs_inc(:), orig_obs_prior_mean(:)
real(r8), intent(inout) :: orig_obs_prior_var(:), qc

real(r8) :: temp(size(obs_inc) + size(orig_obs_prior_mean) + size(orig_obs_prior_var) + 1)

! Pack stuff together into two arrays and use broadcast_recv
call broadcast_recv(owner, obs_prior, temp)
obs_inc = temp(1:size(obs_inc))
orig_obs_prior_mean = temp(size(obs_inc) + 1: size(obs_inc) + size(orig_obs_prior_mean))
orig_obs_prior_var = temp(size(obs_inc) + size(orig_obs_prior_mean) + 1: size(temp) - 1)
qc = temp(size(temp))

end subroutine my_broadcast_recv1

!------------------------------------------------------------------------

subroutine my_broadcast_send2(owner, obs_prior, obs_inc, ss_inflate, ss_inflate_sd, qc)

integer, intent(in) :: owner
real(r8), intent(inout) :: obs_prior(:), obs_inc(:), ss_inflate, ss_inflate_sd, qc

real(r8) :: temp(size(obs_inc) + 3)

! Pack stuff together into two arrays and use broadcast_send
temp(1:size(obs_inc)) = obs_inc
temp(size(obs_inc) + 1) = ss_inflate
temp(size(obs_inc) + 2) = ss_inflate_sd
temp(size(temp)) = qc
call broadcast_send(owner, obs_prior, temp)

end subroutine my_broadcast_send2

!------------------------------------------------------------------------

subroutine my_broadcast_recv2(owner, obs_prior, obs_inc, ss_inflate, ss_inflate_sd, qc)

integer, intent(in) :: owner
real(r8), intent(inout) :: obs_prior(:), obs_inc(:), ss_inflate, ss_inflate_sd, qc

real(r8) :: temp(size(obs_inc) + 3)

call broadcast_recv(owner, obs_prior, temp)

! Pack stuff together into two arrays and use broadcast_recv
obs_inc = temp(1:size(obs_inc))
ss_inflate = temp(size(obs_inc) + 1)
ss_inflate_sd = temp(size(obs_inc) + 2)
qc = temp(size(temp))

end subroutine my_broadcast_recv2

!------------------------------------------------------------------------

subroutine my_broadcast_send3(owner, obs_prior, obs_inc, obs_qc)

integer, intent(in) :: owner
real(r8), intent(inout) :: obs_prior(:), obs_inc(:), obs_qc

real(r8) :: temp(size(obs_inc) + 1)

! Pack the increment and qc
temp(1:size(obs_inc)) = obs_inc
temp(size(temp)) = obs_qc

call broadcast_send(owner, obs_prior, temp)

end subroutine my_broadcast_send3

!------------------------------------------------------------------------

subroutine my_broadcast_recv3(owner, obs_prior, obs_inc, obs_qc)

integer, intent(in) :: owner
real(r8), intent(inout) :: obs_prior(:), obs_inc(:), obs_qc

real(r8) :: temp(size(obs_inc) + 1)

call broadcast_recv(owner, obs_prior, temp)

! Extract the increment and qc
obs_inc = temp(1:size(obs_inc))
obs_qc = temp(size(temp))

end subroutine my_broadcast_recv3

!------------------------------------------------------------------------


!========================================================================
! end module assim_tools_mod
!========================================================================

end module assim_tools_mod
