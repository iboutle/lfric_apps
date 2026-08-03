! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!-------------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!-------------------------------------------------------------------------------

! Cloud Aerosol Interacting Microphysics (CASIM)
! Mechanistic activation of cloud droplets on the aerosol supplied by the host
! model.
!
! These routines are single column versions of the routines held in the UM file
! casim_activation_in_um_mod.F90. The UM works on the (i,j,k) fields of a whole
! domain, whereas an LFRic kernel only ever sees one column, so the horizontal
! dimensions have been dropped and the arrays are simply indexed by level.
!
! The UM holds four near identical copies of the aerosol examination and of the
! activation loop, one for each of the aerosol sources it supports (prognostic
! CASIM tracers, murk, one way UKCA and two way UKCA). The differences between
! the copies are entirely in how the aerosol mass and number are read out of the
! host model. LFRic has already unified that step in
! aerosol_extract_convert_mod, which returns the mass and number of each mode
! whatever the source, so the four copies collapse into a single set of routines
! here:
!
!   examine_aerosol_column     <- the aerosol part of the UM routines
!                                 examine_ukca_aerosol_column,
!                                 examine_aerosol_tracer_column and
!                                 examine_murk_aerosol_column
!   examine_processing_column  <- the l_process part of the UM routines
!                                 examine_aerosol_tracer_column and
!                                 examine_2way_ukca_aerosol_column
!   activate_column            <- the UM routines activate_column_ukca,
!                                 activate_column_murk,
!                                 activate_column_tracers and
!                                 activate_column_2way_ukca
!
! The routines fill the module level derived types held by CASIM itself
! (aerophys, aerochem, aeroact, dustphys, dustchem, dustliq), which the CASIM
! routine activate then reads. Because that state is module level, an invoke
! containing this kernel must not be threaded over cells.

module casim_activation_in_um_mod

use micro_main, only: aerophys, aeroact, aerochem, dustact, dustphys,          &
                      dustchem, dustliq
use um_types,   only: real_umphys

implicit none

character(len=*), parameter, private ::                                        &
  ModuleName='CASIM_ACTIVATION_IN_UM_MOD'

private
public :: examine_aerosol_column, examine_processing_column, activate_column,  &
          sulphate_mode_bk

! Hygroscopicity parameter of ammonium sulphate. The UM routine
! examine_murk_aerosol_column assumes the murk aerosol is ammonium sulphate and
! uses this value; there is no equivalent value available from the murk fields
! themselves, so the caller supplies it for the murk aerosol source.
real(kind=real_umphys), parameter :: sulphate_mode_bk = 0.4_real_umphys

! Values used to initialise the activated aerosol properties
real(kind=real_umphys), parameter, private :: zero       = 0.0_real_umphys
real(kind=real_umphys), parameter, private :: rcrit_high = 999.0_real_umphys

contains

!-------------------------------------------------------------------------------
! Set up the interstitial aerosol seen by CASIM
!-------------------------------------------------------------------------------
! Copies the mass and number of each aerosol mode into the CASIM aerophys and
! aerochem structures, and resets the activated aerosol structures aeroact and
! dustact ready for the activation calculation.
!
! The UM reads the mass and number for each mode straight out of the host model
! tracer arrays; here they have already been extracted and converted by
! aerosol_extract_convert_mod, so they are simply passed in.

subroutine examine_aerosol_column( nlayers,                                    &
                                   aitken_sol_mass,  aitken_sol_number,        &
                                   accum_sol_mass,   accum_sol_number,         &
                                   coarse_sol_mass,  coarse_sol_number,        &
                                   accum_dust_mass,  accum_dust_number,        &
                                   coarse_dust_mass, coarse_dust_number,       &
                                   aitken_sol_bk, accum_sol_bk, coarse_sol_bk, &
                                   l_set_bk )

use mphys_switches,  only: aero_index, l_warm
use lognormal_funcs, only: MNtoRm
use thresholds,      only: ccn_tidy, aeromass_small, aeronumber_small

implicit none

!-------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------

! Number of model levels
integer, intent(in) :: nlayers

! Aitken mode soluble mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(in) :: aitken_sol_mass(nlayers)
real(kind=real_umphys), intent(in) :: aitken_sol_number(nlayers)

! Accumulation mode soluble mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(in) :: accum_sol_mass(nlayers)
real(kind=real_umphys), intent(in) :: accum_sol_number(nlayers)

! Coarse mode soluble mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(in) :: coarse_sol_mass(nlayers)
real(kind=real_umphys), intent(in) :: coarse_sol_number(nlayers)

! Accumulation mode dust mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(in) :: accum_dust_mass(nlayers)
real(kind=real_umphys), intent(in) :: accum_dust_number(nlayers)

! Coarse mode dust mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(in) :: coarse_dust_mass(nlayers)
real(kind=real_umphys), intent(in) :: coarse_dust_number(nlayers)

! Activation parameter Bk for each soluble mode []
real(kind=real_umphys), intent(in) :: aitken_sol_bk(nlayers)
real(kind=real_umphys), intent(in) :: accum_sol_bk(nlayers)
real(kind=real_umphys), intent(in) :: coarse_sol_bk(nlayers)

! Whether the supplied Bk values should be used. The UM only derives Bk from
! the host model aerosol for the UKCA and murk sources; for the prognostic
! CASIM aerosol the values already held by CASIM are left alone.
logical, intent(in) :: l_set_bk

!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------

! Aerosol density [kg m-3]
real(kind=real_umphys) :: density

! Mass, number and hygroscopicity of the mode currently being examined
real(kind=real_umphys) :: mode_m, mode_n, mode_bk

integer :: k     ! Loop counter over levels
integer :: imode ! Loop counter over CASIM aerosol modes

real(kind=real_umphys), parameter :: eps_1 = epsilon(1.0_real_umphys)

!-------------------------------------------------------------------------------
! End of declarations and start of subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------
! Soluble aerosol species
!-------------------------------------------------
do k = 1, nlayers

  density = aerochem(k)%density(1)

  do imode = 1, aero_index%nccn

    mode_n  = zero
    mode_m  = zero
    mode_bk = zero

    if ( imode == aero_index%i_aitken ) then
      mode_n  = aitken_sol_number(k)
      mode_m  = aitken_sol_mass(k)
      mode_bk = aitken_sol_bk(k)
    else if ( imode == aero_index%i_accum ) then
      mode_n  = accum_sol_number(k)
      mode_m  = accum_sol_mass(k)
      mode_bk = accum_sol_bk(k)
    else if ( imode == aero_index%i_coarse ) then
      mode_n  = coarse_sol_number(k)
      mode_m  = coarse_sol_mass(k)
      mode_bk = coarse_sol_bk(k)
    end if

    ! aerosol_extract_convert returns missing data for the hygroscopicity of a
    ! mode which holds no aerosol. That mode will fail the test below, but
    ! guard against a missing value ever reaching CASIM.
    if ( mode_bk < zero ) mode_bk = zero

    if ( mode_n > ccn_tidy .and. mode_m > ccn_tidy*eps_1 ) then

      aerophys(k)%n(imode)  = mode_n
      aerophys(k)%m(imode)  = mode_m
      aerophys(k)%rd(imode) = MNtoRm( mode_m, mode_n, density,                 &
                                      aerophys(k)%sigma(imode) )
      if ( l_set_bk ) aerochem(k)%bk(imode) = mode_bk

    else

      ! Zero the processes.
      aerophys(k)%n(imode)  = zero
      aerophys(k)%m(imode)  = zero
      aerophys(k)%rd(imode) = zero
      if ( l_set_bk ) aerochem(k)%bk(imode) = zero

    end if ! mode_n > ccn_tidy etc

  end do ! imode

  ! Set other aerosol properties for this level
  aeroact(k)%nact       = zero
  aeroact(k)%mact       = zero
  aeroact(k)%rcrit      = rcrit_high
  aeroact(k)%mact_mean  = zero
  aeroact(k)%nact1      = zero
  aeroact(k)%mact1      = zero
  aeroact(k)%rcrit1     = rcrit_high
  aeroact(k)%mact1_mean = zero
  aeroact(k)%nact2      = zero
  aeroact(k)%mact2      = zero
  aeroact(k)%rcrit2     = rcrit_high
  aeroact(k)%mact2_mean = zero
  aeroact(k)%nratio1    = zero
  aeroact(k)%nratio2    = zero

end do ! k

!-------------------------------------------------
! Insoluble aerosol species
!-------------------------------------------------
if ( .not. l_warm ) then

  do k = 1, nlayers

    density = dustchem(k)%density(1)

    ! Examine interstitial dust
    do imode = 1, aero_index%nin

      mode_n = zero
      mode_m = zero

      if ( imode == aero_index%i_coarse_dust ) then
        mode_n = coarse_dust_number(k)
        mode_m = coarse_dust_mass(k)
      end if
      if ( imode == aero_index%i_accum_dust ) then
        mode_n = accum_dust_number(k)
        mode_m = accum_dust_mass(k)
      end if

      if ( mode_m < zero .or. mode_n < zero ) then
        if ( mode_m < zero ) mode_m = aeromass_small
        if ( mode_n < zero ) mode_n = aeronumber_small
      end if

      dustphys(k)%n(imode)  = mode_n
      dustphys(k)%m(imode)  = mode_m
      dustphys(k)%rd(imode) = MNtoRm( mode_m, mode_n, density,                 &
                                      dustphys(k)%sigma(imode) )

    end do ! imode

    ! Examine activated dust - initialise to zero/defaults
    dustact(k)%nact       = zero
    dustact(k)%mact       = zero
    dustact(k)%rcrit      = rcrit_high
    dustact(k)%mact_mean  = zero
    dustact(k)%nact1      = zero
    dustact(k)%nratio1    = zero
    dustact(k)%mact1      = zero
    dustact(k)%rcrit1     = rcrit_high
    dustact(k)%mact1_mean = zero
    dustact(k)%nact2      = zero
    dustact(k)%nratio2    = zero
    dustact(k)%mact2      = zero
    dustact(k)%rcrit2     = rcrit_high
    dustact(k)%mact2_mean = zero
    dustact(k)%nact3      = zero
    dustact(k)%nratio3    = zero
    dustact(k)%mact3      = zero
    dustact(k)%rcrit3     = rcrit_high
    dustact(k)%mact3_mean = zero

  end do ! k

end if ! .not. l_warm

end subroutine examine_aerosol_column

!-------------------------------------------------------------------------------
! Set up the already activated aerosol seen by CASIM
!-------------------------------------------------------------------------------
! For aerosol processing runs the host model carries the mass and number of
! aerosol which is already inside the hydrometeors. This routine partitions that
! aerosol between the cloud droplets and the rain drops and copies it into the
! CASIM aeroact and dustliq structures.
!
! This must be called after examine_aerosol_column, which resets aeroact.

subroutine examine_processing_column( nlayers,                                 &
                                      cloud_number,  cloud_mass,               &
                                      rain_number,   rain_mass,                &
                                      active_sol_liquid, active_sol_rain,      &
                                      active_sol_ice,    active_sol_number,    &
                                      active_insol_liquid, active_insol_ice,   &
                                      active_insol_number )

use mphys_switches,   only: l_separate_rain, l_warm
use mphys_parameters, only: sigma_arc
use lognormal_funcs,  only: MNtoRm
use thresholds,       only: nl_tidy, nr_tidy, qr_tidy
use casim_set_dependent_switches_mod, only: l_passivenumbers,                  &
                                            l_passivenumbers_ice

implicit none

!-------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------

! Number of model levels
integer, intent(in) :: nlayers

! Cloud droplet number [kg-1] and mass [kg kg-1]
real(kind=real_umphys), intent(in) :: cloud_number(nlayers)
real(kind=real_umphys), intent(in) :: cloud_mass(nlayers)

! Rain drop number [kg-1] and mass [kg kg-1]
real(kind=real_umphys), intent(in) :: rain_number(nlayers)
real(kind=real_umphys), intent(in) :: rain_mass(nlayers)

! Soluble aerosol mass activated in cloud, rain and ice [kg kg-1]
real(kind=real_umphys), intent(in) :: active_sol_liquid(nlayers)
real(kind=real_umphys), intent(in) :: active_sol_rain(nlayers)
real(kind=real_umphys), intent(in) :: active_sol_ice(nlayers)

! Soluble aerosol number activated in the hydrometeors [kg-1]
real(kind=real_umphys), intent(in) :: active_sol_number(nlayers)

! Insoluble aerosol mass activated in liquid and ice [kg kg-1]
real(kind=real_umphys), intent(in) :: active_insol_liquid(nlayers)
real(kind=real_umphys), intent(in) :: active_insol_ice(nlayers)

! Insoluble aerosol number activated in the hydrometeors [kg-1]
real(kind=real_umphys), intent(in) :: active_insol_number(nlayers)

!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------

! Aerosol density [kg m-3]
real(kind=real_umphys) :: density

! Activated aerosol masses in cloud, rain and ice [kg kg-1]
real(kind=real_umphys) :: mac, mar, mact, maai
real(kind=real_umphys) :: mad, madl

! Activated aerosol number and the hydrometeor number it is shared over [kg-1]
real(kind=real_umphys) :: ntot, nhtot

! Ratios used to partition the aerosol
real(kind=real_umphys) :: ratio_ali, ratio_dli
real(kind=real_umphys) :: nratio_l, nratio_r

! Mean radius of the activated distribution [m]
real(kind=real_umphys) :: rm_arc

! Working values for the rain category
real(kind=real_umphys) :: rcrit2, mact2

logical :: l_condition   ! There is activated aerosol in the cloud droplets
logical :: l_condition_r ! There is activated aerosol in the rain drops

integer :: k ! Loop counter over levels

real(kind=real_umphys), parameter :: eps_1 = epsilon(1.0_real_umphys)

!-------------------------------------------------------------------------------
! End of declarations and start of subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------
! Activated soluble aerosol
!-------------------------------------------------
do k = 1, nlayers

  density = aerochem(k)%density(1)

  mac = active_sol_liquid(k)

  if ( l_separate_rain ) then
    mar  = active_sol_rain(k)
    mact = mac + mar
    l_condition   = mac > eps_1 .and. cloud_number(k) > nl_tidy
    l_condition_r = mar > eps_1 .and. rain_number(k) > nr_tidy .and.           &
                    rain_mass(k) > qr_tidy
  else
    mar  = zero
    mact = mac
    l_condition   = mac > eps_1 .and.                                          &
                    cloud_number(k) + rain_number(k) > nl_tidy
    l_condition_r = mac > eps_1 .and. rain_number(k) > nr_tidy .and.           &
                    rain_mass(k) > qr_tidy
  end if

  if ( l_warm ) then
    maai = zero
  else
    maai = active_sol_ice(k)
  end if

  ratio_ali = mact / ( mact + maai + eps_1 )

  nhtot = cloud_number(k) + rain_number(k)

  if ( l_passivenumbers ) then
    ntot = active_sol_number(k) * ratio_ali
  else
    ntot = nhtot
  end if

  l_condition = l_condition .and. ntot > eps_1

  ! prevent possible later division by zero
  nhtot    = nhtot + eps_1
  nratio_l = cloud_number(k) / nhtot
  nratio_r = rain_number(k) / nhtot

  if ( l_condition ) then
    aeroact(k)%nact      = ntot
    aeroact(k)%mact      = mac
    aeroact(k)%rcrit     = zero
    aeroact(k)%mact_mean = aeroact(k)%mact / ( aeroact(k)%nact + eps_1 )

    ! Get mean radius of distribution
    rm_arc = MNtoRm( aeroact(k)%mact, aeroact(k)%nact, density, sigma_arc )

    aeroact(k)%rd    = rm_arc
    aeroact(k)%sigma = sigma_arc
  end if

  if ( l_condition_r ) then
    if ( l_separate_rain ) then
      ! Separate rain category
      rcrit2 = zero
      mact2  = mar
    else if ( .not. l_condition ) then
      ! No cloud here
      rcrit2 = zero
      mact2  = mac
    else
      ! Diagnostic partitioning of aerosol between cloud and rain.
      rcrit2 = zero
      mact2  = aeroact(k)%mact * rain_mass(k) / ( cloud_mass(k) + rain_mass(k) )
    end if

    aeroact(k)%nact2      = ntot * nratio_r
    aeroact(k)%rcrit2     = rcrit2
    aeroact(k)%mact2      = mact2
    aeroact(k)%mact2_mean = aeroact(k)%mact2 / ( aeroact(k)%nact2 + eps_1 )
  end if

  aeroact(k)%nact1      = max( zero, aeroact(k)%nact - aeroact(k)%nact2 )
  aeroact(k)%mact1      = max( zero, aeroact(k)%mact - aeroact(k)%mact2 )
  aeroact(k)%rcrit1     = zero
  aeroact(k)%mact1_mean = aeroact(k)%mact1 / ( aeroact(k)%nact1 + eps_1 )

  if ( cloud_number(k) > eps_1 )                                               &
      aeroact(k)%nratio1 = aeroact(k)%nact1 / cloud_number(k)
  if ( rain_number(k) > eps_1 )                                                &
      aeroact(k)%nratio2 = aeroact(k)%nact2 / rain_number(k)

end do ! k

!-------------------------------------------------
! Activated insoluble aerosol in the liquid
!-------------------------------------------------
if ( .not. l_warm ) then

  do k = 1, nlayers

    madl = active_insol_liquid(k)
    mad  = active_insol_ice(k)

    ratio_dli = madl / ( mad + madl + eps_1 )

    nhtot    = cloud_number(k) + rain_number(k) + eps_1
    nratio_l = cloud_number(k) / nhtot
    nratio_r = rain_number(k) / nhtot

    if ( l_passivenumbers_ice ) then
      ntot = active_insol_number(k) * ratio_dli
    else
      ntot = cloud_number(k) + rain_number(k)
    end if

    ! Initialise to zero/defaults
    dustliq(k)%nact       = zero
    dustliq(k)%mact       = zero
    dustliq(k)%rcrit      = rcrit_high
    dustliq(k)%mact_mean  = zero
    dustliq(k)%nact1      = zero
    dustliq(k)%nratio1    = zero
    dustliq(k)%mact1      = zero
    dustliq(k)%rcrit1     = rcrit_high
    dustliq(k)%mact1_mean = zero
    dustliq(k)%nact2      = zero
    dustliq(k)%nratio2    = zero
    dustliq(k)%mact2      = zero
    dustliq(k)%rcrit2     = rcrit_high
    dustliq(k)%mact2_mean = zero
    dustliq(k)%nact3      = zero
    dustliq(k)%nratio3    = zero
    dustliq(k)%mact3      = zero
    dustliq(k)%rcrit3     = rcrit_high
    dustliq(k)%mact3_mean = zero

    if ( madl > eps_1 .and. ntot > nr_tidy ) then

      ! Equal partitioning of aerosol distribution across all liquid species
      dustliq(k)%nact      = ntot
      dustliq(k)%mact      = madl
      dustliq(k)%rcrit     = zero
      dustliq(k)%mact_mean = dustliq(k)%mact / ( dustliq(k)%nact + eps_1 )

      if ( nratio_r > eps_1 ) then
        dustliq(k)%nact2      = ntot * nratio_r
        dustliq(k)%rcrit2     = zero
        dustliq(k)%mact2      = madl * nratio_r
        dustliq(k)%mact2_mean = dustliq(k)%mact2 / ( dustliq(k)%nact2 + eps_1 )
        dustliq(k)%nratio2    = max( zero, min( 1.0_real_umphys,               &
            ntot / ( cloud_number(k) + rain_number(k) ) * nratio_r ) )
      end if

      if ( nratio_l > eps_1 ) then
        dustliq(k)%nact1      = max( zero, dustliq(k)%nact - dustliq(k)%nact2 )
        dustliq(k)%mact1      = max( zero, dustliq(k)%mact - dustliq(k)%mact2 )
        dustliq(k)%rcrit1     = zero
        dustliq(k)%mact1_mean = dustliq(k)%mact1 / ( dustliq(k)%nact1 + eps_1 )
        dustliq(k)%nratio1    = max( zero, min( 1.0_real_umphys,               &
            1.0_real_umphys - dustliq(k)%nratio2 ) )
      end if

    end if ! madl > eps_1 .and. ntot > nr_tidy

  end do ! k

end if ! .not. l_warm

end subroutine examine_processing_column

!-------------------------------------------------------------------------------
! Activate cloud droplets on the aerosol
!-------------------------------------------------------------------------------
! Works down the column calling the CASIM routine activate wherever the liquid
! cloud has been increased by physical processes outside of CASIM, and removing
! droplets again wherever the liquid cloud has been reduced.
!
! examine_aerosol_column, and where processing is active
! examine_processing_column, must have been called first.

subroutine activate_column( nlayers, l_process_aerosol, l_cap_drop_number,     &
                            cloud_mass_post_qtbal, cloud_mass_pre_ap2,         &
                            rho_col, t_col, p_col,                             &
                            cf_liquid, cf_liquid_pre_ap2, w_tke,               &
                            cloudnumber_col, rainnumber_col,                   &
                            aitken_sol_mass,  aitken_sol_number,               &
                            accum_sol_mass,   accum_sol_number,                &
                            coarse_sol_mass,  coarse_sol_number,               &
                            coarse_dust_mass, coarse_dust_number,              &
                            active_sol_liquid,   active_sol_number,            &
                            active_insol_liquid, active_insol_ice,             &
                            active_insol_number )

use timestep_mod,      only: timestep, recip_timestep
use mphys_switches,    only: aero_index, l_bypass_which_mode, iopt_which_mode, &
                             l_warm
use mphys_constants_mod, only: max_drop_casim
use thresholds,        only: ql_tidy, cfliq_small, aeromass_small,             &
                             aeronumber_small
use condensation,      only: dnccn_all, dmac_all, dnccnd_all, dmad_all
use activation,        only: activate
use which_mode_to_use, only: which_mode
use mphys_die,         only: throw_mphys_error, incorrect_opt, std_msg
use casim_set_dependent_switches_mod, only: l_process, l_passivenumbers,       &
                                            l_passivenumbers_ice

implicit none

!-------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------

! Number of model levels
integer, intent(in) :: nlayers

! Whether the aerosol supplied to CASIM is prognostic, and so may be depleted
! by the activation. The UM only does this for the prognostic CASIM aerosol
! and for two way UKCA coupling.
logical, intent(in) :: l_process_aerosol

! Whether to cap the in-cloud droplet number. The UM only does this for the
! murk aerosol source.
logical, intent(in) :: l_cap_drop_number

! Cloud liquid mass after, and before, the non-CASIM physics [kg kg-1]
real(kind=real_umphys), intent(in) :: cloud_mass_post_qtbal(nlayers)
real(kind=real_umphys), intent(in) :: cloud_mass_pre_ap2(nlayers)

! Dry air density [kg m-3]
real(kind=real_umphys), intent(in) :: rho_col(nlayers)

! Temperature [K] and pressure [Pa]
real(kind=real_umphys), intent(in) :: t_col(nlayers)
real(kind=real_umphys), intent(in) :: p_col(nlayers)

! Liquid cloud fraction after, and before, the non-CASIM physics []
real(kind=real_umphys), intent(in) :: cf_liquid(nlayers)
real(kind=real_umphys), intent(in) :: cf_liquid_pre_ap2(nlayers)

! Vertical velocity used for the activation [m s-1]
real(kind=real_umphys), intent(in) :: w_tke(nlayers)

! Cloud droplet number [kg-1]
real(kind=real_umphys), intent(inout) :: cloudnumber_col(nlayers)

! Rain drop number [kg-1]
real(kind=real_umphys), intent(in) :: rainnumber_col(nlayers)

! Soluble aerosol mass [kg kg-1] and number [kg-1] in each mode
real(kind=real_umphys), intent(inout) :: aitken_sol_mass(nlayers)
real(kind=real_umphys), intent(inout) :: aitken_sol_number(nlayers)
real(kind=real_umphys), intent(inout) :: accum_sol_mass(nlayers)
real(kind=real_umphys), intent(inout) :: accum_sol_number(nlayers)
real(kind=real_umphys), intent(inout) :: coarse_sol_mass(nlayers)
real(kind=real_umphys), intent(inout) :: coarse_sol_number(nlayers)

! Coarse mode dust mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(inout) :: coarse_dust_mass(nlayers)
real(kind=real_umphys), intent(inout) :: coarse_dust_number(nlayers)

! Aerosol which is already activated in the liquid
real(kind=real_umphys), intent(inout) :: active_sol_liquid(nlayers)
real(kind=real_umphys), intent(inout) :: active_sol_number(nlayers)
real(kind=real_umphys), intent(inout) :: active_insol_liquid(nlayers)
real(kind=real_umphys), intent(in)    :: active_insol_ice(nlayers)
real(kind=real_umphys), intent(inout) :: active_insol_number(nlayers)

!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------

! Working cloud droplet number for the level being worked on [kg-1]
real(kind=real_umphys) :: w_cloud_number

! Change in the cloud mass and liquid cloud fraction over the non-CASIM physics
real(kind=real_umphys) :: delta_mass, delta_cfliq

! Rates of change of droplet number and of activated aerosol mass
real(kind=real_umphys) :: dnumber, dnumber_a, dnumber_d
real(kind=real_umphys) :: dmac, dmac1, dmac2
real(kind=real_umphys) :: dmad, dnac1, dnac2

! Diagnostic output from activate, not used here
real(kind=real_umphys) :: smax, ait_cdnc, accum_cdnc, tot_cdnc,                &
                          activated_cloud, activated_arg

! Whether the prognostic aerosol is to be depleted by the activation
logical :: l_deplete

integer :: k ! Loop counter over levels

!-------------------------------------------------------------------------------
! End of declarations and start of subroutine
!-------------------------------------------------------------------------------

l_deplete = l_process .and. l_process_aerosol

!------------------------------------------------------------------------------
! Main loop over levels to call the activate routine from within CASIM.
!------------------------------------------------------------------------------
do k = 1, nlayers

  ! Work out if cloud mass or liquid cloud fraction has changed.
  delta_mass  = cloud_mass_post_qtbal(k) - cloud_mass_pre_ap2(k)
  delta_cfliq = cf_liquid(k) - cf_liquid_pre_ap2(k)

  ! Set local working cloud number for this level.
  w_cloud_number = cloudnumber_col(k)

  ! set to zero - check all options are covered below.
  dnumber   = zero
  dnumber_a = zero
  dnumber_d = zero
  dmac      = zero
  dmad      = zero
  dnccn_all(:) = zero
  dmac_all(:)  = zero
  dnccnd_all(aero_index%i_coarse_dust) = zero
  dmad_all(aero_index%i_coarse_dust)   = zero

  if ( delta_mass > zero .and. cloud_mass_post_qtbal(k) > ql_tidy .and.        &
       cf_liquid(k) > cfliq_small ) then

    !---------------------------------------------------------------
    ! There has been an increase in liquid cloud due to non-CASIM
    ! physical processes in the model. Call the CASIM routine activate
    ! (which operates on one grid box) in order to work out what the
    ! new cloud number is based on the aerosol and the cloud mass
    ! added by the other processes.
    !---------------------------------------------------------------
    call activate( timestep, cloud_mass_post_qtbal(k), w_cloud_number,         &
                   w_tke(k), rho_col(k), dnumber, dmac, t_col(k), p_col(k),    &
                   cf_liquid(k), cf_liquid_pre_ap2(k),                         &
                   aerophys(k), aerochem(k), aeroact(k),                       &
                   dustphys(k), dustchem(k), dustliq(k),                       &
                   dnccn_all, dmac_all, dnumber_d, dmad,                       &
                   dnccnd_all, dmad_all,                                       &
                   smax, ait_cdnc, accum_cdnc, tot_cdnc, activated_arg,        &
                   activated_cloud )

    dnumber_a = dnumber

  else if ( cloud_mass_post_qtbal(k) < ql_tidy .or.                            &
            cf_liquid(k)             < cfliq_small ) then

    !---------------------------------------------------------------
    ! Non-CASIM processes have removed all of the liquid cloud or the
    ! liquid cloud fraction, so remove the cloud number.
    !---------------------------------------------------------------
    dnumber = -w_cloud_number * recip_timestep

    w_cloud_number = zero

    !============================
    ! aerosol processing
    !============================
    if ( l_deplete ) then

      dmac = -active_sol_liquid(k) * recip_timestep
      dmad = -active_insol_liquid(k) * recip_timestep

      if ( l_passivenumbers ) then
        dnumber_a = -active_sol_number(k) * recip_timestep
      else
        dnumber_a = dnumber
      end if

      if ( l_passivenumbers_ice ) then
        dnumber_d = -active_insol_number(k) * recip_timestep *                 &
            cloudnumber_col(k) / ( cloudnumber_col(k) + rainnumber_col(k) +    &
                                   epsilon(1.0_real_umphys) ) *                &
            active_insol_liquid(k) / ( active_insol_liquid(k) +                &
                                       active_insol_ice(k) +                   &
                                       epsilon(1.0_real_umphys) )
      else
        dnumber_d = dnumber
      end if

      dmad_all(aero_index%i_coarse_dust)   = dmad
      dnccnd_all(aero_index%i_coarse_dust) = dnumber_d

      call split_between_modes( k, dmac, dnumber_a, dmac2, dnac2 )

    end if ! l_deplete

  else

    if ( delta_cfliq < zero ) then

      ! dcfliq is negative!
      dnumber = w_cloud_number * delta_cfliq / cf_liquid_pre_ap2(k) *          &
                recip_timestep

      if ( l_deplete ) then

        dmac = ( ( active_sol_liquid(k) * delta_cfliq ) /                      &
                 cf_liquid_pre_ap2(k) ) * recip_timestep
        dnumber_a = dnumber

        ! N.B. the UM overwrites the accumulation mode part of the split with
        ! the whole of dmac and dnumber_a immediately after this point, which
        ! double counts the aerosol whenever which_mode has put any of it into
        ! the coarse mode. That is not reproduced here.
        call split_between_modes( k, dmac, dnumber_a, dmac2, dnac2 )

        ! do mass weighting of number and put back into coarse dust
        dnumber_d = ( active_insol_number(k) * active_insol_liquid(k) ) /      &
            ( active_insol_liquid(k) + active_insol_ice(k) +                   &
              epsilon(1.0_real_umphys) ) *                                     &
            ( delta_cfliq / cf_liquid_pre_ap2(k) ) * recip_timestep
        dmad = ( ( active_insol_liquid(k) * delta_cfliq ) /                    &
                 cf_liquid_pre_ap2(k) ) * recip_timestep

        dnccnd_all(aero_index%i_coarse_dust) = dnumber_d
        dmad_all(aero_index%i_coarse_dust)   = dmad

      end if ! l_deplete

    end if ! delta_cfliq < zero

  end if ! delta_mass > zero etc.

  ! The murk aerosol scheme places a cap on the in-cloud droplet number. The
  ! test is only made where there is some liquid cloud fraction to divide by;
  ! the UM does not guard against that division.
  if ( l_cap_drop_number .and. cf_liquid(k) > cfliq_small ) then
    if ( ( w_cloud_number + ( dnumber * timestep ) ) / cf_liquid(k) >          &
         max_drop_casim ) then
      dnumber = ( max_drop_casim * cf_liquid(k) - w_cloud_number ) *           &
                recip_timestep
    end if
  end if

  ! Update cloud number in the column
  w_cloud_number = max( w_cloud_number +                                       &
                        ( ( dnumber + dnumber_d ) * timestep ), zero )
  cloudnumber_col(k) = w_cloud_number

  !----------------------------------------------------------------
  ! Move the activated aerosol out of the interstitial aerosol and
  ! into the activated aerosol.
  !----------------------------------------------------------------
  if ( l_deplete ) then

    if ( l_passivenumbers ) active_sol_number(k) =                             &
        max( zero, active_sol_number(k) + dnumber_a*timestep )
    if ( l_passivenumbers_ice ) active_insol_number(k) =                       &
        max( zero, active_insol_number(k) + dnumber_d*timestep )

    if ( aero_index%i_aitken > 0 ) then
      dmac_all(aero_index%i_aitken) =                                          &
          min( aitken_sol_mass(k),                                             &
               dmac_all(aero_index%i_aitken)*timestep ) * recip_timestep
      active_sol_liquid(k) = max( zero, active_sol_liquid(k) +                 &
          dmac_all(aero_index%i_aitken)*timestep )
      aitken_sol_mass(k) = max( zero, aitken_sol_mass(k) -                     &
          dmac_all(aero_index%i_aitken)*timestep )
      aitken_sol_number(k) = max( zero, aitken_sol_number(k) -                 &
          dnccn_all(aero_index%i_aitken)*timestep )
    end if ! aero_index%i_aitken > 0

    if ( aero_index%i_accum > 0 ) then
      dmac_all(aero_index%i_accum) =                                           &
          min( accum_sol_mass(k),                                              &
               dmac_all(aero_index%i_accum)*timestep ) * recip_timestep
      active_sol_liquid(k) = max( zero, active_sol_liquid(k) +                 &
          dmac_all(aero_index%i_accum)*timestep )
      accum_sol_mass(k) = max( zero, accum_sol_mass(k) -                       &
          dmac_all(aero_index%i_accum)*timestep )
      accum_sol_number(k) = max( zero, accum_sol_number(k) -                   &
          dnccn_all(aero_index%i_accum)*timestep )
      if ( accum_sol_mass(k)   < aeromass_small .or.                           &
           accum_sol_number(k) < aeronumber_small ) then
        accum_sol_mass(k)   = zero
        accum_sol_number(k) = zero
      end if
    end if ! aero_index%i_accum > 0

    if ( aero_index%i_coarse > 0 ) then
      dmac_all(aero_index%i_coarse) =                                          &
          min( coarse_sol_mass(k),                                             &
               dmac_all(aero_index%i_coarse)*timestep ) * recip_timestep
      active_sol_liquid(k) = max( zero, active_sol_liquid(k) +                 &
          dmac_all(aero_index%i_coarse)*timestep )
      coarse_sol_mass(k) = max( zero, coarse_sol_mass(k) -                     &
          dmac_all(aero_index%i_coarse)*timestep )
      coarse_sol_number(k) = max( zero, coarse_sol_number(k) -                 &
          dnccn_all(aero_index%i_coarse)*timestep )
      if ( coarse_sol_mass(k)   < aeromass_small .or.                          &
           coarse_sol_number(k) < aeronumber_small ) then
        coarse_sol_mass(k)   = zero
        coarse_sol_number(k) = zero
      end if
    end if ! aero_index%i_coarse > 0

    if ( .not. l_warm ) then
      ! We may have some dust in the liquid...
      dmad_all(aero_index%i_coarse_dust) =                                     &
          min( coarse_dust_mass(k),                                            &
               dmad_all(aero_index%i_coarse_dust)*timestep ) * recip_timestep
      active_insol_liquid(k) = max( zero, active_insol_liquid(k) +             &
          dmad_all(aero_index%i_coarse_dust)*timestep )
      coarse_dust_mass(k) = max( zero, coarse_dust_mass(k) -                   &
          dmad_all(aero_index%i_coarse_dust)*timestep )
      coarse_dust_number(k) = max( zero, coarse_dust_number(k) -               &
          dnccnd_all(aero_index%i_coarse_dust)*timestep )
    end if ! .not. l_warm

  end if ! l_deplete

end do ! k

contains

!-------------------------------------------------------------------------------
! Split the activated soluble aerosol between the accumulation and coarse modes
!-------------------------------------------------------------------------------

subroutine split_between_modes( level, dmass, dnumber_sol, dmass2, dnum2 )

implicit none

integer,                intent(in)    :: level
real(kind=real_umphys), intent(in)    :: dmass
real(kind=real_umphys), intent(inout) :: dnumber_sol
real(kind=real_umphys), intent(out)   :: dmass2, dnum2

real(kind=real_umphys) :: dmass1, dnum1

dmass1 = zero
dmass2 = zero
dnum1  = zero
dnum2  = zero

if ( aero_index%i_accum > 0 .and. aero_index%i_coarse > 0 ) then

  ! We have both accumulation and coarse modes
  if ( dnumber_sol*dmass <= zero ) then
    dnumber_sol = dmass / 1.0e-18_real_umphys * recip_timestep
  end if

  if ( l_bypass_which_mode ) then
    ! Don't use which_mode - just transfer all to either accum or coarse mode
    if ( iopt_which_mode == 1 ) then
      ! All to accum
      dmass1 = dmass
      dmass2 = zero
      dnum1  = dnumber_sol
      dnum2  = zero
    else if ( iopt_which_mode == 2 ) then
      ! All to coarse
      dmass1 = zero
      dmass2 = dmass
      dnum1  = zero
      dnum2  = dnumber_sol
    else
      write(std_msg, '(A)') "incorrect iopt_which_mode option selected"
      call throw_mphys_error( incorrect_opt, ModuleName, std_msg )
    end if
  else
    call which_mode( dmass, dnumber_sol,                                       &
                     aerophys(level)%rd(aero_index%i_accum),                   &
                     aerophys(level)%rd(aero_index%i_coarse),                  &
                     aerochem(level)%density(aero_index%i_accum),              &
                     aerophys(level)%sigma(aero_index%i_accum),                &
                     dmass1, dmass2, dnum1, dnum2 )
  end if

  ! put the split back into the accumulation and coarse modes
  dmac_all(aero_index%i_accum)   = dmass1
  dnccn_all(aero_index%i_accum)  = dnum1
  dmac_all(aero_index%i_coarse)  = dmass2
  dnccn_all(aero_index%i_coarse) = dnum2

else

  ! Only one of the two modes is available, so it takes all of the aerosol
  if ( aero_index%i_accum > 0 ) then
    dmac_all(aero_index%i_accum)  = dmass
    dnccn_all(aero_index%i_accum) = dnumber_sol
  end if

  if ( aero_index%i_coarse > 0 ) then
    dmac_all(aero_index%i_coarse)  = dmass
    dnccn_all(aero_index%i_coarse) = dnumber_sol
  end if

end if

end subroutine split_between_modes

end subroutine activate_column

end module casim_activation_in_um_mod
