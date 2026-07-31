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
! Convert the aerosol supplied by the host model into the mass and number
! fields expected by CASIM.
!
! These routines are single column versions of the UM routines of the same
! names. The UM works on the (i,j,k) fields of a whole domain, whereas an
! LFRic kernel only ever sees one column, so the horizontal dimensions have
! been dropped and the arrays are simply indexed by level.

module aerosol_extract_convert_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private ::                                        &
  ModuleName='AEROSOL_EXTRACT_CONVERT_MOD'

private
public :: aerosol_extract_convert, aerosol_extract_convert_murk,               &
          aerosol_extract_convert_ft

contains

!-------------------------------------------------------------------------------
! Aerosol from the GLOMAP modal aerosol scheme
!-------------------------------------------------------------------------------
! This routine returns the soluble aitken, accumulation and coarse, and the
! insoluble accumulation and coarse, mass and number.
!
! The UM reads the GLOMAP modes out of its single tracer_ukca array using the
! mmr_index_um and nmr_index_um lookups. LFRic holds each mode component in a
! separately named field instead, so the caller is responsible for gathering
! them into mode_mmr and mode_nmr; the loops over modes and components below
! are otherwise unchanged from the UM.

subroutine aerosol_extract_convert( nlayers, n_cpt, p_theta_levels, t_work,    &
                                    rho_casim, mode_mmr, mode_nmr,             &
                                    aitken_sol_mass,  aitken_sol_number,       &
                                    accum_sol_mass,   accum_sol_number,        &
                                    coarse_sol_mass,  coarse_sol_number,       &
                                    accum_dust_mass,  accum_dust_number,       &
                                    coarse_dust_mass, coarse_dust_number,      &
                                    aitken_sol_bk, accum_sol_bk,               &
                                    coarse_sol_bk )

use ukca_config_specification_mod, only: glomap_variables

use ukca_mode_setup,       only: nmodes,                                       &
                                 mode_ait_sol, mode_acc_sol,                   &
                                 mode_cor_sol, mode_acc_insol,                 &
                                 mode_cor_insol

use ukca_constants,        only: mmw
use water_constants_mod,   only: rho_water
use chemistry_constants_mod, only: boltzmann
use missing_data_mod,      only: rmdi

implicit none

!-------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------

! Number of model levels
integer, intent(in) :: nlayers

! Number of GLOMAP components held by the caller. This must be at least as
! large as glomap_variables%ncp.
integer, intent(in) :: n_cpt

! Pressure at layer centres [Pa]
real(kind=real_umphys), intent(in) :: p_theta_levels(nlayers)

! Local working temperature [K]
real(kind=real_umphys), intent(in) :: t_work(nlayers)

! Air density for CASIM [kg m-3]
real(kind=real_umphys), intent(in) :: rho_casim(nlayers)

! Mass mixing ratio of each component in each GLOMAP mode [kg kg-1]
real(kind=real_umphys), intent(in) :: mode_mmr(nlayers,nmodes,n_cpt)

! Number mixing ratio of each GLOMAP mode [molecules per molecule of air]
real(kind=real_umphys), intent(in) :: mode_nmr(nlayers,nmodes)

! Aitken mode soluble mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: aitken_sol_mass(nlayers)
real(kind=real_umphys), intent(out) :: aitken_sol_number(nlayers)

! Accumulation mode soluble mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: accum_sol_mass(nlayers)
real(kind=real_umphys), intent(out) :: accum_sol_number(nlayers)

! Coarse mode soluble mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: coarse_sol_mass(nlayers)
real(kind=real_umphys), intent(out) :: coarse_sol_number(nlayers)

! Accumulation mode dust mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: accum_dust_mass(nlayers)
real(kind=real_umphys), intent(out) :: accum_dust_number(nlayers)

! Coarse mode dust mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: coarse_dust_mass(nlayers)
real(kind=real_umphys), intent(out) :: coarse_dust_number(nlayers)

! Activation parameter Bk for each soluble mode []
real(kind=real_umphys), intent(out) :: aitken_sol_bk(nlayers)
real(kind=real_umphys), intent(out) :: accum_sol_bk(nlayers)
real(kind=real_umphys), intent(out) :: coarse_sol_bk(nlayers)

!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
logical, pointer :: component(:,:)
real,    pointer :: mm(:)
logical, pointer :: mode(:)
integer, pointer :: ncp
real,    pointer :: no_ions(:)
real,    pointer :: rhocomp(:)

! Volume of each soluble mode, used to weight the Bk values [m3 kg-1]
real(kind=real_umphys) :: aitken_sol_volume(nlayers)
real(kind=real_umphys) :: accum_sol_volume(nlayers)
real(kind=real_umphys) :: coarse_sol_volume(nlayers)

! Number density of air [m-3]
real(kind=real_umphys) :: aird

integer :: k              ! Loop counter over levels
integer :: imode, icp     ! UKCA counters

real(kind=real_umphys), parameter :: real_eps = epsilon(1.0_real_umphys)

!-------------------------------------------------------------------------------
! End of declarations and start of subroutine
!-------------------------------------------------------------------------------

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mm          => glomap_variables%mm
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
no_ions     => glomap_variables%no_ions
rhocomp     => glomap_variables%rhocomp

! Sum up component masses in each mode for total mass and copy number to
! inputs for CASIM.

! Initialise modal mass arrays to zero (will contain sum of mass over all
! components)
do k = 1, nlayers
  aitken_sol_mass(k)   = 0.0_real_umphys
  accum_sol_mass(k)    = 0.0_real_umphys
  coarse_sol_mass(k)   = 0.0_real_umphys
  accum_dust_mass(k)   = 0.0_real_umphys
  coarse_dust_mass(k)  = 0.0_real_umphys
  aitken_sol_bk(k)     = 0.0_real_umphys
  accum_sol_bk(k)      = 0.0_real_umphys
  coarse_sol_bk(k)     = 0.0_real_umphys
  aitken_sol_volume(k) = 0.0_real_umphys
  accum_sol_volume(k)  = 0.0_real_umphys
  coarse_sol_volume(k) = 0.0_real_umphys
end do

! The number is only set for the modes which are active, so make sure the
! remainder are left in a sensible state.
do k = 1, nlayers
  aitken_sol_number(k)  = 0.0_real_umphys
  accum_sol_number(k)   = 0.0_real_umphys
  coarse_sol_number(k)  = 0.0_real_umphys
  accum_dust_number(k)  = 0.0_real_umphys
  coarse_dust_number(k) = 0.0_real_umphys
end do

! Include modes 2 to 6 (Aitsol, accsol, corsol, Aitins, accins, corins)
do imode = mode_ait_sol, mode_cor_insol
  if (mode(imode)) then
    do icp = 1, ncp
      if (component(imode,icp)) then
        ! The code in the CASIM activation scheme that is replaced by this
        ! interface:
        ! start:
        ! Bk=chem%vantHoff(i)*Mw*chem%density(i)/(chem%massMole(i)*rho_water)
        ! no_ions, rho_comp from ukca_mode_setup, massMole is 0.132
        ! which is ammonium sulphate (18*2+96). Density is set to 1777 in
        ! mphys_constants.F90.
        ! end
        ! The Bk value here does not assume ammonium sulphate but is a volume
        ! weighted average of the Bk values that would come from the UKCA
        ! components

        ! For each active mode calculate mass to transfer to CASIM, and Bk
        ! numerator and denominator
        select case (imode)

        case (mode_ait_sol)
          do k = 1, nlayers
            aitken_sol_mass(k)   = aitken_sol_mass(k) + mode_mmr(k,imode,icp)
            aitken_sol_bk(k)     = aitken_sol_bk(k) + no_ions(icp) *           &
                                   mmw * mode_mmr(k,imode,icp) /               &
                                   (rho_water*mm(icp))
            aitken_sol_volume(k) = aitken_sol_volume(k) +                      &
                                   mode_mmr(k,imode,icp) / rhocomp(icp)
          end do

        case (mode_acc_sol)
          do k = 1, nlayers
            accum_sol_mass(k)    = accum_sol_mass(k) + mode_mmr(k,imode,icp)
            accum_sol_bk(k)      = accum_sol_bk(k) + no_ions(icp) *            &
                                   mmw * mode_mmr(k,imode,icp) /               &
                                   (rho_water*mm(icp))
            accum_sol_volume(k)  = accum_sol_volume(k) +                       &
                                   mode_mmr(k,imode,icp) / rhocomp(icp)
          end do

        case (mode_cor_sol)
          do k = 1, nlayers
            coarse_sol_mass(k)   = coarse_sol_mass(k) + mode_mmr(k,imode,icp)
            coarse_sol_bk(k)     = coarse_sol_bk(k) + no_ions(icp) *           &
                                   mmw * mode_mmr(k,imode,icp) /               &
                                   (rho_water*mm(icp))
            coarse_sol_volume(k) = coarse_sol_volume(k) +                      &
                                   mode_mmr(k,imode,icp) / rhocomp(icp)
          end do

        case (mode_acc_insol)
          do k = 1, nlayers
            accum_dust_mass(k)  = accum_dust_mass(k) + mode_mmr(k,imode,icp)
          end do

        case (mode_cor_insol)
          do k = 1, nlayers
            coarse_dust_mass(k) = coarse_dust_mass(k) + mode_mmr(k,imode,icp)
          end do

        end select
      end if ! component(imode,icp)
    end do ! loop over icp

    ! For each active mode calculate number to transfer to CASIM.
    ! aird is required to convert number mixing ratio (from UKCA) to number
    ! density (CASIM uses number per kg hence the division by rho_casim).
    select case (imode)

    case (mode_ait_sol)
      do k = 1, nlayers
        aird = p_theta_levels(k) / ( t_work(k)*boltzmann )
        aitken_sol_number(k)  = mode_nmr(k,imode)*aird/rho_casim(k)
      end do

    case (mode_acc_sol)
      do k = 1, nlayers
        aird = p_theta_levels(k) / ( t_work(k)*boltzmann )
        accum_sol_number(k)   = mode_nmr(k,imode)*aird/rho_casim(k)
      end do

    case (mode_cor_sol)
      do k = 1, nlayers
        aird = p_theta_levels(k) / ( t_work(k)*boltzmann )
        coarse_sol_number(k)  = mode_nmr(k,imode)*aird/rho_casim(k)
      end do

    case (mode_acc_insol)
      do k = 1, nlayers
        aird = p_theta_levels(k) / ( t_work(k)*boltzmann )
        accum_dust_number(k)  = mode_nmr(k,imode)*aird/rho_casim(k)
      end do

    case (mode_cor_insol)
      do k = 1, nlayers
        aird = p_theta_levels(k) / ( t_work(k)*boltzmann )
        coarse_dust_number(k) = mode_nmr(k,imode)*aird/rho_casim(k)
      end do

    end select
  end if ! mode(imode)
end do ! loop over imode

! Calculate activation parameters Bk for soluble Aitken, accumulation and
! coarse modes. These modes are always used in UKCA, even if not always
! used in CASIM.
! If there is no aerosol in a mode then CASIM will not carry out activation
! in that mode. In order to avoid carrying out spurious divides by zero here,
! we set the Bk to missing data if the volume is zero.
! N.B. the accumulation and coarse mode tests are on the Bk rather than the
! volume. This looks like a mistake but is retained to match the UM.
do k = 1, nlayers
  if (abs(aitken_sol_volume(k)) < real_eps) then
    aitken_sol_bk(k) = rmdi
  else
    aitken_sol_bk(k) = aitken_sol_bk(k)/aitken_sol_volume(k)
  end if

  if (abs(accum_sol_bk(k)) < real_eps) then
    accum_sol_bk(k) = rmdi
  else
    accum_sol_bk(k) = accum_sol_bk(k)/accum_sol_volume(k)
  end if

  if (abs(coarse_sol_bk(k)) < real_eps) then
    coarse_sol_bk(k) = rmdi
  else
    coarse_sol_bk(k) = coarse_sol_bk(k)/coarse_sol_volume(k)
  end if
end do

end subroutine aerosol_extract_convert

!-------------------------------------------------------------------------------
! Aerosol from the murk aerosol field
!-------------------------------------------------------------------------------
! All of the murk mass is placed in the soluble accumulation mode and the
! number is derived from it, so all of the other modes are returned as zero.

subroutine aerosol_extract_convert_murk( nlayers, rho_casim, aerosol,          &
                                         aitken_sol_mass,  aitken_sol_number,  &
                                         accum_sol_mass,   accum_sol_number,   &
                                         coarse_sol_mass,  coarse_sol_number,  &
                                         accum_dust_mass,  accum_dust_number,  &
                                         coarse_dust_mass, coarse_dust_number )

use mphys_constants_mod,  only: n0_murk, m0_murk, mu_g_to_kg, min_n_aer
use lsp_autoc_consts_mod, only: power_murk

implicit none

!-------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------

! Number of model levels
integer, intent(in) :: nlayers

! Air density for CASIM [kg m-3]
real(kind=real_umphys), intent(in) :: rho_casim(nlayers)

! Murk aerosol [ug kg-1]
real(kind=real_umphys), intent(in) :: aerosol(nlayers)

! Aitken mode soluble mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: aitken_sol_mass(nlayers)
real(kind=real_umphys), intent(out) :: aitken_sol_number(nlayers)

! Accumulation mode soluble mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: accum_sol_mass(nlayers)
real(kind=real_umphys), intent(out) :: accum_sol_number(nlayers)

! Coarse mode soluble mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: coarse_sol_mass(nlayers)
real(kind=real_umphys), intent(out) :: coarse_sol_number(nlayers)

! Accumulation mode dust mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: accum_dust_mass(nlayers)
real(kind=real_umphys), intent(out) :: accum_dust_number(nlayers)

! Coarse mode dust mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: coarse_dust_mass(nlayers)
real(kind=real_umphys), intent(out) :: coarse_dust_number(nlayers)

!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------

real(kind=real_umphys) :: n_aer ! Local working aerosol number

real(kind=real_umphys), parameter :: zero = 0.0_real_umphys

integer :: k ! Loop counter over levels

!-------------------------------------------------------------------------------
! End of declarations and start of subroutine
!-------------------------------------------------------------------------------

do k = 1, nlayers

  ! Follow lsp_taper_ndrop to get aerosol number from murk mass:
  ! use Haywood et al (2008) formulae to get droplet number

  n_aer = max( aerosol(k) / m0_murk*mu_g_to_kg, min_n_aer)
  ! mu_g_to_kg converts from ug/kg to kg/kg

  !-----------------------------------------------
  ! Set aerosol mass to accumulation mode
  !-----------------------------------------------
  accum_sol_mass(k) = aerosol(k)*mu_g_to_kg

  !-----------------------------------------------
  ! Calculation of the aerosol number
  !-----------------------------------------------
  n_aer = n0_murk * n_aer ** power_murk
  accum_sol_number(k) = n_aer / rho_casim(k)

  !----------------------------------------------
  ! Set all other aerosol species to zero
  !----------------------------------------------
  aitken_sol_mass(k)    = zero
  aitken_sol_number(k)  = zero
  coarse_sol_mass(k)    = zero
  coarse_sol_number(k)  = zero
  accum_dust_mass(k)    = zero
  accum_dust_number(k)  = zero
  coarse_dust_mass(k)   = zero
  coarse_dust_number(k) = zero
end do

end subroutine aerosol_extract_convert_murk

!-------------------------------------------------------------------------------
! Fixed aerosol, or aerosol from the prognostic CASIM aerosol fields
!-------------------------------------------------------------------------------
! The UM reads the prognostic aerosol out of its free tracer array. LFRic
! holds each species in a separately named field, so they are passed in
! individually here.

subroutine aerosol_extract_convert_ft( nlayers,                                &
                                       aitken_sol_mass_in,                     &
                                       aitken_sol_number_in,                   &
                                       accum_sol_mass_in,                      &
                                       accum_sol_number_in,                    &
                                       coarse_sol_mass_in,                     &
                                       coarse_sol_number_in,                   &
                                       accum_dust_mass_in,                     &
                                       accum_dust_number_in,                   &
                                       coarse_dust_mass_in,                    &
                                       coarse_dust_number_in,                  &
                                       aitken_sol_mass,  aitken_sol_number,    &
                                       accum_sol_mass,   accum_sol_number,     &
                                       coarse_sol_mass,  coarse_sol_number,    &
                                       accum_dust_mass,  accum_dust_number,    &
                                       coarse_dust_mass, coarse_dust_number )

use casim_switches,   only: l_fix_aerosol, l_tracer_aerosol
use mphys_inputs_mod, only: fixed_accum_sol_mass, fixed_accum_sol_number

implicit none

!-------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------

! Number of model levels
integer, intent(in) :: nlayers

! Prognostic aerosol mass [kg kg-1] and number [kg-1] fields
real(kind=real_umphys), intent(in) :: aitken_sol_mass_in(nlayers)
real(kind=real_umphys), intent(in) :: aitken_sol_number_in(nlayers)
real(kind=real_umphys), intent(in) :: accum_sol_mass_in(nlayers)
real(kind=real_umphys), intent(in) :: accum_sol_number_in(nlayers)
real(kind=real_umphys), intent(in) :: coarse_sol_mass_in(nlayers)
real(kind=real_umphys), intent(in) :: coarse_sol_number_in(nlayers)
real(kind=real_umphys), intent(in) :: accum_dust_mass_in(nlayers)
real(kind=real_umphys), intent(in) :: accum_dust_number_in(nlayers)
real(kind=real_umphys), intent(in) :: coarse_dust_mass_in(nlayers)
real(kind=real_umphys), intent(in) :: coarse_dust_number_in(nlayers)

! Aitken mode soluble mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: aitken_sol_mass(nlayers)
real(kind=real_umphys), intent(out) :: aitken_sol_number(nlayers)

! Accumulation mode soluble mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: accum_sol_mass(nlayers)
real(kind=real_umphys), intent(out) :: accum_sol_number(nlayers)

! Coarse mode soluble mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: coarse_sol_mass(nlayers)
real(kind=real_umphys), intent(out) :: coarse_sol_number(nlayers)

! Accumulation mode dust mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: accum_dust_mass(nlayers)
real(kind=real_umphys), intent(out) :: accum_dust_number(nlayers)

! Coarse mode dust mass [kg kg-1] and number [kg-1]
real(kind=real_umphys), intent(out) :: coarse_dust_mass(nlayers)
real(kind=real_umphys), intent(out) :: coarse_dust_number(nlayers)

!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------

real(kind=real_umphys), parameter :: zero = 0.0_real_umphys

integer :: k ! Loop counter over levels

!-------------------------------------------------------------------------------
! End of declarations and start of subroutine
!-------------------------------------------------------------------------------

if ( l_fix_aerosol ) then
  ! No transport of aerosol, set CASIM aerosol accumulation mode to a
  ! fixed value, all others set to 0.0

  do k = 1, nlayers

    aitken_sol_mass(k)    = zero
    aitken_sol_number(k)  = zero
    coarse_sol_mass(k)    = zero
    coarse_sol_number(k)  = zero
    accum_dust_mass(k)    = zero
    accum_dust_number(k)  = zero
    coarse_dust_mass(k)   = zero
    coarse_dust_number(k) = zero

    accum_sol_mass(k)     = fixed_accum_sol_mass
    accum_sol_number(k)   = fixed_accum_sol_number

  end do ! k

else if ( l_tracer_aerosol ) then

  do k = 1, nlayers

    ! Set soluble aerosol dependent on namelist setting
    aitken_sol_mass(k)    = aitken_sol_mass_in(k)
    aitken_sol_number(k)  = aitken_sol_number_in(k)

    accum_sol_mass(k)     = accum_sol_mass_in(k)
    accum_sol_number(k)   = accum_sol_number_in(k)

    coarse_sol_mass(k)    = coarse_sol_mass_in(k)
    coarse_sol_number(k)  = coarse_sol_number_in(k)

    ! Set the insoluble aerosol dependent on namelist
    accum_dust_mass(k)    = accum_dust_mass_in(k)
    accum_dust_number(k)  = accum_dust_number_in(k)

    coarse_dust_mass(k)   = coarse_dust_mass_in(k)
    coarse_dust_number(k) = coarse_dust_number_in(k)

  end do ! k

end if ! l_fix_aerosol / l_tracer_aerosol

end subroutine aerosol_extract_convert_ft

end module aerosol_extract_convert_mod
