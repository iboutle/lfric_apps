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
! Return the aerosol which CASIM has taken out of, or given back to, the
! interstitial aerosol to whichever aerosol scheme supplied it.
!
! These routines are single column versions of the UM routines
! aerosol_convert_return_solinsol and aerosol_convert_return_ft. The UM works
! on the (i,j,k) fields of a whole domain, whereas an LFRic kernel only ever
! sees one column, so the horizontal dimensions have been dropped and the
! arrays are simply indexed by level.
!
! The aerosol is held at or above zero on the way back. The UM does not need
! to do this here because its callers have already limited the depletion to
! what the aerosol holds, but it does the same clipping inline in
! activate_column_2way_ukca.

module aerosol_convert_return_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private ::                                        &
  ModuleName='AEROSOL_CONVERT_RETURN_MOD'

private
public :: aerosol_convert_return, aerosol_convert_return_ft

contains

!-------------------------------------------------------------------------------
! Apply the CASIM aerosol increments to the GLOMAP modal aerosol
!-------------------------------------------------------------------------------
! The increments are supplied in the units CASIM works in, that is a mass
! mixing ratio and a number per kilogram of air. The GLOMAP mass is also a mass
! mixing ratio, but the GLOMAP number is a number mixing ratio, so the number
! increment is converted using the number density of air in the same way as
! aerosol_extract_convert converts it on the way in.
!
! As in the UM this routine is written for the simplified soluble/insoluble
! GLOMAP setup, in which each mode carries a single component: the soluble
! modes are all sulphate and the insoluble modes are all dust. The mass taken
! from a mode is therefore returned to the sulphate component of the soluble
! modes and to the dust component of the insoluble modes.
!
! The UM reads and writes the GLOMAP modes straight out of its single
! tracer_ukca array using the mmr_index_um and nmr_index_um lookups. LFRic
! holds each mode component in a separately named field instead, so the caller
! is responsible for gathering them into mode_mmr and mode_nmr and for
! scattering the updated values back afterwards.

subroutine aerosol_convert_return( nlayers, n_cpt, p_theta_levels, t_work,     &
                                   rho_casim,                                  &
                                   d_aitken_sol_mass,  d_aitken_sol_number,    &
                                   d_accum_sol_mass,   d_accum_sol_number,     &
                                   d_coarse_sol_mass,  d_coarse_sol_number,    &
                                   d_accum_dust_mass,  d_accum_dust_number,    &
                                   d_coarse_dust_mass, d_coarse_dust_number,   &
                                   mode_mmr, mode_nmr )

use ukca_config_specification_mod, only: glomap_variables

use ukca_mode_setup,       only: nmodes,                                       &
                                 mode_ait_sol, mode_acc_sol,                   &
                                 mode_cor_sol, mode_acc_insol,                 &
                                 mode_cor_insol, cp_su, cp_du

use chemistry_constants_mod, only: boltzmann

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

! Change in the soluble mass [kg kg-1] and number [kg-1] of each mode
real(kind=real_umphys), intent(in) :: d_aitken_sol_mass(nlayers)
real(kind=real_umphys), intent(in) :: d_aitken_sol_number(nlayers)
real(kind=real_umphys), intent(in) :: d_accum_sol_mass(nlayers)
real(kind=real_umphys), intent(in) :: d_accum_sol_number(nlayers)
real(kind=real_umphys), intent(in) :: d_coarse_sol_mass(nlayers)
real(kind=real_umphys), intent(in) :: d_coarse_sol_number(nlayers)

! Change in the dust mass [kg kg-1] and number [kg-1] of each mode
real(kind=real_umphys), intent(in) :: d_accum_dust_mass(nlayers)
real(kind=real_umphys), intent(in) :: d_accum_dust_number(nlayers)
real(kind=real_umphys), intent(in) :: d_coarse_dust_mass(nlayers)
real(kind=real_umphys), intent(in) :: d_coarse_dust_number(nlayers)

! Mass mixing ratio of each component in each GLOMAP mode [kg kg-1]
real(kind=real_umphys), intent(inout) :: mode_mmr(nlayers,nmodes,n_cpt)

! Number mixing ratio of each GLOMAP mode [molecules per molecule of air]
real(kind=real_umphys), intent(inout) :: mode_nmr(nlayers,nmodes)

!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
logical, pointer :: component(:,:)
logical, pointer :: mode(:)
logical, pointer :: soluble(:)

! The component each mode returns its mass to
integer :: this_cp(nmodes)

! Number density of air [m-3]
real(kind=real_umphys) :: aird(nlayers)

integer :: k          ! Loop counter over levels
integer :: imode, icp ! UKCA counters

real(kind=real_umphys), parameter :: zero = 0.0_real_umphys

!-------------------------------------------------------------------------------
! End of declarations and start of subroutine
!-------------------------------------------------------------------------------

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component => glomap_variables%component
mode      => glomap_variables%mode
soluble   => glomap_variables%soluble

! Split the target component by solubility (as in SOL/INSOL). Currently set to
! sulphate (cp_su) for the soluble modes and dust (cp_du) for the insoluble
! modes.
do imode = 1, nmodes
  if (soluble(imode)) then
    this_cp(imode) = cp_su
  else
    this_cp(imode) = cp_du
  end if
end do

! aird is required to convert the number per kilogram used by CASIM back into
! the number mixing ratio used by UKCA.
do k = 1, nlayers
  aird(k) = p_theta_levels(k) / ( t_work(k)*boltzmann )
end do

! Include modes 2 to 7 (Aitsol, accsol, corsol, Aitins, accins, corins)
do imode = mode_ait_sol, mode_cor_insol
  if (mode(imode)) then

    ! For each active mode return the mass to the UKCA component, if that
    ! component is active. The Aitken insoluble mode has no increment to
    ! return, so it is absent from the select constructs below.
    icp = this_cp(imode)
    if (component(imode,icp)) then
      select case (imode)

      case (mode_ait_sol)
        do k = 1, nlayers
          mode_mmr(k,imode,icp) = max( zero,                                   &
                                       mode_mmr(k,imode,icp) +                 &
                                       d_aitken_sol_mass(k) )
        end do

      case (mode_acc_sol)
        do k = 1, nlayers
          mode_mmr(k,imode,icp) = max( zero,                                   &
                                       mode_mmr(k,imode,icp) +                 &
                                       d_accum_sol_mass(k) )
        end do

      case (mode_cor_sol)
        do k = 1, nlayers
          mode_mmr(k,imode,icp) = max( zero,                                   &
                                       mode_mmr(k,imode,icp) +                 &
                                       d_coarse_sol_mass(k) )
        end do

      case (mode_acc_insol)
        do k = 1, nlayers
          mode_mmr(k,imode,icp) = max( zero,                                   &
                                       mode_mmr(k,imode,icp) +                 &
                                       d_accum_dust_mass(k) )
        end do

      case (mode_cor_insol)
        do k = 1, nlayers
          mode_mmr(k,imode,icp) = max( zero,                                   &
                                       mode_mmr(k,imode,icp) +                 &
                                       d_coarse_dust_mass(k) )
        end do

      end select
    end if ! component(imode,icp)

    ! For each active mode return the number to UKCA
    select case (imode)

    case (mode_ait_sol)
      do k = 1, nlayers
        mode_nmr(k,imode) = max( zero, mode_nmr(k,imode) +                     &
            d_aitken_sol_number(k)*rho_casim(k)/aird(k) )
      end do

    case (mode_acc_sol)
      do k = 1, nlayers
        mode_nmr(k,imode) = max( zero, mode_nmr(k,imode) +                     &
            d_accum_sol_number(k)*rho_casim(k)/aird(k) )
      end do

    case (mode_cor_sol)
      do k = 1, nlayers
        mode_nmr(k,imode) = max( zero, mode_nmr(k,imode) +                     &
            d_coarse_sol_number(k)*rho_casim(k)/aird(k) )
      end do

    case (mode_acc_insol)
      do k = 1, nlayers
        mode_nmr(k,imode) = max( zero, mode_nmr(k,imode) +                     &
            d_accum_dust_number(k)*rho_casim(k)/aird(k) )
      end do

    case (mode_cor_insol)
      do k = 1, nlayers
        mode_nmr(k,imode) = max( zero, mode_nmr(k,imode) +                     &
            d_coarse_dust_number(k)*rho_casim(k)/aird(k) )
      end do

    end select
  end if ! mode(imode)
end do ! loop over imode

end subroutine aerosol_convert_return

!-------------------------------------------------------------------------------
! Apply the CASIM aerosol increments to the prognostic CASIM aerosol
!-------------------------------------------------------------------------------
! The prognostic aerosol is held in the units CASIM works in, so there is
! nothing to convert and the increments are simply added on. The UM passes the
! pressure, temperature and density in here as well, but only uses them in a
! commented out conversion which is not wanted while the aerosol is not shared
! with UKCA, so they are left out.
!
! The UM reads and writes the prognostic aerosol straight out of its free
! tracer array using the i_AitkenSolMass and similar lookups. LFRic holds each
! species in a separately named field instead, so the caller is responsible for
! gathering them into the arrays below and for scattering the updated values
! back afterwards.

subroutine aerosol_convert_return_ft( nlayers,                                 &
                                   d_aitken_sol_mass,  d_aitken_sol_number,    &
                                   d_accum_sol_mass,   d_accum_sol_number,     &
                                   d_coarse_sol_mass,  d_coarse_sol_number,    &
                                   d_accum_dust_mass,  d_accum_dust_number,    &
                                   d_coarse_dust_mass, d_coarse_dust_number,   &
                                   aitken_sol_mass,  aitken_sol_number,        &
                                   accum_sol_mass,   accum_sol_number,         &
                                   coarse_sol_mass,  coarse_sol_number,        &
                                   accum_dust_mass,  accum_dust_number,        &
                                   coarse_dust_mass, coarse_dust_number )

use casim_switches, only: l_tracer_aerosol

implicit none

!-------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------

! Number of model levels
integer, intent(in) :: nlayers

! Change in the soluble mass [kg kg-1] and number [kg-1] of each mode
real(kind=real_umphys), intent(in) :: d_aitken_sol_mass(nlayers)
real(kind=real_umphys), intent(in) :: d_aitken_sol_number(nlayers)
real(kind=real_umphys), intent(in) :: d_accum_sol_mass(nlayers)
real(kind=real_umphys), intent(in) :: d_accum_sol_number(nlayers)
real(kind=real_umphys), intent(in) :: d_coarse_sol_mass(nlayers)
real(kind=real_umphys), intent(in) :: d_coarse_sol_number(nlayers)

! Change in the dust mass [kg kg-1] and number [kg-1] of each mode
real(kind=real_umphys), intent(in) :: d_accum_dust_mass(nlayers)
real(kind=real_umphys), intent(in) :: d_accum_dust_number(nlayers)
real(kind=real_umphys), intent(in) :: d_coarse_dust_mass(nlayers)
real(kind=real_umphys), intent(in) :: d_coarse_dust_number(nlayers)

! Prognostic soluble mass [kg kg-1] and number [kg-1] of each mode
real(kind=real_umphys), intent(inout) :: aitken_sol_mass(nlayers)
real(kind=real_umphys), intent(inout) :: aitken_sol_number(nlayers)
real(kind=real_umphys), intent(inout) :: accum_sol_mass(nlayers)
real(kind=real_umphys), intent(inout) :: accum_sol_number(nlayers)
real(kind=real_umphys), intent(inout) :: coarse_sol_mass(nlayers)
real(kind=real_umphys), intent(inout) :: coarse_sol_number(nlayers)

! Prognostic dust mass [kg kg-1] and number [kg-1] of each mode
real(kind=real_umphys), intent(inout) :: accum_dust_mass(nlayers)
real(kind=real_umphys), intent(inout) :: accum_dust_number(nlayers)
real(kind=real_umphys), intent(inout) :: coarse_dust_mass(nlayers)
real(kind=real_umphys), intent(inout) :: coarse_dust_number(nlayers)

!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------

integer :: k ! Loop counter over levels

real(kind=real_umphys), parameter :: zero = 0.0_real_umphys

!-------------------------------------------------------------------------------
! End of declarations and start of subroutine
!-------------------------------------------------------------------------------

if ( l_tracer_aerosol ) then

  do k = 1, nlayers
    aitken_sol_mass(k)    = max( zero,                                         &
                                 aitken_sol_mass(k)  + d_aitken_sol_mass(k) )
    aitken_sol_number(k)  = max( zero,                                         &
                                 aitken_sol_number(k) + d_aitken_sol_number(k) )
    accum_sol_mass(k)     = max( zero,                                         &
                                 accum_sol_mass(k)   + d_accum_sol_mass(k) )
    accum_sol_number(k)   = max( zero,                                         &
                                 accum_sol_number(k) + d_accum_sol_number(k) )
    coarse_sol_mass(k)    = max( zero,                                         &
                                 coarse_sol_mass(k)  + d_coarse_sol_mass(k) )
    coarse_sol_number(k)  = max( zero,                                         &
                                 coarse_sol_number(k) + d_coarse_sol_number(k) )
    accum_dust_mass(k)    = max( zero,                                         &
                                 accum_dust_mass(k)  + d_accum_dust_mass(k) )
    accum_dust_number(k)  = max( zero,                                         &
                                 accum_dust_number(k) + d_accum_dust_number(k) )
    coarse_dust_mass(k)   = max( zero,                                         &
                                 coarse_dust_mass(k) + d_coarse_dust_mass(k) )
    coarse_dust_number(k) = max( zero,                                         &
                               coarse_dust_number(k) + d_coarse_dust_number(k) )
  end do

end if ! l_tracer_aerosol

end subroutine aerosol_convert_return_ft

end module aerosol_convert_return_mod
