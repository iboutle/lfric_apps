!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!-------------------------------------------------------------------------------
!> @brief Interface to the CASIM mechanistic cloud droplet activation.

module casim_aerosol_act_kernel_mod

use argument_mod,      only: arg_type,                  &
                             GH_FIELD, GH_REAL,         &
                             GH_READ, GH_READWRITE,     &
                             CELL_COLUMN
use fs_continuity_mod, only: WTHETA
use kernel_mod,        only: kernel_type
use aerosol_config_mod, only: murk_prognostic

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: casim_aerosol_act_kernel_type
  private
  type(arg_type) :: meta_args(51) = (/                                      &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! m_cl
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! m_cl_pre_fast
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! m_r
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! cf_liq
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! cf_liq_pre_fast
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! theta_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! exner_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! rho_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! w_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! wvar
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! nl_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! nr_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! murk
       ! Prognostic CASIM aerosol
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! aitken_sol_mass_in
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! aitken_sol_number_in
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! accum_sol_mass_in
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! accum_sol_number_in
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! coarse_sol_mass_in
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! coarse_sol_number_in
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! accum_dust_mass_in
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! accum_dust_number_in
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! coarse_dust_mass_in
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! coarse_dust_number_in
       ! CASIM aerosol processing prognostics
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! active_sol_liquid
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! active_sol_rain
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! active_insol_ice
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! active_sol_ice
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! active_insol_liquid
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! active_sol_number
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! active_insol_number
       ! GLOMAP modal aerosol
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! n_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! ait_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! ait_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! ait_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! n_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! acc_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! acc_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! acc_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! acc_sol_ss
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! n_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! cor_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! cor_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! cor_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! cor_sol_ss
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! n_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! ait_ins_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! ait_ins_om
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! n_acc_ins
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! acc_ins_du
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! n_cor_ins
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA)                    & ! cor_ins_du
       /)
   integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: casim_aerosol_act_code
end type

public :: casim_aerosol_act_code

contains

!> @brief Mechanistic activation of cloud droplets on the aerosol
!> @details Reproduces the aerosol activation branch of the UM routine
!>          casim_step_cloud. Where liquid cloud has been created by physical
!>          processes outside of CASIM, the CASIM activation scheme is called to
!>          work out how many cloud droplets that cloud should hold, given the
!>          aerosol available. Where liquid cloud has been removed, the droplets
!>          are removed with it. For aerosol processing runs, the activated
!>          aerosol is moved out of the interstitial aerosol and into the
!>          activated aerosol prognostics. Where the interstitial aerosol came
!>          from the CASIM tracer prognostics it is depleted in place; where it
!>          came from GLOMAP and two way coupling has been asked for, the
!>          depletion is converted back into the GLOMAP units and applied to
!>          the GLOMAP modes.
!> @param[in]     nlayers               Number of layers
!> @param[in]     m_cl                  Cloud liquid mixing ratio in wth
!> @param[in]     m_cl_pre_fast         Cloud liquid mixing ratio as it was
!>                                      before the fast physics was called
!> @param[in]     m_r                   Rain mixing ratio in wth
!> @param[in]     cf_liq                Liquid cloud fraction
!> @param[in]     cf_liq_pre_fast       Liquid cloud fraction as it was
!>                                      before the fast physics was called
!> @param[in]     theta_in_wth          Potential temperature field
!> @param[in]     exner_in_wth          Exner pressure in wth
!> @param[in]     rho_in_wth            Dry density in wth
!> @param[in]     w_in_wth              'Vertical' wind in wth
!> @param[in]     wvar                  Vertical velocity variance from the BL
!> @param[in,out] nl_mphys              CASIM cloud-droplet number concentration
!> @param[in]     nr_mphys              CASIM rain-drop number concentration
!> @param[in]     murk                  Murk aerosol
!> @param[in,out] aitken_sol_mass_in    Prognostic soluble Aitken mode mass
!> @param[in,out] aitken_sol_number_in  Prognostic soluble Aitken number
!> @param[in,out] accum_sol_mass_in     Prognostic soluble accum mode mass
!> @param[in,out] accum_sol_number_in   Prognostic soluble accum mode number
!> @param[in,out] coarse_sol_mass_in    Prognostic soluble coarse mode mass
!> @param[in,out] coarse_sol_number_in  Prognostic soluble coarse mode number
!> @param[in]     accum_dust_mass_in    Prognostic accum mode dust mass
!> @param[in]     accum_dust_number_in  Prognostic accum mode dust number
!> @param[in,out] coarse_dust_mass_in   Prognostic coarse mode dust mass
!> @param[in,out] coarse_dust_number_in Prognostic coarse mode dust number
!> @param[in,out] active_sol_liquid     Soluble aerosol mass activated in liquid
!> @param[in]     active_sol_rain       Soluble aerosol mass activated in rain
!> @param[in]     active_insol_ice      Insoluble aerosol mass activated in ice
!> @param[in]     active_sol_ice        Soluble aerosol mass activated in ice
!> @param[in,out] active_insol_liquid   Insoluble aerosol mass in liquid
!> @param[in,out] active_sol_number     Soluble aerosol number in cloud
!> @param[in,out] active_insol_number   Insoluble aerosol number in cloud
!> @param[in,out] n_ait_sol             GLOMAP Aitken soluble number mr
!> @param[in,out] ait_sol_su            GLOMAP Aitken soluble sulphate mmr
!> @param[in]     ait_sol_bc            GLOMAP Aitken soluble black carbon mmr
!> @param[in]     ait_sol_om            GLOMAP Aitken soluble organic carbon mmr
!> @param[in,out] n_acc_sol             GLOMAP accum soluble number mr
!> @param[in,out] acc_sol_su            GLOMAP accum soluble sulphate mmr
!> @param[in]     acc_sol_bc            GLOMAP accum soluble black carbon mmr
!> @param[in]     acc_sol_om            GLOMAP accum soluble organic carbon mmr
!> @param[in]     acc_sol_ss            GLOMAP accum soluble sea salt mmr
!> @param[in,out] n_cor_sol             GLOMAP coarse soluble number mr
!> @param[in,out] cor_sol_su            GLOMAP coarse soluble sulphate mmr
!> @param[in]     cor_sol_bc            GLOMAP coarse soluble black carbon mmr
!> @param[in]     cor_sol_om            GLOMAP coarse soluble organic carbon mmr
!> @param[in]     cor_sol_ss            GLOMAP coarse soluble sea salt mmr
!> @param[in]     n_ait_ins             GLOMAP Aitken insoluble number mr
!> @param[in]     ait_ins_bc            GLOMAP Aitken insoluble black carbon mmr
!> @param[in]     ait_ins_om            GLOMAP Aitken insoluble organic carbon mmr
!> @param[in,out] n_acc_ins             GLOMAP accum insoluble number mr
!> @param[in,out] acc_ins_du            GLOMAP accum insoluble dust mmr
!> @param[in,out] n_cor_ins             GLOMAP coarse insoluble number mr
!> @param[in,out] cor_ins_du            GLOMAP coarse insoluble dust mmr
!> @param[in]     ndf_wth               Number of degrees of freedom per cell for
!!                                       potential temperature space
!> @param[in]     undf_wth              Number unique of degrees of freedom for
!!                                       potential temperature space
!> @param[in]     map_wth               Dofmap for the cell at the base of the
!!                                       column for potential temperature space

subroutine casim_aerosol_act_code( nlayers,                     &
                                   m_cl, m_cl_pre_fast, m_r,    &
                                   cf_liq, cf_liq_pre_fast,     &
                                   theta_in_wth, exner_in_wth,  &
                                   rho_in_wth, w_in_wth, wvar,  &
                                   nl_mphys, nr_mphys, murk,    &
                                   aitken_sol_mass_in,          &
                                   aitken_sol_number_in,        &
                                   accum_sol_mass_in,           &
                                   accum_sol_number_in,         &
                                   coarse_sol_mass_in,          &
                                   coarse_sol_number_in,        &
                                   accum_dust_mass_in,          &
                                   accum_dust_number_in,        &
                                   coarse_dust_mass_in,         &
                                   coarse_dust_number_in,       &
                                   active_sol_liquid,           &
                                   active_sol_rain,             &
                                   active_insol_ice,            &
                                   active_sol_ice,              &
                                   active_insol_liquid,         &
                                   active_sol_number,           &
                                   active_insol_number,         &
                                   n_ait_sol, ait_sol_su,       &
                                   ait_sol_bc, ait_sol_om,      &
                                   n_acc_sol, acc_sol_su,       &
                                   acc_sol_bc, acc_sol_om,      &
                                   acc_sol_ss,                  &
                                   n_cor_sol, cor_sol_su,       &
                                   cor_sol_bc, cor_sol_om,      &
                                   cor_sol_ss,                  &
                                   n_ait_ins, ait_ins_bc,       &
                                   ait_ins_om,                  &
                                   n_acc_ins, acc_ins_du,       &
                                   n_cor_ins, cor_ins_du,       &
                                   ndf_wth, undf_wth, map_wth   )

    use constants_mod,              only: r_def, i_def, r_um, i_um

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use planet_constants_mod,       only: p_zero, kappa
    use casim_switches,             only: l_fix_aerosol,                &
                                          l_tracer_aerosol,             &
                                          l_ukca_aerosol,               &
                                          l_ukca_feeding_out,           &
                                          no_aerosol_modes,             &
                                          no_processing,                &
                                          l_mp_activesolliquid,         &
                                          l_mp_activesolrain,           &
                                          l_mp_activeinsolice,          &
                                          l_mp_activesolice,            &
                                          l_mp_activeinsolliquid,       &
                                          l_mp_activesolnumber,         &
                                          l_mp_activeinsolnumber
    use mphys_inputs_mod,           only: casim_aerosol_option,         &
                                          casim_aerosol_process_level,  &
                                          wvarfac
    use casim_set_dependent_switches_mod,                               &
                                    only: l_process
    use aerosol_extract_convert_mod,                                    &
                                    only: aerosol_extract_convert,      &
                                          aerosol_extract_convert_murk, &
                                          aerosol_extract_convert_ft
    use aerosol_convert_return_mod, only: aerosol_convert_return
    use casim_activation_in_um_mod, only: examine_aerosol_column,       &
                                          examine_processing_column,    &
                                          activate_column,              &
                                          sulphate_mode_bk
    use ukca_mode_setup,            only: nmodes,                       &
                                          mode_ait_sol, mode_acc_sol,   &
                                          mode_cor_sol, mode_ait_insol, &
                                          mode_acc_insol,               &
                                          mode_cor_insol,               &
                                          cp_su, cp_bc, cp_oc, cp_cl,   &
                                          cp_du

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth

    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth

    real(kind=r_def), intent(in),    dimension(undf_wth) :: m_cl
    real(kind=r_def), intent(in),    dimension(undf_wth) :: m_cl_pre_fast
    real(kind=r_def), intent(in),    dimension(undf_wth) :: m_r
    real(kind=r_def), intent(in),    dimension(undf_wth) :: cf_liq
    real(kind=r_def), intent(in),    dimension(undf_wth) :: cf_liq_pre_fast
    real(kind=r_def), intent(in),    dimension(undf_wth) :: theta_in_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: exner_in_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: rho_in_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: w_in_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: wvar
    real(kind=r_def), intent(inout), dimension(undf_wth) :: nl_mphys
    real(kind=r_def), intent(in),    dimension(undf_wth) :: nr_mphys
    real(kind=r_def), intent(in),    dimension(undf_wth) :: murk

    real(kind=r_def), intent(inout), dimension(undf_wth) :: aitken_sol_mass_in
    real(kind=r_def), intent(inout), dimension(undf_wth) :: aitken_sol_number_in
    real(kind=r_def), intent(inout), dimension(undf_wth) :: accum_sol_mass_in
    real(kind=r_def), intent(inout), dimension(undf_wth) :: accum_sol_number_in
    real(kind=r_def), intent(inout), dimension(undf_wth) :: coarse_sol_mass_in
    real(kind=r_def), intent(inout), dimension(undf_wth) :: coarse_sol_number_in
    real(kind=r_def), intent(in),    dimension(undf_wth) :: accum_dust_mass_in
    real(kind=r_def), intent(in),    dimension(undf_wth) :: accum_dust_number_in
    real(kind=r_def), intent(inout), dimension(undf_wth) :: coarse_dust_mass_in
    real(kind=r_def), intent(inout), dimension(undf_wth) :: coarse_dust_number_in

    real(kind=r_def), intent(inout), dimension(undf_wth) :: active_sol_liquid
    real(kind=r_def), intent(in),    dimension(undf_wth) :: active_sol_rain
    real(kind=r_def), intent(in),    dimension(undf_wth) :: active_insol_ice
    real(kind=r_def), intent(in),    dimension(undf_wth) :: active_sol_ice
    real(kind=r_def), intent(inout), dimension(undf_wth) :: active_insol_liquid
    real(kind=r_def), intent(inout), dimension(undf_wth) :: active_sol_number
    real(kind=r_def), intent(inout), dimension(undf_wth) :: active_insol_number

    real(kind=r_def), intent(inout), dimension(undf_wth) :: n_ait_sol
    real(kind=r_def), intent(inout), dimension(undf_wth) :: ait_sol_su
    real(kind=r_def), intent(in),    dimension(undf_wth) :: ait_sol_bc
    real(kind=r_def), intent(in),    dimension(undf_wth) :: ait_sol_om
    real(kind=r_def), intent(inout), dimension(undf_wth) :: n_acc_sol
    real(kind=r_def), intent(inout), dimension(undf_wth) :: acc_sol_su
    real(kind=r_def), intent(in),    dimension(undf_wth) :: acc_sol_bc
    real(kind=r_def), intent(in),    dimension(undf_wth) :: acc_sol_om
    real(kind=r_def), intent(in),    dimension(undf_wth) :: acc_sol_ss
    real(kind=r_def), intent(inout), dimension(undf_wth) :: n_cor_sol
    real(kind=r_def), intent(inout), dimension(undf_wth) :: cor_sol_su
    real(kind=r_def), intent(in),    dimension(undf_wth) :: cor_sol_bc
    real(kind=r_def), intent(in),    dimension(undf_wth) :: cor_sol_om
    real(kind=r_def), intent(in),    dimension(undf_wth) :: cor_sol_ss
    real(kind=r_def), intent(in),    dimension(undf_wth) :: n_ait_ins
    real(kind=r_def), intent(in),    dimension(undf_wth) :: ait_ins_bc
    real(kind=r_def), intent(in),    dimension(undf_wth) :: ait_ins_om
    real(kind=r_def), intent(inout), dimension(undf_wth) :: n_acc_ins
    real(kind=r_def), intent(inout), dimension(undf_wth) :: acc_ins_du
    real(kind=r_def), intent(inout), dimension(undf_wth) :: n_cor_ins
    real(kind=r_def), intent(inout), dimension(undf_wth) :: cor_ins_du

    !-------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------

    ! Column working arrays, on the UM physics precision
    real(r_um), dimension(nlayers) ::                                          &
         t_col, p_col, rho_col, w_tke, murk_col,                               &
         qcl_now_col, qcl_pre_ap2_col, cf_liquid_col, cf_liq_pre_ap2_col,      &
         cloudnumber_col, rainnumber_col, rainmass_col,                        &
         ait_sol_mass_um,  ait_sol_num_um,                                     &
         acc_sol_mass_um,  acc_sol_num_um,                                     &
         cor_sol_mass_um,  cor_sol_num_um,                                     &
         acc_dust_mass_um, acc_dust_num_um,                                    &
         cor_dust_mass_um, cor_dust_num_um,                                    &
         ait_sol_bk_um,    acc_sol_bk_um,   cor_sol_bk_um,                     &
         prog_ait_sol_mass,  prog_ait_sol_num,                                 &
         prog_acc_sol_mass,  prog_acc_sol_num,                                 &
         prog_cor_sol_mass,  prog_cor_sol_num,                                 &
         prog_acc_dust_mass, prog_acc_dust_num,                                &
         prog_cor_dust_mass, prog_cor_dust_num,                                &
         act_sol_liq_col,   act_sol_rain_col,  act_sol_ice_col,                &
         act_sol_num_col,   act_insol_liq_col, act_insol_ice_col,              &
         act_insol_num_col

    ! Change made to the aerosol by the activation, used to return the
    ! depleted aerosol to GLOMAP
    real(r_um), dimension(nlayers) ::                                          &
         d_ait_sol_mass,  d_ait_sol_num,                                       &
         d_acc_sol_mass,  d_acc_sol_num,                                       &
         d_cor_sol_mass,  d_cor_sol_num,                                       &
         d_acc_dust_mass, d_acc_dust_num,                                      &
         d_cor_dust_mass, d_cor_dust_num

    ! At a later date we would like to obtain this as a variable via the api.
    ! This value is hard coded for the time being to work with
    ! ukca_mode_sussbcocdu_7mode, matching glomap_aerosol_kernel_mod.
    integer(i_um), parameter :: ncp_lfric = 6

    ! GLOMAP mass and number mixing ratios gathered by mode and component
    real(r_um), dimension(nlayers,nmodes,ncp_lfric) :: mode_mmr_um
    real(r_um), dimension(nlayers,nmodes)           :: mode_nmr_um

    ! Whether the hygroscopicity supplied by the extraction is to be used
    logical :: l_set_bk

    ! Whether the aerosol supplied to CASIM came from murk, which caps the
    ! droplet number
    logical :: l_murk_source

    ! Whether the aerosol supplied to CASIM came from GLOMAP
    logical :: l_ukca_source

    ! Whether the aerosol supplied to CASIM is prognostic, and so may be
    ! depleted by the activation
    logical :: l_process_aerosol

    ! Whether the aerosol processing prognostics are in use
    logical :: l_processing_on

    integer(i_um) :: k

    ! Bounds applied to the vertical velocity used for the activation, as in
    ! the UM routine casim_step_cloud
    real(r_um), parameter :: min_velocity = 0.01_r_um
    real(r_um), parameter :: max_velocity = 4.00_r_um

    real(r_um), parameter :: eps = epsilon(1.0_r_um)

    !-------------------------------------------------------------------------
    ! End of Declarations
    !-------------------------------------------------------------------------

    l_processing_on   = ( casim_aerosol_process_level > no_processing )
    l_ukca_source     = ( l_ukca_aerosol .and.                                 &
                          casim_aerosol_option > no_aerosol_modes )
    l_murk_source     = .false.

    ! The aerosol may only be depleted where the model has somewhere to put
    ! the depleted aerosol back, that is where it came from the CASIM tracer
    ! prognostics or where GLOMAP has been set up to take it back.
    l_process_aerosol = l_tracer_aerosol .or.                                  &
                        ( l_ukca_source .and. l_ukca_feeding_out )

    !-----------------------------------------------------------------------
    ! Gather the column
    !-----------------------------------------------------------------------
    do k = 1, nlayers

      t_col(k)   = theta_in_wth(map_wth(1) + k) * exner_in_wth(map_wth(1) + k)
      p_col(k)   = p_zero * ( exner_in_wth(map_wth(1) + k) )**(1.0_r_um/kappa)
      rho_col(k) = rho_in_wth(map_wth(1) + k)

      qcl_now_col(k)         = m_cl(map_wth(1) + k)
      qcl_pre_ap2_col(k)     = m_cl_pre_fast(map_wth(1) + k)
      cf_liquid_col(k)       = cf_liq(map_wth(1) + k)
      cf_liq_pre_ap2_col(k)  = cf_liq_pre_fast(map_wth(1) + k)
      cloudnumber_col(k)     = nl_mphys(map_wth(1) + k)
      rainnumber_col(k)      = nr_mphys(map_wth(1) + k)
      rainmass_col(k)        = m_r(map_wth(1) + k)

      ! Add the sub-grid contribution from the boundary layer to the vertical
      ! velocity seen by the activation, then bound it.
      if ( wvar(map_wth(1) + k) > eps ) then
        w_tke(k) = w_in_wth(map_wth(1) + k) +                                  &
                   wvarfac * sqrt( wvar(map_wth(1) + k) )
      else
        w_tke(k) = w_in_wth(map_wth(1) + k)
      end if
      w_tke(k) = min( max_velocity, max( min_velocity, w_tke(k) ) )

    end do

    !-----------------------------------------------------------------------
    ! Extract and convert the aerosol which is to be supplied to CASIM.
    ! This follows the source selection made in casim_kernel_mod.
    !-----------------------------------------------------------------------
    if ( ( l_fix_aerosol .or. l_tracer_aerosol ) .and.                         &
         casim_aerosol_option > no_aerosol_modes ) then

      if ( l_tracer_aerosol ) then
        do k = 1, nlayers
          prog_ait_sol_mass(k)  = aitken_sol_mass_in(map_wth(1) + k)
          prog_ait_sol_num(k)   = aitken_sol_number_in(map_wth(1) + k)
          prog_acc_sol_mass(k)  = accum_sol_mass_in(map_wth(1) + k)
          prog_acc_sol_num(k)   = accum_sol_number_in(map_wth(1) + k)
          prog_cor_sol_mass(k)  = coarse_sol_mass_in(map_wth(1) + k)
          prog_cor_sol_num(k)   = coarse_sol_number_in(map_wth(1) + k)
          prog_acc_dust_mass(k) = accum_dust_mass_in(map_wth(1) + k)
          prog_acc_dust_num(k)  = accum_dust_number_in(map_wth(1) + k)
          prog_cor_dust_mass(k) = coarse_dust_mass_in(map_wth(1) + k)
          prog_cor_dust_num(k)  = coarse_dust_number_in(map_wth(1) + k)
        end do
      else
        prog_ait_sol_mass(:)  = 0.0_r_um
        prog_ait_sol_num(:)   = 0.0_r_um
        prog_acc_sol_mass(:)  = 0.0_r_um
        prog_acc_sol_num(:)   = 0.0_r_um
        prog_cor_sol_mass(:)  = 0.0_r_um
        prog_cor_sol_num(:)   = 0.0_r_um
        prog_acc_dust_mass(:) = 0.0_r_um
        prog_acc_dust_num(:)  = 0.0_r_um
        prog_cor_dust_mass(:) = 0.0_r_um
        prog_cor_dust_num(:)  = 0.0_r_um
      end if

      call aerosol_extract_convert_ft( int(nlayers),                           &
               prog_ait_sol_mass, prog_ait_sol_num,                            &
               prog_acc_sol_mass, prog_acc_sol_num,                            &
               prog_cor_sol_mass, prog_cor_sol_num,                            &
               prog_acc_dust_mass, prog_acc_dust_num,                          &
               prog_cor_dust_mass, prog_cor_dust_num,                          &
               ait_sol_mass_um, ait_sol_num_um,                                &
               acc_sol_mass_um, acc_sol_num_um,                                &
               cor_sol_mass_um, cor_sol_num_um,                                &
               acc_dust_mass_um, acc_dust_num_um,                              &
               cor_dust_mass_um, cor_dust_num_um )

      ! The prognostic aerosol carries no hygroscopicity, so the values
      ! already held by CASIM are left alone.
      ait_sol_bk_um(:) = 0.0_r_um
      acc_sol_bk_um(:) = 0.0_r_um
      cor_sol_bk_um(:) = 0.0_r_um
      l_set_bk = .false.

    else if ( l_ukca_source ) then

      ! Extract the GLOMAP aerosol. The UM reads the modes straight out of
      ! its single UKCA tracer array; LFRic holds each component in its own
      ! field so they are gathered by mode and component here.
      mode_mmr_um(:,:,:) = 0.0_r_um
      mode_nmr_um(:,:)   = 0.0_r_um

      do k = 1, nlayers
        mode_nmr_um(k,mode_ait_sol)       = n_ait_sol(map_wth(1) + k)
        mode_mmr_um(k,mode_ait_sol,cp_su) = ait_sol_su(map_wth(1) + k)
        mode_mmr_um(k,mode_ait_sol,cp_bc) = ait_sol_bc(map_wth(1) + k)
        mode_mmr_um(k,mode_ait_sol,cp_oc) = ait_sol_om(map_wth(1) + k)

        mode_nmr_um(k,mode_acc_sol)       = n_acc_sol(map_wth(1) + k)
        mode_mmr_um(k,mode_acc_sol,cp_su) = acc_sol_su(map_wth(1) + k)
        mode_mmr_um(k,mode_acc_sol,cp_bc) = acc_sol_bc(map_wth(1) + k)
        mode_mmr_um(k,mode_acc_sol,cp_oc) = acc_sol_om(map_wth(1) + k)
        mode_mmr_um(k,mode_acc_sol,cp_cl) = acc_sol_ss(map_wth(1) + k)

        mode_nmr_um(k,mode_cor_sol)       = n_cor_sol(map_wth(1) + k)
        mode_mmr_um(k,mode_cor_sol,cp_su) = cor_sol_su(map_wth(1) + k)
        mode_mmr_um(k,mode_cor_sol,cp_bc) = cor_sol_bc(map_wth(1) + k)
        mode_mmr_um(k,mode_cor_sol,cp_oc) = cor_sol_om(map_wth(1) + k)
        mode_mmr_um(k,mode_cor_sol,cp_cl) = cor_sol_ss(map_wth(1) + k)

        mode_nmr_um(k,mode_ait_insol)       = n_ait_ins(map_wth(1) + k)
        mode_mmr_um(k,mode_ait_insol,cp_bc) = ait_ins_bc(map_wth(1) + k)
        mode_mmr_um(k,mode_ait_insol,cp_oc) = ait_ins_om(map_wth(1) + k)

        mode_nmr_um(k,mode_acc_insol)       = n_acc_ins(map_wth(1) + k)
        mode_mmr_um(k,mode_acc_insol,cp_du) = acc_ins_du(map_wth(1) + k)

        mode_nmr_um(k,mode_cor_insol)       = n_cor_ins(map_wth(1) + k)
        mode_mmr_um(k,mode_cor_insol,cp_du) = cor_ins_du(map_wth(1) + k)
      end do

      ! LFRic has no prognostic dust in the soluble accumulation and coarse
      ! modes, so those components are left at zero.

      call aerosol_extract_convert( int(nlayers), int(ncp_lfric),              &
               p_col, t_col, rho_col,                                          &
               mode_mmr_um, mode_nmr_um,                                       &
               ait_sol_mass_um, ait_sol_num_um,                                &
               acc_sol_mass_um, acc_sol_num_um,                                &
               cor_sol_mass_um, cor_sol_num_um,                                &
               acc_dust_mass_um, acc_dust_num_um,                              &
               cor_dust_mass_um, cor_dust_num_um,                              &
               ait_sol_bk_um, acc_sol_bk_um, cor_sol_bk_um )

      l_set_bk = .true.

    else if ( murk_prognostic ) then

      ! Convert murk into soluble accumulation mode mass and number
      do k = 1, nlayers
        murk_col(k) = murk(map_wth(1) + k)
      end do

      call aerosol_extract_convert_murk( int(nlayers), rho_col,                &
               murk_col,                                                       &
               ait_sol_mass_um, ait_sol_num_um,                                &
               acc_sol_mass_um, acc_sol_num_um,                                &
               cor_sol_mass_um, cor_sol_num_um,                                &
               acc_dust_mass_um, acc_dust_num_um,                              &
               cor_dust_mass_um, cor_dust_num_um )

      ! The murk aerosol is assumed to be ammonium sulphate, following the UM
      ! routine examine_murk_aerosol_column.
      ait_sol_bk_um(:) = 0.0_r_um
      acc_sol_bk_um(:) = sulphate_mode_bk
      cor_sol_bk_um(:) = 0.0_r_um
      l_set_bk      = .true.
      l_murk_source = .true.

    else

      ! No aerosol supplied to CASIM (failsafe option)
      ait_sol_mass_um(:)  = 0.0_r_um
      ait_sol_num_um(:)   = 0.0_r_um
      acc_sol_mass_um(:)  = 0.0_r_um
      acc_sol_num_um(:)   = 0.0_r_um
      cor_sol_mass_um(:)  = 0.0_r_um
      cor_sol_num_um(:)   = 0.0_r_um
      acc_dust_mass_um(:) = 0.0_r_um
      acc_dust_num_um(:)  = 0.0_r_um
      cor_dust_mass_um(:) = 0.0_r_um
      cor_dust_num_um(:)  = 0.0_r_um
      ait_sol_bk_um(:)    = 0.0_r_um
      acc_sol_bk_um(:)    = 0.0_r_um
      cor_sol_bk_um(:)    = 0.0_r_um
      l_set_bk = .false.

    end if ! aerosol source

    !-----------------------------------------------------------------------
    ! Set up the aerosol held inside CASIM
    !-----------------------------------------------------------------------
    call examine_aerosol_column( int(nlayers),                                 &
             ait_sol_mass_um, ait_sol_num_um,                                  &
             acc_sol_mass_um, acc_sol_num_um,                                  &
             cor_sol_mass_um, cor_sol_num_um,                                  &
             acc_dust_mass_um, acc_dust_num_um,                                &
             cor_dust_mass_um, cor_dust_num_um,                                &
             ait_sol_bk_um, acc_sol_bk_um, cor_sol_bk_um,                      &
             l_set_bk )

    !-----------------------------------------------------------------------
    ! For aerosol processing runs, supply the aerosol processing prognostics
    ! to CASIM.
    !-----------------------------------------------------------------------
    act_sol_liq_col(:)   = 0.0_r_um
    act_sol_rain_col(:)  = 0.0_r_um
    act_sol_ice_col(:)   = 0.0_r_um
    act_sol_num_col(:)   = 0.0_r_um
    act_insol_liq_col(:) = 0.0_r_um
    act_insol_ice_col(:) = 0.0_r_um
    act_insol_num_col(:) = 0.0_r_um

    if ( l_processing_on ) then

      if ( l_mp_activesolliquid ) then
        do k = 1, nlayers
          act_sol_liq_col(k) = active_sol_liquid(map_wth(1) + k)
        end do
      end if

      if ( l_mp_activesolrain ) then
        do k = 1, nlayers
          act_sol_rain_col(k) = active_sol_rain(map_wth(1) + k)
        end do
      end if

      if ( l_mp_activesolice ) then
        do k = 1, nlayers
          act_sol_ice_col(k) = active_sol_ice(map_wth(1) + k)
        end do
      end if

      if ( l_mp_activesolnumber ) then
        do k = 1, nlayers
          act_sol_num_col(k) = active_sol_number(map_wth(1) + k)
        end do
      end if

      if ( l_mp_activeinsolliquid ) then
        do k = 1, nlayers
          act_insol_liq_col(k) = active_insol_liquid(map_wth(1) + k)
        end do
      end if

      if ( l_mp_activeinsolice ) then
        do k = 1, nlayers
          act_insol_ice_col(k) = active_insol_ice(map_wth(1) + k)
        end do
      end if

      if ( l_mp_activeinsolnumber ) then
        do k = 1, nlayers
          act_insol_num_col(k) = active_insol_number(map_wth(1) + k)
        end do
      end if

      call examine_processing_column( int(nlayers),                            &
               cloudnumber_col, qcl_now_col,                                   &
               rainnumber_col, rainmass_col,                                   &
               act_sol_liq_col, act_sol_rain_col,                              &
               act_sol_ice_col, act_sol_num_col,                               &
               act_insol_liq_col, act_insol_ice_col,                           &
               act_insol_num_col )

    end if ! l_processing_on

    ! Take a copy of the aerosol so that the change made by the activation can
    ! be worked out afterwards. The accumulation mode dust is included for
    ! completeness even though the activation cannot alter it.
    do k = 1, nlayers
      d_ait_sol_mass(k)  = -ait_sol_mass_um(k)
      d_ait_sol_num(k)   = -ait_sol_num_um(k)
      d_acc_sol_mass(k)  = -acc_sol_mass_um(k)
      d_acc_sol_num(k)   = -acc_sol_num_um(k)
      d_cor_sol_mass(k)  = -cor_sol_mass_um(k)
      d_cor_sol_num(k)   = -cor_sol_num_um(k)
      d_acc_dust_mass(k) = -acc_dust_mass_um(k)
      d_acc_dust_num(k)  = -acc_dust_num_um(k)
      d_cor_dust_mass(k) = -cor_dust_mass_um(k)
      d_cor_dust_num(k)  = -cor_dust_num_um(k)
    end do

    !-----------------------------------------------------------------------
    ! Activate the cloud droplets
    !-----------------------------------------------------------------------
    call activate_column( int(nlayers), l_process_aerosol, l_murk_source,       &
             qcl_now_col, qcl_pre_ap2_col,                                     &
             rho_col, t_col, p_col,                                            &
             cf_liquid_col, cf_liq_pre_ap2_col, w_tke,                         &
             cloudnumber_col, rainnumber_col,                                  &
             ait_sol_mass_um, ait_sol_num_um,                                  &
             acc_sol_mass_um, acc_sol_num_um,                                  &
             cor_sol_mass_um, cor_sol_num_um,                                  &
             cor_dust_mass_um, cor_dust_num_um,                                &
             act_sol_liq_col, act_sol_num_col,                                 &
             act_insol_liq_col, act_insol_ice_col,                             &
             act_insol_num_col )

    do k = 1, nlayers
      d_ait_sol_mass(k)  = d_ait_sol_mass(k)  + ait_sol_mass_um(k)
      d_ait_sol_num(k)   = d_ait_sol_num(k)   + ait_sol_num_um(k)
      d_acc_sol_mass(k)  = d_acc_sol_mass(k)  + acc_sol_mass_um(k)
      d_acc_sol_num(k)   = d_acc_sol_num(k)   + acc_sol_num_um(k)
      d_cor_sol_mass(k)  = d_cor_sol_mass(k)  + cor_sol_mass_um(k)
      d_cor_sol_num(k)   = d_cor_sol_num(k)   + cor_sol_num_um(k)
      d_acc_dust_mass(k) = d_acc_dust_mass(k) + acc_dust_mass_um(k)
      d_acc_dust_num(k)  = d_acc_dust_num(k)  + acc_dust_num_um(k)
      d_cor_dust_mass(k) = d_cor_dust_mass(k) + cor_dust_mass_um(k)
      d_cor_dust_num(k)  = d_cor_dust_num(k)  + cor_dust_num_um(k)
    end do

    !-----------------------------------------------------------------------
    ! Scatter the results back into the fields
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      nl_mphys(map_wth(1) + k) = cloudnumber_col(k)
    end do
    nl_mphys(map_wth(1)) = nl_mphys(map_wth(1) + 1)

    ! The interstitial aerosol is only depleted where it is prognostic.
    if ( l_process .and. l_process_aerosol ) then

      ! Where the aerosol came from the CASIM tracer prognostics the depleted
      ! aerosol goes straight back into them.
      if ( l_tracer_aerosol ) then

        do k = 1, nlayers
          aitken_sol_mass_in(map_wth(1) + k)    = ait_sol_mass_um(k)
          aitken_sol_number_in(map_wth(1) + k)  = ait_sol_num_um(k)
          accum_sol_mass_in(map_wth(1) + k)     = acc_sol_mass_um(k)
          accum_sol_number_in(map_wth(1) + k)   = acc_sol_num_um(k)
          coarse_sol_mass_in(map_wth(1) + k)    = cor_sol_mass_um(k)
          coarse_sol_number_in(map_wth(1) + k)  = cor_sol_num_um(k)
          coarse_dust_mass_in(map_wth(1) + k)   = cor_dust_mass_um(k)
          coarse_dust_number_in(map_wth(1) + k) = cor_dust_num_um(k)
        end do
        aitken_sol_mass_in(map_wth(1))    = aitken_sol_mass_in(map_wth(1) + 1)
        aitken_sol_number_in(map_wth(1))  = aitken_sol_number_in(map_wth(1) + 1)
        accum_sol_mass_in(map_wth(1))     = accum_sol_mass_in(map_wth(1) + 1)
        accum_sol_number_in(map_wth(1))   = accum_sol_number_in(map_wth(1) + 1)
        coarse_sol_mass_in(map_wth(1))    = coarse_sol_mass_in(map_wth(1) + 1)
        coarse_sol_number_in(map_wth(1))  = coarse_sol_number_in(map_wth(1) + 1)
        coarse_dust_mass_in(map_wth(1))   = coarse_dust_mass_in(map_wth(1) + 1)
        coarse_dust_number_in(map_wth(1)) =                                    &
                                        coarse_dust_number_in(map_wth(1) + 1)

      end if ! l_tracer_aerosol

      ! Where the aerosol came from GLOMAP the change is converted back out of
      ! the CASIM units and applied to the GLOMAP modes. The UM does this
      ! inline in activate_column_2way_ukca; LFRic keeps the activation itself
      ! free of the aerosol source and converts afterwards, in the same way as
      ! the UM routine aerosol_convert_return_solinsol.
      if ( l_ukca_source ) then

        call aerosol_convert_return( int(nlayers), int(ncp_lfric),             &
                 p_col, t_col, rho_col,                                        &
                 d_ait_sol_mass,  d_ait_sol_num,                               &
                 d_acc_sol_mass,  d_acc_sol_num,                               &
                 d_cor_sol_mass,  d_cor_sol_num,                               &
                 d_acc_dust_mass, d_acc_dust_num,                              &
                 d_cor_dust_mass, d_cor_dust_num,                              &
                 mode_mmr_um, mode_nmr_um )

        do k = 1, nlayers
          n_ait_sol(map_wth(1) + k)  = mode_nmr_um(k,mode_ait_sol)
          ait_sol_su(map_wth(1) + k) = mode_mmr_um(k,mode_ait_sol,cp_su)

          n_acc_sol(map_wth(1) + k)  = mode_nmr_um(k,mode_acc_sol)
          acc_sol_su(map_wth(1) + k) = mode_mmr_um(k,mode_acc_sol,cp_su)

          n_cor_sol(map_wth(1) + k)  = mode_nmr_um(k,mode_cor_sol)
          cor_sol_su(map_wth(1) + k) = mode_mmr_um(k,mode_cor_sol,cp_su)

          n_acc_ins(map_wth(1) + k)  = mode_nmr_um(k,mode_acc_insol)
          acc_ins_du(map_wth(1) + k) = mode_mmr_um(k,mode_acc_insol,cp_du)

          n_cor_ins(map_wth(1) + k)  = mode_nmr_um(k,mode_cor_insol)
          cor_ins_du(map_wth(1) + k) = mode_mmr_um(k,mode_cor_insol,cp_du)
        end do
        n_ait_sol(map_wth(1))  = n_ait_sol(map_wth(1) + 1)
        ait_sol_su(map_wth(1)) = ait_sol_su(map_wth(1) + 1)
        n_acc_sol(map_wth(1))  = n_acc_sol(map_wth(1) + 1)
        acc_sol_su(map_wth(1)) = acc_sol_su(map_wth(1) + 1)
        n_cor_sol(map_wth(1))  = n_cor_sol(map_wth(1) + 1)
        cor_sol_su(map_wth(1)) = cor_sol_su(map_wth(1) + 1)
        n_acc_ins(map_wth(1))  = n_acc_ins(map_wth(1) + 1)
        acc_ins_du(map_wth(1)) = acc_ins_du(map_wth(1) + 1)
        n_cor_ins(map_wth(1))  = n_cor_ins(map_wth(1) + 1)
        cor_ins_du(map_wth(1)) = cor_ins_du(map_wth(1) + 1)

      end if ! l_ukca_source

      if ( l_mp_activesolliquid ) then
        do k = 1, nlayers
          active_sol_liquid(map_wth(1) + k) = act_sol_liq_col(k)
        end do
        active_sol_liquid(map_wth(1)) = active_sol_liquid(map_wth(1) + 1)
      end if

      if ( l_mp_activesolnumber ) then
        do k = 1, nlayers
          active_sol_number(map_wth(1) + k) = act_sol_num_col(k)
        end do
        active_sol_number(map_wth(1)) = active_sol_number(map_wth(1) + 1)
      end if

      if ( l_mp_activeinsolliquid ) then
        do k = 1, nlayers
          active_insol_liquid(map_wth(1) + k) = act_insol_liq_col(k)
        end do
        active_insol_liquid(map_wth(1)) = active_insol_liquid(map_wth(1) + 1)
      end if

      if ( l_mp_activeinsolnumber ) then
        do k = 1, nlayers
          active_insol_number(map_wth(1) + k) = act_insol_num_col(k)
        end do
        active_insol_number(map_wth(1)) = active_insol_number(map_wth(1) + 1)
      end if

    end if ! l_process .and. l_process_aerosol

end subroutine casim_aerosol_act_code

end module casim_aerosol_act_kernel_mod
