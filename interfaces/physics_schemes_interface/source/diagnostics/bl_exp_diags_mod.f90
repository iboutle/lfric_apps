!-------------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Processes diagnostics for bl_exp_alg

module bl_exp_diags_mod

  use constants_mod,       only: l_def
  use field_mod,           only: field_type
  use integer_field_mod,   only: integer_field_type
  use timing_mod,          only: start_timing, stop_timing, tik, LPROF
  use initialise_diagnostics_mod,     only : init_diag => init_diagnostic_field

  implicit none

  private

  ! Logical indicating whether diagnostics are requested
  logical( l_def ) :: zht_flag, oblen_flag
  logical( l_def ) :: weparm_flag, weparm_dsc_flag, dbtop_flag,               &
                      dbtop_dsc_flag, dftop_flag, dftop_dsc_flag

  public :: initialise_diags_for_bl_exp
  public :: output_diags_for_bl_exp

contains

  !> @brief Initialise fields for locally-computed diagnostics
  !> @param[inout] zht                Turbulent mixing height
  !> @param[inout] oblen              Obukhov length
  !> @param[inout] weparm             SML-top entrainment rate
  !> @param[inout] weparm_dsc         BL-top entrainment rate
  !> @param[inout] dbtop              Buoyancy jump across SML inversion
  !> @param[inout] dbtop_dsc          Buoyancy jump across BL-top inversion
  !> @param[inout] dftop              Radiative flux jump across SML top / cp
  !> @param[inout] dftop_dsc          Radiative flux jump across BL top / cp
  subroutine initialise_diags_for_bl_exp(zht, oblen, weparm, weparm_dsc,       &
                                         dbtop, dbtop_dsc, dftop, dftop_dsc)

    implicit none

    type( field_type ), intent(inout) :: zht
    type( field_type ), intent(inout) :: oblen
    type( field_type ), intent(inout) :: weparm, weparm_dsc
    type( field_type ), intent(inout) :: dbtop, dbtop_dsc
    type( field_type ), intent(inout) :: dftop, dftop_dsc
    integer( tik )                    :: id

    if ( LPROF ) call start_timing( id, 'diags.bl_exp' )

    zht_flag = init_diag(zht, 'turbulence__zht')
    oblen_flag = init_diag(oblen, 'turbulence__oblen')
    weparm_flag = init_diag(weparm, 'turbulence__weparm')
    weparm_dsc_flag = init_diag(weparm_dsc, 'turbulence__weparm_dsc')
    dbtop_flag = init_diag(dbtop, 'turbulence__dbtop')
    dbtop_dsc_flag = init_diag(dbtop_dsc, 'turbulence__dbtop_dsc')
    dftop_flag = init_diag(dftop, 'turbulence__dftop')
    dftop_dsc_flag = init_diag(dftop_dsc, 'turbulence__dftop_dsc')

    if ( LPROF ) call stop_timing( id, 'diags.bl_exp' )

  end subroutine initialise_diags_for_bl_exp

  !> @brief Output diagnostics from bl_exp_alg
  !> @param[in] ntml              Number of turbulently mixed levels
  !> @param[in] cumulus           Cumulus flag (true/false)
  !> @param[in] bl_type_ind       Diagnosed BL types
  !> @param[in] wvar              Vertical velocity variance in wth
  !> @param[in] dsldzm            Liquid potential temperature gradient in wth
  !> @param[in] mix_len_bm        Turb length-scale for bimodal in wth
  !> @param[in] gradrinr          Gradient Richardson number in wth
  !> @param[in] rhokh_bl          Heat eddy diffusivity on BL levels
  !> @param[in] tke_bl            Turbulent kinetic energy (m2 s-2)
  !> @param[in] dtrdz_tq_bl       dt/(rho*r*r*dz) in wth
  !> @param[in] rdz_tq_bl         1/dz in w3
  !> @param[in] zhsc              Height of decoupled layer top
  !> @param[in] level_ent         Level of surface mixed layer inversion
  !> @param[in] level_ent_dsc     Level of decoupled stratocumulus inversion
  !> @param[in] ent_we_lim        Rho * entrainment rate at surface ML inversion (kg m-2 s-1)
  !> @param[in] ent_t_frac        Fraction of time surface ML inversion is above level
  !> @param[in] ent_zrzi          Level height as fraction of DSC inversion height above DSC ML base
  !> @param[in] ent_we_lim_dsc    Rho * entrainment rate at DSC inversion (kg m-2 s-1)
  !> @param[in] ent_t_frac_dsc    Fraction of time DSC inversion is above level
  !> @param[in] ent_zrzi_dsc      Level height as fraction of DSC inversion height above DSC ML base
  !> @param[in] zht               Turbulent mixing height
  !> @param[in] oblen             Obukhov length
  !> @param[in] bl_weight_1dbl    Blending weight to 1D BL scheme in the BL
  !> @param[in] inv_depth         Inversion thickness
  !> @param[in] weparm            SML-top entrainment rate
  !> @param[in] weparm_dsc        BL-top entrainment rate
  !> @param[in] dbtop             Buoyancy jump across SML inversion
  !> @param[in] dbtop_dsc         Buoyancy jump across BL-top inversion
  !> @param[in] dftop             Radiative flux jump across SML top / cp
  !> @param[in] dftop_dsc         Radiative flux jump across BL top / cp
  subroutine output_diags_for_bl_exp(ntml, cumulus, bl_type_ind,               &
                                     wvar, dsldzm, mix_len_bm,                 &
                                     gradrinr, rhokh_bl, tke_bl, dtrdz_tq_bl,  &
                                     rdz_tq_bl, zhsc, level_ent, level_ent_dsc,&
                                     ent_we_lim, ent_t_frac, ent_zrzi,         &
                                     ent_we_lim_dsc, ent_t_frac_dsc,           &
                                     ent_zrzi_dsc, zht, oblen,                 &
                                     bl_weight_1dbl, inv_depth,                &
                                     weparm, weparm_dsc, dbtop, dbtop_dsc,     &
                                     dftop, dftop_dsc)

    implicit none

    ! Prognostic fields to output
    type( field_type ), intent(in)    :: wvar, dsldzm, mix_len_bm, gradrinr,   &
                                         rhokh_bl, tke_bl,                     &
                                         dtrdz_tq_bl, rdz_tq_bl, zhsc,         &
                                         zht,                                  &
                                         ent_we_lim, ent_t_frac, ent_zrzi,     &
                                         ent_we_lim_dsc, ent_t_frac_dsc,       &
                                         ent_zrzi_dsc, oblen, bl_weight_1dbl,  &
                                         inv_depth, weparm, weparm_dsc,        &
                                         dbtop, dbtop_dsc, dftop, dftop_dsc
    type(integer_field_type), intent(in) :: level_ent, level_ent_dsc, ntml,    &
                                            cumulus, bl_type_ind
    integer( tik )  :: id

    if ( LPROF ) call start_timing( id, 'diags.bl_exp' )

    ! Prognostic fields from turbulence collection
    call ntml%write_field('turbulence__ntml')
    call cumulus%write_field('turbulence__cumulus')
    call bl_type_ind%write_field('turbulence__bl_type_ind')
    call wvar%write_field('turbulence__wvar')
    call dsldzm%write_field('turbulence__dsldzm')
    call mix_len_bm%write_field('turbulence__mix_len_bm')
    call gradrinr%write_field('turbulence__gradrinr')
    call rhokh_bl%write_field('turbulence__rhokh')
    call tke_bl%write_field('turbulence__tke')
    call dtrdz_tq_bl%write_field('turbulence__dtrdz_tq')
    call rdz_tq_bl%write_field('turbulence__rdz_tq')
    call zhsc%write_field('turbulence__zhsc')
    call level_ent%write_field('turbulence__level_ent')
    call level_ent_dsc%write_field('turbulence__level_ent_dsc')
    call ent_we_lim%write_field('turbulence__ent_we_lim')
    call ent_t_frac%write_field('turbulence__ent_t_frac')
    call ent_zrzi%write_field('turbulence__ent_zrzi')
    call ent_we_lim_dsc%write_field('turbulence__ent_we_lim_dsc')
    call ent_t_frac_dsc%write_field('turbulence__ent_t_frac_dsc')
    call ent_zrzi_dsc%write_field('turbulence__ent_zrzi_dsc')
    call bl_weight_1dbl%write_field('turbulence__bl_weight_1dbl')
    call inv_depth%write_field('turbulence__inv_depth')

    ! Diagnostics computed in the kernel
    if (zht_flag) call zht%write_field()
    if (oblen_flag) call oblen%write_field()
    if (weparm_flag) call weparm%write_field()
    if (weparm_dsc_flag) call weparm_dsc%write_field()
    if (dbtop_flag) call dbtop%write_field()
    if (dbtop_dsc_flag) call dbtop_dsc%write_field()
    if (dftop_flag) call dftop%write_field()
    if (dftop_dsc_flag) call dftop_dsc%write_field()

    if ( LPROF ) call stop_timing( id, 'diags.bl_exp' )

  end subroutine output_diags_for_bl_exp
end module bl_exp_diags_mod
