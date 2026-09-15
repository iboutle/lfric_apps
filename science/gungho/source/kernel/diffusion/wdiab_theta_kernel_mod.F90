!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to diabatic mixing of potential temperature
!>
module wdiab_theta_kernel_mod

  use argument_mod,          only : arg_type,                     &
                                    GH_FIELD, GH_REAL,            &
                                    GH_READ, GH_WRITE,            &
                                    CELL_COLUMN, ANY_SPACE_1
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : Wtheta, W2
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : B

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: wdiab_theta_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                       &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta),      & ! dtheta_out
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),      & ! theta
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),      & ! rho
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2),          & ! dA
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), & ! dAsh
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta)       & ! dtheta_in
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: wdiab_theta_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: wdiab_theta_code

contains

!> @details A diabatic heat source of potential temperature will naturally
!>          induce a vertical velocity in response to this. Rather than allow
!>          the solver to produce this vertical velocity, this kernel
!>          pre-mixes the theta increment in response to the vertical velocity
!>          it would have produced.
!> @param[in]     nlayers       Number of layers in the mesh
!> @param[in,out] dtheta_out    Output theta increment after diabatic mixing
!> @param[in]     theta         Input theta before diabatic heating
!> @param[in]     rho           Denisty at wtheta points
!> @param[in]     dA            Cell face areas on w2 mesh
!> @param[in]     dAsh          Cell face areas on shifted w2 mesh
!> @param[in]     dtheta_in     Input theta increment which generates mixing
!> @param[in]     ndf_wth       Number of degrees of freedom per cell
!> @param[in]     undf_wth      Number of unique degrees of freedom per cell
!> @param[in]     map_wth       Cell dofmap for wth fields
!> @param[in]     ndf_w2        Number of degrees of freedom per cell
!> @param[in]     undf_w2       Number of unique degrees of freedom per cell
!> @param[in]     map_w2        Cell dofmap for wth fields
!> @param[in]     ndf_w2sh      Number of degrees of freedom per cell
!> @param[in]     undf_w2sh     Number of unique degrees of freedom per cell
!> @param[in]     map_w2sh      Cell dofmap for wth fields
subroutine wdiab_theta_code( nlayers,                               &
                             dtheta_out,                            &
                             theta,                                 &
                             rho,                                   &
                             dA,                                    &
                             dAsh,                                  &
                             dtheta_in,                             &
                             ndf_wth,                               &
                             undf_wth,                              &
                             map_wth,                               &
                             ndf_w2,                                &
                             undf_w2,                               &
                             map_w2,                                &
                             ndf_w2sh,                              &
                             undf_w2sh,                             &
                             map_w2sh)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
  integer(kind=i_def), dimension(ndf_wth),  intent(in)  :: map_wth
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), dimension(ndf_w2),  intent(in)  :: map_w2
  integer(kind=i_def), intent(in) :: ndf_w2sh, undf_w2sh
  integer(kind=i_def), dimension(ndf_w2sh),  intent(in)  :: map_w2sh

  real(kind=r_def), dimension(undf_wth),  intent(inout) :: dtheta_out
  real(kind=r_def), dimension(undf_wth),  intent(in)    :: theta, dtheta_in, &
                                                           rho
  real(kind=r_def), dimension(undf_w2), intent(in) :: dA
  real(kind=r_def), dimension(undf_w2sh), intent(in) :: dAsh

  ! Internal variables
  integer(kind=i_def) :: k
  integer(kind=i_def), parameter :: lev_pad = 1

  real(kind=r_def) :: lapse, mix_coeff, rho_flux
  real(kind=r_def) :: flux(nlayers)

  do k = 1, nlayers
    flux(k) = 0.0_r_def
  end do

  ! run from levels 2 to model_levels-1
  do k = 1+lev_pad, nlayers-lev_pad

    ! if theta increment is a heating
    if (dtheta_in(map_wth(1)+k) > 0.0_r_def) then

      ! calculate lapse rate
      lapse = theta(map_wth(1)+k+1) - theta(map_wth(1)+k)

      ! calculate mixing coefficient
      if (lapse > 2.0_r_def * dtheta_in(map_wth(1)+k)) then
        mix_coeff = dtheta_in(map_wth(1)+k) / lapse
      else
        mix_coeff = 0.5_r_def
      end if

      ! calculate flux into level above
      rho_flux = 0.5_r_def * (rho(map_wth(1)+k) + rho(map_wth(1)+k+1))
      flux(k+1) = flux(k+1) + mix_coeff * dtheta_in(map_wth(1)+k) &
                * rho_flux * dAsh(map_w2sh(B)+k+1)

    ! if theta increment is a cooling
    else if (dtheta_in(map_wth(1)+k) < 0.0_r_def) then

      ! calculate lapse rate
      lapse = theta(map_wth(1)+k-1) - theta(map_wth(1)+k)

      ! calculate mixing coefficient
      if (lapse < 2.0_r_def * dtheta_in(map_wth(1)+k)) then
        mix_coeff = dtheta_in(map_wth(1)+k) / lapse
      else
        mix_coeff = 0.5_r_def
      end if

      ! calculate flux into level below
      rho_flux = 0.5_r_def * (rho(map_wth(1)+k-1) + rho(map_wth(1)+k))
      flux(k) = flux(k) - mix_coeff * dtheta_in(map_wth(1)+k) &
              * rho_flux * dAsh(map_w2sh(B)+k)

    end if

  end do

  ! calculate mixing increment from flux
  do k = 0, lev_pad
    dtheta_out(map_wth(1)+k) = - flux(k+1) & ! flux(k) = 0
                             / (rho(map_wth(1)+k) * dA(map_w2(B)+k))
  end do
  do k = 1+lev_pad, nlayers-lev_pad
    dtheta_out(map_wth(1)+k) = (flux(k) - flux(k+1)) &
                             / (rho(map_wth(1)+k) * dA(map_w2(B)+k))
  end do
  do k = nlayers-lev_pad+1, nlayers
    dtheta_out(map_wth(1)+k) = flux(k) & ! flux (k+1) = 0
                             / (rho(map_wth(1)+k) * dA(map_w2(B)+k))
  end do

end subroutine wdiab_theta_code

end module wdiab_theta_kernel_mod
