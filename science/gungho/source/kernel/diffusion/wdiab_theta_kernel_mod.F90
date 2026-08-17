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
                                    CELL_COLUMN
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : Wtheta
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: wdiab_theta_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                    &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta),                   &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                   &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta)                    &
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
!> @param[in]     dtheta_in     Input theta increment which generates mixing
!> @param[in]     ndf_wth       Number of degrees of freedom per cell
!> @param[in]     undf_wth      Number of unique degrees of freedom per cell
!> @param[in]     map_wth       Cell dofmap for wth fields
subroutine wdiab_theta_code( nlayers,                               &
                             dtheta_out,                            &
                             theta,                                 &
                             dtheta_in,                             &
                             ndf_wth,                               &
                             undf_wth,                              &
                             map_wth )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
  integer(kind=i_def), dimension(ndf_wth),  intent(in)  :: map_wth

  real(kind=r_def), dimension(undf_wth),  intent(inout) :: dtheta_out
  real(kind=r_def), dimension(undf_wth),  intent(in)    :: theta, dtheta_in

  ! Internal variables
  integer(kind=i_def) :: k
  integer(kind=i_def), parameter :: lev_pad = 1

  real(kind=r_def) :: lapse, mix_coeff

  ! levels 0 and 1 output increment is just input increment but can mix
  ! down into here
  do k = 0, lev_pad
    dtheta_out(map_wth(1)+k) = dtheta_in(map_wth(1)+k)
  end do
  ! top level output increment is just input increment but can mix up
  ! into here
  do k = nlayers-lev_pad+1, nlayers
    dtheta_out(map_wth(1)+k) = dtheta_in(map_wth(1)+k)
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

      ! calculate mixing increment at level above
      dtheta_out(map_wth(1)+k+1) = dtheta_out(map_wth(1)+k+1) &
                                      + mix_coeff * dtheta_in(map_wth(1)+k)

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

      ! calculate mixing increment at level below
      dtheta_out(map_wth(1)+k-1) = dtheta_out(map_wth(1)+k-1) &
                                      + mix_coeff * dtheta_in(map_wth(1)+k)

    else

      mix_coeff = 0.0_r_def

    end if

    ! calculate mixing increment at current level
    dtheta_out(map_wth(1)+k) = dtheta_out(map_wth(1)+k) &
                                  + (1.0_r_def - mix_coeff) &
                                  * dtheta_in(map_wth(1)+k)

  end do

end subroutine wdiab_theta_code

end module wdiab_theta_kernel_mod
