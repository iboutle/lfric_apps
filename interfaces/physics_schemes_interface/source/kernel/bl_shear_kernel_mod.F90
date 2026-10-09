!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculate the vertical wind shear on wtheta levels for use in the
!>        boundary layer scheme.
module bl_shear_kernel_mod

  use argument_mod,      only : arg_type,                  &
                                GH_FIELD, GH_REAL,         &
                                GH_READ, GH_WRITE,         &
                                CELL_COLUMN,               &
                                ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only : i_def, r_def, eps
  use fs_continuity_mod, only : Wtheta, W2
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  type, public, extends(kernel_type) :: bl_shear_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta), &! shear
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),     &! u_physics
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta)  &! dz_wth
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: bl_shear_code
  end type

  public :: bl_shear_code

contains

  !> @brief Compute the vertical shear of the horizontal wind on wtheta levels
  !> @param[in]     nlayers    Number of layers
  !> @param[in,out] shear      Vertical wind shear on wtheta levels
  !> @param[in]     u_physics  3D physical wind in W2 space
  !> @param[in]     dz_wth     Layer depths at wtheta points
  !> @param[in]     ndf_wth    Number of DOFs per cell for wtheta space
  !> @param[in]     undf_wth   Number of unique DOFs for wtheta space
  !> @param[in]     map_wth    Dofmap for the cell at the base of the column for wtheta space
  !> @param[in]     ndf_w2     Number of DOFs per cell for W2 space
  !> @param[in]     undf_w2    Number of unique DOFs for W2 space
  !> @param[in]     map_w2     Dofmap for the cell at the base of the column for W2 space
  subroutine bl_shear_code(nlayers,                    &
                           shear,                      &
                           u_physics,                  &
                           dz_wth,                     &
                           ndf_wth, undf_wth, map_wth, &
                           ndf_w2, undf_w2, map_w2)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth)
    integer(kind=i_def), intent(in) :: map_w2(ndf_w2)

    real(kind=r_def), dimension(undf_wth), intent(inout) :: shear
    real(kind=r_def), dimension(undf_w2),  intent(in)    :: u_physics
    real(kind=r_def), dimension(undf_wth), intent(in)    :: dz_wth

    ! Local variables
    integer(kind=i_def) :: k
    real(kind=r_def)    :: idz, ssq

    shear(map_wth(1)) = 0.0_r_def

    do k = 1, nlayers - 1
      idz = (1.0_r_def / dz_wth(map_wth(1) + k))**2
      ssq = ((u_physics(map_w2(1) + k) - u_physics(map_w2(1) + k-1))**2 + &
             (u_physics(map_w2(2) + k) - u_physics(map_w2(2) + k-1))**2 + &
             (u_physics(map_w2(3) + k) - u_physics(map_w2(3) + k-1))**2 + &
             (u_physics(map_w2(4) + k) - u_physics(map_w2(4) + k-1))**2 ) &
            / 2.0_r_def
      shear(map_wth(1) + k) = sqrt(ssq * idz + eps)
    end do

    ! Horizontal W2 dofs only exist for layers 0 to nlayers-1
    shear(map_wth(1) + nlayers) = 0.0_r_def

  end subroutine bl_shear_code

end module bl_shear_kernel_mod
