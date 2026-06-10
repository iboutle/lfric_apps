!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief
!>
module multi_to_w2_kernel_mod

  use argument_mod,      only: arg_type,          &
                               GH_FIELD, GH_REAL, &
                               GH_WRITE, GH_READ, &
                               CELL_COLUMN,       &
                               ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only: r_def, i_def
  use fs_continuity_mod, only: W2
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: multi_to_w2_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,  W2)      &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: multi_to_w2_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: multi_to_w2_code

contains

!> @brief The subroutine which is called directly by the psy layer
!! @param[in] nlayers     Integer the number of layers
!! @param[in] multi_u     Multi-data wind field
!! @param[in,out] w2_u    W2 wind field
!! @param[in] ndf_mult    The number of degrees of freedom per cell for mult
!! @param[in] undf_mult   The number of unique degrees of freedom for mult
!! @param[in] map_mult    Integer array holding the dofmap for the cell at the
!>                        base of the column for mult
!! @param[in] ndf_w2      The number of degrees of freedom per cell for w2
!! @param[in] undf_w2     The number of unique degrees of freedom for w2
!! @param[in] map_w2      Integer array holding the dofmap for the cell at the
!>                        base of the column for w2
subroutine multi_to_w2_code( nlayers,                       &
                             multi_u,                       &
                             w2_u,                          &
                             ndf_mult, undf_mult, map_mult, &
                             ndf_w2, undf_w2, map_w2        &
                             )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_mult, undf_mult
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2

  real(kind=r_def), dimension(undf_mult), intent(in) :: multi_u
  real(kind=r_def), dimension(undf_w2),  intent(inout)    :: w2_u
  integer(kind=i_def), dimension(ndf_mult),  intent(in)    :: map_mult
  integer(kind=i_def), dimension(ndf_w2),   intent(in)    :: map_w2

  ! Internal variables
  integer(kind=i_def) :: k, df

  do k = 0, nlayers-1
    do df = 1, 6

      w2_u(map_w2(df) + k ) = multi_u(map_mult(1) + (df-1)*nlayers + k)

    end do
  end do

end subroutine multi_to_w2_code

end module multi_to_w2_kernel_mod
