!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates Uzig
!> @details Uzig is given a value of 1 if the value of u (horizontal velocity)
!!          changes value by > 50m/s over 2 vertical levels. Seeing this
!!          happen at a range of vertical levels can show oscillations
!!          in the wind field, which are not physical
module calc_uzig_kernel_mod

  use argument_mod,      only: arg_type,          &
                               GH_FIELD, GH_REAL, &
                               GH_READ,           &
                               GH_READWRITE,      &
                               CELL_COLUMN
  use constants_mod,     only: r_def, i_def
  use fs_continuity_mod, only: W3, W2
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: calc_uzig_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),      & !u
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  W3) & !uzig
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: calc_uzig_code
  end type calc_uzig_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: calc_uzig_code

contains
!> @brief The identifies values of u_in_w2 greater than 50 and sets them to 1
!>        in the uzig field.
!! @param[in] nlayers Integer the number of layers
!! @param[in] u_in_w2 Real array, w component of u_physics
!! @param[in,out] uzig - the uzig field
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] undf_w3 The number of unique degrees of freedom for w3
!! @param[in] map_w3 Integer array holding the dofmap for the cell at the
!>            base of the column for w3
subroutine calc_uzig_code(nlayers,                 &
                          u_in_w2,                 &
                          uzig,                    &
                          ndf_w2, undf_w2, map_w2, &
                          ndf_w3, undf_w3, map_w3  &
                          )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3

  real(kind=r_def), dimension(undf_w2), intent(in)    :: u_in_w2
  real(kind=r_def), dimension(undf_w3), intent(inout) :: uzig
  integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3

  ! Internal variables
  integer(kind=i_def) :: k, df

  do k = 1, nlayers-2
    do df = 1, 4
      if ( abs(u_in_w2(map_w2(df)+k+1)-u_in_w2(map_w2(df)+k-1)) > 50.0_r_def ) then

        uzig(map_w3(1)+k) = uzig(map_w3(1)+k) + 1.0_r_def

      else

        uzig(map_w3(1)+k) = uzig(map_w3(1)+k) + 0.0_r_def

      end if
    end do
  end do

end subroutine calc_uzig_code

end module calc_uzig_kernel_mod
