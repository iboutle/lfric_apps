!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Activate number due to convection
!>
module casim_simple_act_kernel_mod

  use argument_mod,       only : arg_type,              &
                                 GH_FIELD, GH_REAL,     &
                                 GH_READ, GH_READWRITE, &
                                 GH_WRITE, CELL_COLUMN, &
                                 GH_SCALAR

  use constants_mod,      only : r_def, r_double, i_def, i_um, r_um
  use fs_continuity_mod,  only : Wtheta
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: casim_simple_act_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                   &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),   & ! dmx_conv
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),   & ! nx_mphys
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),   & ! cf_x
         arg_type(GH_SCALAR, GH_REAL, GH_READ),           & ! fixed_num
         arg_type(GH_SCALAR, GH_REAL, GH_READ)            & ! tidy_num
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: casim_simple_act_code
  end type

  public :: casim_simple_act_code

contains

  !> @details Activate number to a fixed in-cloud value
  !> @param[in]     nlayers       Number of layers
  !> @param[in]     dmx_conv      Increment from convection
  !> @param[in,out] nx_mphys      CASIM number concentration
  !> @param[in]     cf_x          Cloud fraction
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth      Number unique of degrees of freedom for potential temperature space
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column for potential temperature space
  subroutine casim_simple_act_code(nlayers,      &
                                   dmx_conv,     &
                                   nx_mphys,     &
                                   cf_x,         &
                                   fixed_num,    &
                                   tidy_num,     &
                                   ndf_wth,      &
                                   undf_wth,     &
                                   map_wth)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)     :: nlayers
    integer(kind=i_def), intent(in)     :: ndf_wth
    integer(kind=i_def), intent(in)     :: undf_wth

    integer(kind=i_def), intent(in),    dimension(ndf_wth)  :: map_wth

    real(kind=r_def),    intent(in), dimension(undf_wth) :: dmx_conv
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: nx_mphys
    real(kind=r_def),    intent(in), dimension(undf_wth) :: cf_x
    real(kind=r_def),    intent(in) :: fixed_num, tidy_num

    ! Local variables for the kernel
    integer(i_um) :: k

    do k = 1, nlayers

      ! If there is mass but no number
      if (dmx_conv(map_wth(1) + k) > 0.0_r_def .and. &
           nx_mphys(map_wth(1) + k) <= tidy_num) then

        ! scale to grid-box mean
        nx_mphys( map_wth(1) + k) = fixed_num &
                                    * cf_x(map_wth(1) + k)
      end if
    end do

    ! Save value of nx_mphys at level 1 for level 0 increment
    nx_mphys( map_wth(1) + 0) = nx_mphys( map_wth(1) + 1)

  end subroutine casim_simple_act_code

end module casim_simple_act_kernel_mod
