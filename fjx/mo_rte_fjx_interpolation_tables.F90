module mo_rte_fjx_interpolation_tables 
  use mo_rte_kind, only: wp, wl
  implicit none
  private 

  public :: ty_rte_fjx_interp_table, & 
            interp_x, make_rte_fjx_table

  type ty_rte_fjx_interp_table
   !private
    character(len=32) :: title 
    character(len=32) :: reaction
    character(len=64) :: source
    ! Will be allocated as ngpt x {1, 2, or 3}
    real(wp), dimension(:,:), allocatable :: cross_section
    ! Will  be allocated as 1, 2, or 3 
    real(wp), dimension(:),   allocatable :: coordinate !temperature or pressure
    logical(wl) :: if_p !true if interpolation is over pressure coordinates
    logical(wl) :: if_x
  end type ty_rte_fjx_interp_table

contains
  ! --------------------------------------------------------------------------------------------------
  function interp_x(x, ty_table) result(values)
    real(wp),                   intent(in) :: x ! Pressure or temperature at which values are desired 
    type(ty_rte_fjx_interp_table), intent(in) :: ty_table
    real(wp), dimension(size(ty_table%cross_section,1)) :: values 

    select case(size(ty_table%cross_section,2))
      case(1)
         values(:) = ty_table%cross_section(:,1)
      case(2)
        values =linear_interp(x, ty_table%coordinate(1), ty_table%coordinate(2), ty_table%cross_section)
      case(3)
        if (x < ty_table%coordinate(2)) then 
          values =linear_interp(x, ty_table%coordinate(1), ty_table%coordinate(2), ty_table%cross_section(:, 1:2))
        else 
          values =linear_interp(x, ty_table%coordinate(2), ty_table%coordinate(3), ty_table%cross_section(:, 2:3))
        end if
    end select
  end function interp_x

  function linear_interp(x, x1, x2, table)
    real(wp),                 intent(in) :: x, x1, x2 
    real(wp), dimension(:,:), intent(in) :: table
    real(wp), dimension(size(table,1))   :: linear_interp

    real(wp)                             :: fact

    ! Compute weighting, provide do interpolation
    fact=max(0._wp, min(1._wp, (x-x1)/(x2-x1)))
    linear_interp=table(:,1)+fact*(table(:,2)-table(:,1))
  end function linear_interp
  ! --------------------------------------------------------------------------------------------------
  function make_rte_fjx_table(table,title,reaction,source,cross_section,coordinate,if_p,if_x)  result(error_msg)
    type(ty_rte_fjx_interp_table), intent(inout)  :: table
    character(len=32), intent(in) :: title, reaction
    character(len=64), intent(in) :: source
    ! Will be allocated as ngpt x {1, 2, or 3}
    real(wp), dimension(:,:), intent(in) :: cross_section
    ! Will  be allocated as 1, 2, or 3 
    real(wp), dimension(:),   intent(in) :: coordinate !temperature or pressure
    logical(wl), intent(in) :: if_p !true if interpolation is over pressure coordinates
    logical(wl), intent(in) :: if_x
    character(len=32) :: error_msg  ! intent(out) gives compiler error, symbol is not a dummy variable

    error_msg="" !redundant, I guess.

    table%title=title
    table%reaction=reaction
    table%source=source
    table%cross_section=cross_section
    table%coordinate=coordinate
    table%if_p=if_p
    table%if_x=if_x

  end function make_rte_fjx_table
  ! --------------------------------------------------------------------------------------------------
end module mo_rte_fjx_interpolation_tables
