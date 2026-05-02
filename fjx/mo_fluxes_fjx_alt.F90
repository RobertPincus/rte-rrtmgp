! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and George Milly
!
! Copyright 2022-  Trustees of Columbia University. All rights reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
!> ## Compute output quantities from spectrally-resolved flux profiles
!>
!>    This module contains an alternate version of mo_fluxes_fjx configured to receive rather than calculate
!>    the actinic flux values used in the calculation of photolysis rates.  This is primarily to maintain
!>    compatibility with testfluxes, which tests the photolysis calculation using CloudJ77Prather fff data.
!
! -------------------------------------------------------------------------------------------------
module mo_fluxes_fjx_alt
  use mo_rte_kind,       only: wp
  use mo_rte_config,     only: check_extents
  use mo_rte_util_array_validation, only: extents_are
  use mo_optical_props,  only: ty_optical_props
  use mo_fluxes,         only: ty_fluxes
  use mo_gas_optics_util_string
  use mo_rte_fjx_interpolation_tables
  implicit none
  private

  ! -----------------------------------------------------------------------------------------------
  type, extends(ty_fluxes), public :: ty_fluxes_fjx_alt
    private
    !probably shouldn't be private.. these are the "output" of reduce function, but not in arg list
    !our outputs will be the photolysis rates
    !this is our reduce function

    !extra input args needed for jratet
    real(wp), dimension(:,:), allocatable, public :: play ,tlay
    !extra output arg needed for jratet
    real(wp), dimension(:,:,:), pointer, public :: valjl => NULL()
    !real(wp), dimension(:,:,:), allocatable, public :: valjl

    ! QQQ: Supplied cross sections in each wavelength bin (cm2)
    type(ty_rte_fjx_interp_table), dimension(:), allocatable :: QQQ  ! cross-sections, also name and t coords
    integer :: o1d_loc, o3_loc

  contains
    procedure, public :: reduce      => calc_jrate
    procedure, public :: are_desired => are_desired_fjx
    procedure, public :: load
    procedure, public :: set_temp_press
    procedure, public :: get_njx
    procedure, public :: get_react_names
    procedure, private :: jratet
  end type ty_fluxes_fjx_alt
  ! -----------------------------------------------------------------------------------------------

  ! -----------------------------------------------------------------------------------------------
  !
contains
  ! --------------------------------------------------------------------------------------
  !
  !> calculate jrate
  !
  ! --------------------------------------------------------------------------------------
  function calc_jrate(this, gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir) result(error_msg)
    class(ty_fluxes_fjx_alt),          intent(inout) :: this
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_up ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_dn ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    class(ty_optical_props),           intent(in   ) :: spectral_disc  !< derived type with spectral information
    logical,                           intent(in   ) :: top_at_1
    real(kind=wp), dimension(:,:,:), optional, &
                                       intent(in   ) :: gpt_flux_dn_dir! Direct flux down
    character(len=128)                               :: error_msg
    ! ------
    integer :: ncol, nlev, ngpt, C,L
    real(wp), dimension(size(gpt_flux_up,dim=1),size(gpt_flux_up,dim=2),size(gpt_flux_up,dim=3)) :: fff
    real(wp) :: localpress
    ! ------
    ncol = size(gpt_flux_up, DIM=1)
    nlev = size(gpt_flux_up, DIM=2)
    ngpt = size(gpt_flux_up, DIM=3)
    error_msg = ""

    if(.not. associated(this%valjl)) then
      error_msg="must provide variable to hold valjl (photolysis rate) output"
      return
    endif
    if(check_extents) then
      !
      ! Check array sizes
      !  Input arrays
      !
      if(.not. extents_are(gpt_flux_dn, ncol, nlev, ngpt)) &
        error_msg = "reduce: gpt_flux_dn array incorrectly sized"

      if(present(gpt_flux_dn_dir)) then
        if(.not. extents_are(gpt_flux_dn_dir, ncol, nlev, ngpt)) &
          error_msg = "reduce: gpt_flux_dn_dir array incorrectly sized"
      end if
      !
      ! Output arrays
      !
      print *, "shape gpt_flux up:", shape(gpt_flux_up)
      if(.not. extents_are(this%valjl, ncol, nlev-1, this%get_njx())) &
        error_msg = "reduce: valjl array incorrectly sized"

      if(error_msg /= "") return
    end if

    !fff is actinic flux, some combination of gpt_flux_up & gpt_flux_down
    !play and tlay are "new" inputs needed, passed in to class above
    !will add valjl to class too, since seems like reduce function is supposed to give ultimate output
    !chose to just put fff in this function, since it will be calculated here
    !could still basically put jratet inside of this function, right? not sure best organization

    !for now..
    fff=gpt_flux_up
    !possibly... fff=2._wp*(gpt_flux_up+gpt_flux_down)+gpt_flux_dn_dir if this is what Manners is saying?

    do L=1,nlev-1  !be consistent, lev or lay... also, should this be nlay-1?
      do C=1,ncol
        localpress=this%play(C,L)
        !must zero bin-11 (216-222 & 287-291 nm) below 100 hPa since O2 e-fold is too weak
        if (localpress .gt. 100._wp) then
          fff(C,L,11) = 0._wp
        endif
      enddo
    enddo

    error_msg=this%jratet(this%play,this%tlay,fff)
    if(error_msg /="") return

  end function calc_jrate
  ! --------------------------------------------------------------------------------------
  !
  !> Are any fluxes desired from this set of g-point fluxes? We can tell because memory will
  !>   be allocated for output
  !
  ! --------------------------------------------------------------------------------------
  function are_desired_fjx(this)
    class(ty_fluxes_fjx_alt), intent(in   )   :: this
    logical                                   :: are_desired_fjx

   !are_desired_fjx = any( [associated(this%flux_up),     &
   !                              associated(this%flux_dn),     &
   !                              associated(this%flux_dn_dir), &
   !                              associated(this%flux_net)] )
  end function are_desired_fjx
  ! --------------------------------------------------------------------------------------
  ! calculate photolysis rates
  ! this should be its own function, right?
  ! is spectral_reducer just used to get flux in form desired by jratet,
  !  from form it leaves something like mo_rte_sw.F90 in?
  ! --------------------------------------------------------------------------------------
  function jratet(this,play,tlay,fff)  result(error_msg)
    class(ty_fluxes_fjx_alt), intent(inout)       :: this
    real(kind=wp), dimension(:,:), intent(in   ) :: play,tlay
    real(kind=wp), dimension(:,:,:), intent(in   ) :: fff
    character(len=128)           :: error_msg     !! Empty if successful 
    ! Local variables
    integer :: ncol,nlay,L,J,K,C
    real(wp) :: localtemp,localpress
    real(wp), dimension(this%get_njx()) :: valj
    real(wp), dimension(size(fff,dim=3)) :: QQQTvals,XQO2vals,XQO3vals,XQ1Dvals,QO31D
    ncol = size(play,dim=1)
    nlay = size(play,dim=2)
    error_msg = ""

!-----------------------------------------------------------------------
! in:
!        PPJ(L_+1) = pressure profile at edges
!        TTJ(L_+1) = = temperatures at mid-level
!        FFF(K=1:NW, L=1:L_) = mean actinic flux
! out:
!        VALJL(L_,JX_)  JX_ = no of dimensioned J-values in CTM code
!-----------------------------------------------------------------------

    do L=1,nlay-1
      do C=1,ncol
        ! need temperature, pressure, and density at mid-layer (for some quantum yields):
        ! Prather calculated mid-layer pressure here, but we do so elsewhere
        localtemp=tlay(C,L)
        localpress=play(C,L)

        do J = 1,this%get_njx()
          if (J==this%o1d_loc) then
            XQO3vals=interp_x(localtemp, this%QQQ(this%o3_loc))
            XQ1Dvals=interp_x(localtemp, this%QQQ(J))
            QO31D=XQ1Dvals*XQO3vals
            VALJ(J) = dot_product(QO31D,FFF(C,L,:))
          else
            ! need to allow for Pressure interpolation
            if (this%QQQ(J)%if_p) then
              QQQTvals=interp_x(localpress, this%QQQ(J))
            else
             !print *, 'shape this QQQ(J)', shape(this%QQQ(J)) not yet defined, I think
              QQQTvals=interp_x(localtemp, this%QQQ(J))
            endif
            VALJ(J) = dot_product(QQQTvals,FFF(C,L,:))
          endif
          this%valjl(C,L,J)=valj(J)  ! should we eliminate valj and just fill valjl directly?
        enddo

      enddo
    enddo
    if(error_msg /="") return

  end function jratet
  ! --------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------
  !
  !
  ! Initialization
  !
  !--------------------------------------------------------------------------------------------------------------------

  function load (this, in_table, reacts_list) result(error_msg)
    class(ty_fluxes_fjx_alt), intent(inout) :: this
    type(ty_rte_fjx_interp_table), dimension(:), intent(in) :: in_table
    character(len=*), dimension(:), optional, target, intent(in) :: reacts_list
    character(len=128)           :: error_msg     !! Empty if successful

    character(len=32), dimension(:), pointer :: reactions
    integer nreacts, n, nn
    error_msg = ""

    ! should we add a check that no reaction is duplicated in reacts_list?

    if (present(reacts_list)) then
      ! check that if o1d is requested that o3 is as well
      if (string_in_array('O3(1D)',reacts_list) .and. .not. string_in_array('O3',reacts_list)) then
        error_msg='use of O3(1D) requires use of O3'
        return
      endif
      reactions => reacts_list
    else
      allocate (reactions(size(in_table)))
      do n=1,size(in_table)
        reactions(n)=in_table(n)%title
      enddo
    endif
    nreacts=size(reactions)
    allocate (this%QQQ(nreacts)) !is this consistent with shape(this%QQQ(J)) being empty in jratet above? 

    do n=1,nreacts
      if(present(reacts_list)) then
        ! check that elements of reacts_list are in in_table titles
        if (.not. string_in_array(reactions(n),in_table%title)) then
          error_msg=reactions(n) // " not found in available reactions tables"
          return
        endif
        ! find location of reaction n in in_table - nn from string_loc_in_array
        ! case insensitive due to string_loc_in_array use of lower_case
        nn=string_loc_in_array(reactions(n),in_table%title) 
        this%QQQ(n)=in_table(nn)
      else
        this%QQQ(n)=in_table(n)
      end if
      if (reactions(n)=='O3(1D)') this%o1d_loc=n
      if (reactions(n)=='O3') this%o3_loc=n
    end do    

  end function load 
  !--------------------------------------------------------------------------------------------------------------------

  function set_temp_press(this,temp,press) result(error_msg)
    class(ty_fluxes_fjx_alt), intent(inout) :: this
    real(wp), dimension(:,:), intent(in) :: temp,press
    character(len=128)           :: error_msg     !! Empty if successful 

    error_msg = "" ! jratet crashed without initizializing this, but this didn't

    if (minval(temp) < 0._wp) then 
      error_msg="all temperature values must be greater than zero"
      return
    endif

    if (minval(press) < 0._wp) then
      error_msg="all pressure values must be greater than zero"
      return
    endif

    this%tlay = temp    ! allocate on assignment 
    this%play = press

   !print *, "shape this%tlay", shape(this%tlay)

  end function set_temp_press
  !--------------------------------------------------------------------------------------------------------------------

  !
  ! Inquiry functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  
  pure function get_njx (this)
    class(ty_fluxes_fjx_alt), intent(in) :: this
    integer :: get_njx

    if (.not. allocated(this%QQQ)) then
      get_njx=0
    else  
      get_njx=size(this%QQQ)
    end if

  end function get_njx        
  !--------------------------------------------------------------------------------------------------------------------

  function get_react_names(this)
    class(ty_fluxes_fjx_alt), intent(in) :: this
    character(len=32), dimension(:), allocatable :: get_react_names

    if (.not. allocated(this%QQQ)) then
      get_react_names=''
    else
      get_react_names=this%QQQ%title
    end if

  end function get_react_names
  !--------------------------------------------------------------------------------------------------------------------

end module mo_fluxes_fjx_alt

