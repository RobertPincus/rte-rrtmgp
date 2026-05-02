! This code is part of rte-fjx
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
!>    This module contains an abstract class and a broadband implmentation that sums over all spectral points
!>    The abstract base class defines the routines that extenstions must implement: `reduce()` and `are_desired()`
!>    The intent is for users to extend it as required, using mo_flxues_broadband as an example
!
! -------------------------------------------------------------------------------------------------
module mo_fluxes_fjx
  use mo_rte_kind,       only: wp
  use mo_rte_config,     only: check_extents
  use mo_rte_util_array_validation, only: extents_are
  use mo_optical_props,  only: ty_optical_props, ty_optical_props_arry
  use mo_fluxes,         only: ty_fluxes
  use mo_gas_optics_util_string
  use mo_rte_fjx_interpolation_tables
  implicit none
  private

  ! -----------------------------------------------------------------------------------------------
  type, extends(ty_fluxes), public :: ty_fluxes_fjx
    private
    !probably shouldn't be private.. these are the "output" of reduce function, but not in arg list
    !our outputs will be the photolysis rates
    !this is our reduce function

    !extra input args needed for jratet
    real(wp), dimension(:,:), allocatable :: play ,tlay
    !extra output arg needed for jratet
    real(wp), dimension(:,:,:), pointer, public :: valjl => NULL()
    !real(wp), dimension(:,:,:), allocatable, public :: valjl


    !extra input needed for actinic flux calc
    real(wp), dimension(:,:,:), allocatable :: tau
    !mimic mo_fluxes_bygpoint to pass actinic flux out
    real(wp), dimension(:,:,:), pointer, public :: gpt_act_flux => NULL()  !(ncol,nlay,nband)

    !---- Variables in file 'FJX_spec.dat' (RD_XXX)
    !> will we use any of these?
    ! WL: Centres of wavelength bins - 'effective wavelength'  (nm)
   !real(wp)  WL(S_)
    ! WBIN: Boundaries of wavelength bins                  (microns)
   !real(wp)  WBIN(S_+1)
    ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
    ! FW: Solar flux in W/m2
    ! FP: PAR quantum action spectrum
   !real(wp)  FL(S_),FW(S_),FP(S_)

    ! QQQ: Supplied cross sections in each wavelength bin (cm2)
    type(ty_rte_fjx_interp_table), dimension(:), allocatable :: QQQ  ! cross-sections, also name and t coords
    integer :: o1d_loc, o3_loc

  contains
    procedure, public :: reduce      => calc_jrate
    procedure, public :: are_desired => are_desired_fjx
    procedure, public :: load
    procedure, public :: set_temp_press
    procedure, public :: set_tau
    procedure, public :: get_njx
    procedure, public :: get_react_names
    procedure, private :: jratet
  end type ty_fluxes_fjx
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
    class(ty_fluxes_fjx),              intent(inout) :: this
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_up ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_dn ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    class(ty_optical_props),           intent(in   ) :: spectral_disc  !< derived type with spectral information
    logical,                           intent(in   ) :: top_at_1
    real(kind=wp), dimension(:,:,:), optional, &
                                       intent(in   ) :: gpt_flux_dn_dir! Direct flux down
    character(len=128)                               :: error_msg
    ! ------
    integer :: ncol, nlev, ngpt, C,L,K
    real(wp), dimension(size(gpt_flux_up,dim=1),size(gpt_flux_up,dim=2)-1,size(gpt_flux_up,dim=3)) :: fff
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
    if(.not. associated(this%gpt_act_flux)) then ! thinking for now that fff is only output for testing
      error_msg="must provide variable to hold fff (actinic flux) output"
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
      if(.not. extents_are(this%gpt_act_flux, ncol, nlev-1, ngpt)) &
        error_msg = "reduce: gpt_act_flux array incorrectly sized"

      if(error_msg /= "") return
    end if

    !fff is actinic flux
    !calculate from energetic fluxes
    do L=1,nlev-1  !this is nlay
      do C=1,ncol
        localpress=this%play(C,L)
        do K=1,ngpt
          if (this%tau(C,L,K) .lt. sqrt(epsilon(1._wp))) then
            !Manners' Method 1
            fff(C,L,K) = 2._wp*(gpt_flux_up(C,L,K)+gpt_flux_dn(C,L,K))-gpt_flux_dn_dir(C,L,K) !dn includes dir, remove factor of 2
          else !assuming for now toa_at_1 false, as it is in our test data
            fff(C,L,K) = ( (gpt_flux_dn(C,L+1,K)-gpt_flux_up(C,L+1,K)) - & 
                           (gpt_flux_dn(C,L  ,K)-gpt_flux_up(C,L,  K)) ) / &
              this%tau(C,L,K) !Manners' Method 2, total flux divergence div by tau in layer
          endif
        enddo
        !must zero bin-11 (216-222 & 287-291 nm) below 100 hPa since O2 e-fold is too weak
        if (localpress .gt. 100._wp) fff(C,L,11) = 0._wp
      enddo
    enddo
    print *, gpt_flux_up(1,:,1)
    print *, "shape gpt_flux dn:", shape(gpt_flux_dn)
    print *, gpt_flux_dn(1,:,1)
    print *, "shape gpt_flux dn_dir:", shape(gpt_flux_dn_dir)
    print *, gpt_flux_dn_dir(1,:,1)

    this%gpt_act_flux=fff

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
    class(ty_fluxes_fjx), intent(in   )       :: this
    logical                                   :: are_desired_fjx

   !are_desired_fjx = any( [associated(this%flux_up),     &
   !                              associated(this%flux_dn),     &
   !                              associated(this%flux_dn_dir), &
   !                              associated(this%flux_net)] )
  !are_desired_fjx=.true.
  end function are_desired_fjx
  ! --------------------------------------------------------------------------------------
  ! calculate photolysis rates
  ! --------------------------------------------------------------------------------------
  function jratet(this,play,tlay,fff)  result(error_msg)
    class(ty_fluxes_fjx), intent(inout)       :: this
    real(kind=wp), dimension(:,:), intent(in   ) :: play,tlay
    real(kind=wp), dimension(:,:,:), intent(in   ) :: fff
    character(len=128)           :: error_msg     !! Empty if successful 
    ! Local variables
    integer :: ncol,nlay,L,J,C
    real(wp) :: localtemp,localpress
    real(wp), dimension(this%get_njx()) :: valj
    real(wp), dimension(size(fff,dim=3)) :: QQQTvals,XQO3vals,XQ1Dvals,QO31D
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
!     real*8, intent(in)  ::  PPJ(LU+1),TTJ(LU+1)
!     real*8, intent(inout)  ::  FFF(W_,LU)
!     real*8, intent(out), dimension(LU,NJXU) ::  VALJL

! we can just use mid-layer pressure, as rte-fjx already requires this elsewhere (calculated in testcase.F90)
    do L=1,nlay-1
      do C=1,ncol
        ! need temperature, pressure, and density at mid-layer (for some quantum yields):
        localtemp=tlay(C,L)
        localpress=play(C,L)
!       if (L .eq. 1) then
!         PP = PPJ(1)
!       else
!         PP  = (PPJ(L)+PPJ(L+1))*0.5d0
!       endif

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
    class(ty_fluxes_fjx), intent(inout) :: this
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
        reactions(n)=in_table(n)%title  ! will this ever be used?
      enddo
    endif
    nreacts=size(reactions)
    allocate (this%QQQ(nreacts))

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
    class(ty_fluxes_fjx), intent(inout) :: this
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

  function set_tau(this,tau) result(error_msg)
    class(ty_fluxes_fjx), intent(inout)     :: this
    real(wp), dimension(:,:,:), intent(in)  :: tau
    character(len=128)                      :: error_msg

    error_msg = ""

    !do some check

    this%tau = tau

  end function set_tau
  !--------------------------------------------------------------------------------------------------------------------

  !
  ! Inquiry functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  
  pure function get_njx (this)
    class(ty_fluxes_fjx), intent(in) :: this
    integer :: get_njx

    if (.not. allocated(this%QQQ)) then
      get_njx=0
    else  
      get_njx=size(this%QQQ)
    end if

  end function get_njx        
  !--------------------------------------------------------------------------------------------------------------------

  function get_react_names(this)
    class(ty_fluxes_fjx), intent(in) :: this
    character(len=32), dimension(:), allocatable :: get_react_names

    if (.not. allocated(this%QQQ)) then
      get_react_names=''
    else
      get_react_names=this%QQQ%title
    end if

  end function get_react_names
end module mo_fluxes_fjx

