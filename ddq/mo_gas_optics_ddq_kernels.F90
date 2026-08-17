module mo_gas_optics_ddq_kernels
  use mo_rte_kind,             only: wp, wl
  use mo_gas_optics_constants, only: boltzmann_k, lightspeed, planck_h

  implicit none
  private
  public :: tau_absorption_from_fits, add_tau_rayleigh
  integer, parameter, public :: fax_norder = 2, fax_nterms = 3, xsec_nterms = 3

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Compute absorption optical depth from second-order polynomial approximations
  !    for absorption cross-section
  !
  subroutine tau_absorption_from_fits(ncol, nlay, nnu, ngas, &
                  nus, &
                  play, tlay, dry_num, vmrs, &
                  fax_ngas, fax_num_index, fax_a, fax_b, fax_T0, fax_c, fax_p0, fax_sigma0, fax_S, &
                  xsec_ngas, xsec_num_index, xsec_p, &
                  mtckd_ngas, mtckd_num_index, mtckd_cself, mtckd_cfrgn, mtckd_n, mtckd_T0, mtckd_p0, &
                  tau) bind(C, name="ddq_compute_tau_absorption")
    integer,  intent(in)  :: ncol, nlay, nnu, ngas
    real(wp), intent(in)  :: nus(nnu)
    real(wp), dimension(     ncol, nlay), &
              intent(in)  :: play, tlay, dry_num
! VMRs start at ngas = 0, with element 0 having value 0, to stand in for gases users haven't supplied.
! Column dimension first so the vectorized icol loops below read with unit stride.
    real(wp), intent(in)  :: vmrs(ncol, nlay, 0:ngas)

    ! Functional approximations to cross-sections
    integer,  intent(in) :: fax_ngas
    integer,  intent(in) :: fax_num_index(fax_ngas)
    real(wp), dimension(0:2, fax_ngas, nnu), &
              intent(in)  :: fax_a, fax_b
    real(wp), intent(in)  :: fax_c(0:3, fax_ngas, nnu)
    real(wp), intent(in)  :: fax_sigma0(fax_ngas, nnu)
    real(wp), dimension(fax_ngas) &
                          :: fax_S, fax_T0, fax_p0
    ! Cross-section fits
    integer,  intent(in) :: xsec_ngas
    integer,  intent(in) :: xsec_num_index(xsec_ngas)
    real(wp), intent(in) :: xsec_p(0:3, xsec_ngas, nnu)

    integer,  intent(in) :: mtckd_ngas
    integer,  intent(in) :: mtckd_num_index(mtckd_ngas)
    real(wp), dimension(mtckd_ngas, nnu), &
              intent(in) :: mtckd_cself, mtckd_cfrgn, mtckd_n
    real(wp), intent(in) :: mtckd_T0, mtckd_p0

    real(wp), intent(out) :: tau(ncol, nlay, nnu)
    ! -----------------
    integer  :: igas, icol, ilay, inu
    real(wp) :: vmr
    real(wp) :: x, P_scale, T_scale, delta_T
    real(wp) :: cself, cfrgn, R ! MT_CKD
    ! Per-(igas,inu) coefficients hoisted to scalars
    real(wp) :: c0, c1, c2, xh, a0, a1, a2, b0, b1, b2, sig0
    real(wp) :: q0, q1, q2, q3
    real(wp) :: cs, cf, en, nu_c
    ! Per-layer invariants, computed once per (icol,igas) instead of once
    ! per (icol,igas,inu): all logs and divisions below are nu-independent.
    real(wp) :: acc    (ncol)
    real(wp) :: fax_x  (ncol, fax_ngas), fax_dT(ncol, fax_ngas), fax_w(ncol, fax_ngas)
    real(wp) :: xsec_w (ncol, xsec_ngas)
    real(wp) :: mt_tr  (ncol, mtckd_ngas), mt_logr(ncol, mtckd_ngas), &
                mt_ps  (ncol, mtckd_ngas), mt_pf  (ncol, mtckd_ngas), &
                mt_w   (ncol, mtckd_ngas)
    real(wp) :: inv_2kT(ncol)

    do ilay = 1, nlay
      !
      ! Hoist all nu-independent quantities out of the spectral loop
      !
      do igas = 1, fax_ngas
        do icol = 1, ncol
          vmr = vmrs(icol, ilay, fax_num_index(igas))
          ! Increase pressure to account for self-broadening
          fax_x (icol, igas) = log(play(icol, ilay) * (1 + vmr * fax_S(igas)) / fax_p0(igas))
          fax_dT(icol, igas) = tlay(icol, ilay) - fax_T0(igas)
          fax_w (icol, igas) = vmr * dry_num(icol, ilay)   ! Integrated number density [mol/m**2]
        end do
      end do
      do igas = 1, xsec_ngas
        do icol = 1, ncol
          xsec_w(icol, igas) = vmrs(icol, ilay, xsec_num_index(igas)) * dry_num(icol, ilay)
        end do
      end do
      do igas = 1, mtckd_ngas
        do icol = 1, ncol
          vmr = vmrs(icol, ilay, mtckd_num_index(igas))
          mt_tr  (icol, igas) = mtckd_T0/tlay(icol, ilay)
          mt_logr(icol, igas) = log(mt_tr(icol, igas))
          mt_ps  (icol, igas) = (play(icol, ilay)/mtckd_p0) * vmr
          mt_pf  (icol, igas) = (play(icol, ilay)/mtckd_p0) * (1._wp - vmr)
          mt_w   (icol, igas) = vmr * dry_num(icol, ilay)
        end do
      end do
      do icol = 1, ncol
        inv_2kT(icol) = (planck_h * lightspeed * 100._wp) &
                      / (2._wp * boltzmann_k * tlay(icol, ilay))
      end do

      do inu = 1, nnu
        acc(1:ncol) = 0._wp
        !
        ! Functional approximation to cross-sections
        !
        do igas = 1, fax_ngas
          ! fax_c(3,:,:) is the hinge point x_h
          c0 = fax_c(0, igas, inu); c1 = fax_c(1, igas, inu)
          c2 = fax_c(2, igas, inu); xh = fax_c(3, igas, inu)
          a0 = fax_a(0, igas, inu); a1 = fax_a(1, igas, inu); a2 = fax_a(2, igas, inu)
          b0 = fax_b(0, igas, inu); b1 = fax_b(1, igas, inu); b2 = fax_b(2, igas, inu)
          sig0 = fax_sigma0(igas, inu)
          do icol = 1, ncol
            x       = fax_x (icol, igas)
            delta_T = fax_dT(icol, igas)
            P_scale = c0 + c1 * x + (c2 - c1) * max(x - xh, 0._wp)
            T_scale = (a0 + a1*delta_T + a2*delta_T**2) &
                    / (b0 + b1*delta_T + b2*delta_T**2)
            acc(icol) = acc(icol) &
              + (sig0 * exp(P_scale + T_scale)) & ! cross-section [m**2/mol]
              * fax_w(icol, igas)
          end do
        end do
        !
        ! Cross-sections pressure and temperature dependence following doi:10.1029/2022MS003239
        !
        do igas = 1, xsec_ngas
          q0 = xsec_p(0, igas, inu); q1 = xsec_p(1, igas, inu)
          q2 = xsec_p(2, igas, inu); q3 = xsec_p(3, igas, inu)
          do icol = 1, ncol
            acc(icol) = acc(icol) &
              + (q0 + q1 * tlay(icol, ilay)    &
                    + q2 * tlay(icol, ilay)**2 &
                    + q3 * play(icol, ilay))   &
              * xsec_w(icol, igas)
          end do
        end do
        !
        ! MT_CKD continuum
        !
        do igas = 1, mtckd_ngas
          cs   = mtckd_cself(igas, inu)
          cf   = mtckd_cfrgn(igas, inu)
          en   = 1._wp + mtckd_n(igas, inu)
          nu_c = nus(inu)
          do icol = 1, ncol
            ! (T0/T)**(1+n) written as exp((1+n)*log(T0/T)) with the log hoisted
            cself = exp(en * mt_logr(icol, igas)) * mt_ps(icol, igas) * cs
            cfrgn = mt_tr(icol, igas)             * mt_pf(icol, igas) * cf
            ! nu supplied in kaysers (cm^-1); convert to MKS for tanh argument;
            ! R needs to be in units cm^-1 as cself and cfrgn are in of m^2/molecule cm
            ! R * (cself + cfrgn) is in units of m^2/molecule
            R = nu_c * tanh(nu_c * inv_2kT(icol))
            acc(icol) = acc(icol)   &
              + R * (cself + cfrgn) &
              * mt_w(icol, igas)
          end do
        end do
        do icol = 1, ncol
          tau(icol, ilay, inu) = MAX(0._wp, acc(icol))
        end do
      end do
    end do
  end subroutine tau_absorption_from_fits
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Compute absorption optical depth from second-order polynomial approximations
  !    for absorption cross-section
  !
  subroutine add_tau_rayleigh(ncol, nlay, nnu,  &
                              dry_num,          &
                              rayleigh_xsec,    &
                              tau, ssa) bind(C, name="ddq_add_tau_rayleigh")
    integer,  intent(in)    :: ncol, nlay, nnu
    real(wp), intent(in)    :: dry_num(ncol, nlay)
    real(wp), intent(in)    :: rayleigh_xsec(nnu)
    real(wp), intent(inout) :: tau(ncol, nlay, nnu), ssa(ncol, nlay, nnu)
    ! -----------------
    integer  :: icol, ilay, inu
    real(wp) :: t, t_r

    do inu = 1, nnu
      do ilay = 1, nlay
        do icol = 1, ncol
          t_r = dry_num(icol, ilay) * rayleigh_xsec(inu)
          t = tau(icol, ilay, inu)
          tau(icol, ilay, inu) = t + t_r
          ssa(icol, ilay, inu) = t_r/(t + t_r)
        end do
      end do
    end do
  end subroutine add_tau_rayleigh
  !--------------------------------------------------------------------------------------------------------------------
end module mo_gas_optics_ddq_kernels
