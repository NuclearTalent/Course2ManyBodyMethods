!
!
!  TALENT school, course #2 -- GANIL 6-24 July 2015
!
!
!  "Green's function calculations of spectralfucntion for a pairing model"
!
!  This codes used the GF method to solve the toy pairing model with 4 levels 
! and 4 spin-1/2 particles that was introduced during the first week of lectures
! (see http://nucleartalent.github.io/Course2ManyBodyMethods/doc/web/course.html ).
!
!  Calculation canbe done at the different levels of theory by using the ADC(2)
! approach at second order, by adding all order resummations of the two-particle
! and two-hole ladders (the so-called extended-ADC(2) ), and in the full ADC(3).
! approximation.
!
!  It also iterates the energy independent self-energy (i.e. the 'correlated HF'
! diagram with the dressed propagator to achieve self-consistency at the sc0
! level.
!
!
!
!  The ground state is then calculated in FCI to obtain the exact fundamental
! energy and correlation energy, for comparison with the Koltun Sr result. It
! is up to you, as an excercise, to do FCI for the N=3 and N=5 particle systems
! and to extract the exact spectral function to compare with the results from
! the ADC(n) method.
!
!  The coupled cluster equations are also solved to resum the pp and hh ladders.
! These diagrams are a subset of the full CCD method that could be compared to the
! correlations added in the extended-ADC(2).
!
!  The particular pairing Hamiltonian used in this excercise allows us to use
! the following two simplifications:
!
!   - The spin-up and spin-down parts of the Dyson equation completly decouple
!    and are equal because of symmetry. Hence, we need to calculate only one.
!
!   - The particle-hole contribution from the pairing interactions is zero in
!    the ring summation. So, we don't need to worry about rings.
!
!
!
!   You will need to link to the LAPACK library to compile thi sprogram.
!
!
!  C. Barbieri,  Surrey/GANIL,   July 2015.
!
!

module Bases

  implicit none


  integer, parameter :: MaxVectsDim = 16 ! well redundant...

  integer :: N_1b_mdsp
  integer ::   isp_lv(MaxVectsDim),  isp_sp(MaxVectsDim),  isp_oc(MaxVectsDim)
  real*8  ::   xsp_en(MaxVectsDim)


  integer :: N_1b_bas
  integer :: i_1b_bas(MaxVectsDim),i_1b_occ(MaxVectsDim)

  integer :: N_2p1h_bas, N_2h1p_bas
  integer ::  innk_n1(MaxVectsDim), innk_n2(MaxVectsDim), innk_k3(MaxVectsDim)
  integer ::  ikkn_k1(MaxVectsDim), ikkn_k2(MaxVectsDim), ikkn_n3(MaxVectsDim)
  real*8  ::  xnnk_en(MaxVectsDim), xkkn_en(MaxVectsDim)

  integer :: N_2p_bas, N_2h_bas
  integer ::  i_2p_n1(MaxVectsDim), i_2h_k1(MaxVectsDim)
  integer ::  i_2p_n2(MaxVectsDim), i_2h_k2(MaxVectsDim)
  real*8  ::  x_2p_en(MaxVectsDim), x_2h_en(MaxVectsDim)

  !
  !  We will use the following to store the poles and overlap
  ! amplitudes of the propagator
  integer, parameter :: MaxDySols = 20 ! well redundant...
  integer :: N_bk_poles,  N_fw_poles
  real*8 ::  E_bk_poles(MaxVectsDim),            E_fw_poles(MaxVectsDim)
  real*8 :: SF_bk_poles(MaxVectsDim),           SF_fw_poles(MaxVectsDim)
  real*8 ::  Y_bk_poles(MaxVectsDim,MaxDySols),  X_fw_poles(MaxVectsDim,MaxDySols)
  real*8 :: E_Fermi


  logical :: ADC2ex, ADC3, TsqLadd, FullCCD

contains


subroutine Build_bases(xdlt)
  implicit none
real*8, intent(in) :: xdlt


  integer :: i, i1, i2, i3, k

  integer :: n1,n2


  !real*8  :: x1, x2


  N_1b_bas = 0
  N_2p_bas = 0
  N_2h_bas = 0
  N_2p1h_bas = 0
  N_2h1p_bas = 0

  write(6,*)
  write(6,*) ' Building the single particle basis...'

  N_1b_mdsp=8

  i1 = 0
  do i = 1 , N_1b_mdsp

    isp_oc(i) = 0
    if (i <= 4) isp_oc(i) = 1
    isp_lv(i) = (i-1)/2
    isp_sp(i) = -1 + 2*MOD(i,2)
    xsp_en(i) = isp_lv(i) * xdlt

    write(6,*) i, isp_lv(i), isp_SP(i),isp_oc(i),  xsp_en(i)

    if (0 < isp_sp(i)) then
     i1=i1+1
     i_1b_bas(i1) = i
     i_1b_occ(i1) = isp_oc(i)
    end if

  end do

  N_1b_bas = i1


  i1 = 0
  i2 = 0
  do n1 =  1 , N_1b_mdsp
  do n2 = n1+1 , N_1b_mdsp
   do k  = 1 , N_1b_mdsp

    if (isp_lv(n1) /=  isp_lv(n2)) cycle
    if (isp_sp(n1) /= -isp_sp(n2)) cycle
    if (isp_sp(k)  > 0) cycle

    if ( (0 == isp_oc(n1)) .AND. (0 == isp_oc(n2)) .AND. (1 == isp_oc(k)) ) then
      i1=i1+1
      innk_n1(i1) = n1
      innk_n2(i1) = n2
      innk_k3(i1) = k
      xnnk_en(i1) = xsp_en(n1) + xsp_en(n2) - xsp_en(k)
    end if

    if ( (1 == isp_oc(n1)) .AND. (1 == isp_oc(n2)) .AND. (0 == isp_oc(k)) ) then
      i2=i2+1
      ikkn_k1(i2) = n1
      ikkn_k2(i2) = n2
      ikkn_n3(i2) = k
      xkkn_en(i2) = xsp_en(n1) + xsp_en(n2) - xsp_en(k)
    end if

   end do
  end do
  end do

  N_2p1h_bas = i1
  N_2h1p_bas = i2



  write(6,*)
  write(6,*)
  write(6,*) ' Single particle states in the spin-up channel...'
  write(6,*)

  do i = 1 , N_1b_bas
    i1 = i_1b_bas(i)
    write(6,*) i,i_1b_bas(i),i_1b_occ(i),isp_lv(i1),isp_sp(i1),isp_oc(i1),xsp_en(i1)
  end do

  write(6,*)
  write(6,*)
  write(6,*) ' 2p1h basis...'
  write(6,*)
  do i = 1 , N_2p1h_bas
    i1 = innk_n1(i)
    i2 = innk_n2(i)
    i3 = innk_k3(i)
    write(6,100) i,i1,i2,i3,isp_lv(i1),isp_sp(i1),isp_lv(i2),isp_sp(i2),isp_lv(i3),isp_sp(i3),xnnk_en(i)
  end do
100 format(i4,2x,3I3,' - (',i2,',',i2') (',i2,',',i2') (',i2,',',i2')',f10.4)

  write(6,*)
  write(6,*)
  write(6,*) ' 2h1p basis...'
  write(6,*)
  do i = 1 , N_2h1p_bas
    i1 = ikkn_k1(i)
    i2 = ikkn_k2(i)
    i3 = ikkn_n3(i)
    write(6,100) i,i1,i2,i3,isp_lv(i1),isp_sp(i1),isp_lv(i2),isp_sp(i2),isp_lv(i3),isp_sp(i3),xkkn_en(i)
  end do


  write(6,*)




  !
  !  Build the 2qp and 2qh bases. These are used to calculate the
  ! corrections for the ADC(3) coupling and the CCD(ladder) resummations
  !if (ADC3) then

    i1 = 0
    i2 = 0
    do n1 =  1 , N_1b_mdsp
     do n2 = n1+1 , N_1b_mdsp

      if (isp_lv(n1) /=  isp_lv(n2)) cycle
      if (isp_sp(n1) /= -isp_sp(n2)) cycle

      if ( (0 == isp_oc(n1)) .AND. (0 == isp_oc(n2)) ) then
        i1=i1+1
        i_2p_n1(i1) = n1
        i_2p_n2(i1) = n2
        x_2p_en(i1) = xsp_en(n1) + xsp_en(n2)
      end if

      if ( (1 == isp_oc(n1)) .AND. (1 == isp_oc(n2)) ) then
        i2=i2+1
        i_2h_k1(i2) = n1
        i_2h_k2(i2) = n2
        x_2h_en(i2) = xsp_en(n1) + xsp_en(n2)
      end if

     end do
    end do
    N_2p_bas = i1
    N_2h_bas = i2

    write(6,*)
    write(6,*)
    write(6,*) ' 2p basis...'
    write(6,*)
    do i = 1 , N_2p_bas
      i1 = i_2p_n1(i)
      i2 = i_2p_n2(i)
      write(6,105) i,i1,i2,isp_lv(i1),isp_sp(i1),isp_lv(i2),isp_sp(i2),x_2p_en(i)
    end do
105 format(i4,2x,2I3,' - (',i2,',',i2') (',i2,',',i2')',f10.4)

    write(6,*)
    write(6,*)
    write(6,*) ' 2h basis...'
    write(6,*)
    do i = 1 , N_2h_bas
      i1 = i_2h_k1(i)
      i2 = i_2h_k2(i)
      write(6,105) i,i1,i2,isp_lv(i1),isp_sp(i1),isp_lv(i2),isp_sp(i2),x_2h_en(i)
    end do

  !end if




end subroutine Build_bases


end module Bases


program pairing

use Bases
implicit none



integer :: n1,n2


real*8  :: xdlt, xg, x1, x2, x3, x4

integer :: n_sc_itrs

real*8, external :: Exact_energy, Dyson, SolveCCD



  xdlt = 1.d0
  xg = -1.0d0

  write(6,*)
  write(6,*) ' Value of g ? '
  read(5,*) xg


  ADC2ex = .false.
  ADC3   = .false.
  TsqLadd = .false.
  FullCCD = .false.

  n1 = -1
  do while(n1 < 0)
    write(6,*)
    write(6,*) ' type 1 for ADC(2),  2 for ext-ADC(2) or 3 for ADC(3)? '
    read(5,*) n2
    if (1 == n2) then
      n1 = +1
      ADC2ex = .false.
      ADC3   = .false.
    end if
    if (2 == n2) then
      n1 = +1
      ADC2ex = .true.
      ADC3   = .false.
    end if
    if (3 == n2) then
      n1 = +1
      ADC2ex = .true.
      ADC3   = .true.
    end if
  end do

  n1 = -1
  do while(n1 < 0)
    write(6,*)
    write(6,*) ' type 1 for CC(ladder)(2),  2 adding the TVT(ladders) or 3 for Full CCD ? '
    read(5,*) n2
    if (1 == n2) then
      n1 = +1
      TsqLadd = .false.
      FullCCD = .false.
    end if
    if (2 == n2) then
      n1 = +1
      TsqLadd = .true.
      FullCCD = .false.
    end if
    if (3 == n2) then
      n1 = +1
      TsqLadd = .true.
      FullCCD = .true.
    end if
  end do


  !write(6,*)
  !write(6,*) ' Add the ADC(3) coupling corrections (T/F) ? '
  !read(5,*)  ADC3


  call Build_bases(xdlt)

  do n1 = 0, 40
   xg = -1.d0 + n1* 0.05d0

  n_sc_itrs = 0   ! Just do one dyagonalization with the unperturbed HF diagram
  !n_sc_itrs = 10  ! Iterate the sc0 scheme a given number of times
  x1 = Dyson(xdlt,xg,n_sc_itrs)

  x3 = Exact_energy(xdlt,xg)

  x4 = SolveCCD(xdlt,xg,x2)


  !
  ! If one want to plot the results to a file:
  !
  write(8,*) xg,x1+xg-2.0, x3+xg-2.0, x4, x2

  end do

end program pairing


real*8 function Dyson(xdlt,xg,n_sc_itrs)
  use Bases
  implicit none
  real*8,  intent(in) :: xdlt, xg
  integer, intent(in) :: n_sc_itrs

  ! The stuff needed by the LAPACK eigenvalue pakage:
  ! For 'DSYEV':
  integer, parameter :: LDA = 40
  integer, parameter :: LWORK  = 2+(6+2*LDA)*(LDA)  ! For 'DSYEVD' only...
  integer, parameter :: LIWORK = 3+5*LDA              ! For 'DSYEVD' only...
  real*8  :: A(LDA,LDA), B(LDA,LDA)
  real*8  :: W(LDA), WORK(LWORK)
  integer :: IWORK(LIWORK)

  integer :: INFO, LWopt, LIWopt

  integer :: i, i1, i2, i3, i4, j, k, n1,   ih1, ih2

  real*8  :: x1, x2, xA, xSF, xNorm, xKolt, xA_run, xKolt_run

  integer :: Ntot, itr, max_itrs


  real*8, external :: Vpair


  write(6,*)
  write(6,*) ' --- --- ---'
  write(6,*)
  write(6,*) '  Start calculations with the GF method...'
  write(6,*)



  !
  !  Seek for an estimation of the Fermi energy,  to be used later
  ! for separating qps and qhs
  x1 = -1.d20
  x2 = +1.d20
  do  i = 1 , N_1b_bas
    k = i_1b_bas(i)
    if (1 == isp_oc(k)) x1 = MAX(x1,xsp_en(k))
    if (0 == isp_oc(k)) x2 = MIN(x2,xsp_en(k))
  end do
  E_Fermi = (x1+x2)/2.d0
  write(6,'(A,f10.4)') ' The Fermi energy will be taken to be ',E_Fermi

  B(:,:) = 0.d0
  Ntot = N_1b_bas + N_2p1h_bas + N_2h1p_bas
  do i = 1 , N_1b_bas

    i4 = i_1b_bas(i)

    !
    ! This will be added during the sc0 iterations
    !do j = 1 , N_1b_bas
    !  i2 = i_1b_bas(j)
    !  x1 = 0.d0
    !  do k = 1 , N_1b_mdsp
    !    if (1 == isp_oc(k)) x1 = x1 + Vpair(i4,k,i2,k,xg)
    !  end do
    !  B(i,j) = x1
    !end do
    !
    B(i,i) = B(i,i) + xsp_en(i4)

    do j = 1 , N_2p1h_bas
      x1 = Vpair(i4,innk_k3(j),innk_n1(j),innk_n2(j),xg)
      B(i,N_1b_bas+j) = x1
      B(N_1b_bas+j,i) = x1
    end do

    do j = 1 , N_2h1p_bas
      x1 = Vpair(ikkn_k1(j),ikkn_k2(j),i4,ikkn_n3(j),xg)
      B(i,N_1b_bas+N_2p1h_bas+j) = x1
      B(N_1b_bas+N_2p1h_bas+j,i) = x1
    end do

  end do

  do j = 1 , N_2p1h_bas
    B(N_1b_bas+j, N_1b_bas+j) = xnnk_en(j)
  end do
  do j = 1 , N_2h1p_bas
    B(N_1b_bas+N_2p1h_bas+j, N_1b_bas+N_2p1h_bas+j) = xkkn_en(j)
  end do


  if (ADC2ex) then
    !
    !  Add the pp and hh ladders that define the extended-ADC(2) approximation
    ! and are also included in the ADC(3) method.
    !
    !  A comment: diagonalizing these sub matrices separately, would let you
    ! calculate the very same pp or hh ladders that were done during the
    ! coupled cluster CCD exercises.
    !
    n1 = N_1b_bas
    do i = 1 , N_2p1h_bas
     do j = 1 , N_2p1h_bas
      if (innk_k3(i) /= innk_k3(j)) cycle
      x1 = Vpair(innk_n1(i),innk_n2(i),innk_n1(j),innk_n2(j),xg)
      B( n1 + i , n1 + j ) = B( n1 + i , n1 + j ) + x1
     end do
    end do
    !
    n1 = N_1b_bas + N_2p1h_bas
    do i = 1 , N_2h1p_bas
     do j = 1 , N_2h1p_bas
      if (ikkn_n3(i) /= ikkn_n3(j)) cycle
      x1 = Vpair(ikkn_k1(i),ikkn_k2(i),ikkn_k1(j),ikkn_k2(j),xg)
      B( n1 + i , n1 + j ) = B( n1 + i , n1 + j ) - x1
     end do
    end do
    !
  end if

  if (ADC3) then
   !
   !  Add the ADC(3) corrections for the \alpha-2p1h and \alpha-2h1p
   ! coupling temrs.
   !  ONLY the pp/hh ladder corrections are included here because the 
   ! rings do not contribute to this model...(?)
   !
   !
    do i = 1 , N_1b_bas
     i4 = i_1b_bas(i)
     do j = 1 ,  N_2p1h_bas

      x1 = 0.d0
      do k = 1 , N_2h_bas
        ih1 = i_2h_k1(k)
        ih2 = i_2h_k2(k)
        x1 = x1 + Vpair(i4,innk_k3(j),ih1,ih2,xg) * Vpair(ih1,ih2,innk_n1(j),innk_n2(j),xg) /   &
            (x_2h_en(k) - xsp_en(innk_n1(j)) - xsp_en(innk_n2(j)) )
      end do

      B(i,N_1b_bas+j) = B(i,N_1b_bas+j) + x1
      B(N_1b_bas+j,i) = B(N_1b_bas+j,i) + x1

     end do
    end do
    !
    do i = 1 , N_1b_bas
     i4 = i_1b_bas(i)
     do j = 1 , N_2h1p_bas

      x1 = 0.d0
      do k = 1 , N_2p_bas
        ih1 = i_2p_n1(k)
        ih2 = i_2p_n2(k)
        x1 = x1 + Vpair(ikkn_k1(j),ikkn_k2(j),ih1,ih2,xg) * Vpair(ih1,ih2,i4,ikkn_n3(j),xg) /   &
            (xsp_en(ikkn_k1(j)) + xsp_en(ikkn_k2(j)) - x_2p_en(k))
      end do

      B(i,N_1b_bas+N_2p1h_bas+j) = B(i,N_1b_bas+N_2p1h_bas+j) + x1
      B(N_1b_bas+N_2p1h_bas+j,i) = B(N_1b_bas+N_2p1h_bas+j,i) + x1
     end do
    end do
    !
  end if



  max_itrs = MAX(0, n_sc_itrs)
  do itr = 0 , max_itrs

    ! load the part fothe Dyson matrix that will not change with iterations
    A = B

    !
    ! Add the HF(cHF) self energy:
    !
    do i = 1 , N_1b_bas
     i1 = i_1b_bas(i)
     do j = 1 , N_1b_bas
      i2 = i_1b_bas(j)

      if (0 == itr) then
        !
        !  At the first iterations, just use the 1st
        ! order HF diagram
        x1 = 0.d0
        do k = 1 , N_1b_mdsp
          if (1 == isp_oc(k)) x1 = x1 + Vpair(i1,k,i2,k,xg)
        end do
      else
        !
        !  Now we have a propagator and thus we can do
        ! the fully correlated 1
        x1 = 0.d0
        do i3 = 1 , N_1b_bas
         do i4 = 1 , N_1b_bas
          do k  = 1 , N_bk_poles
            ! 
            !  We are using a trick here: in principle we should be calculating also the
            ! spin down part of the spectral function (which we are skipping because it is
            ! exactly the same as the the spin up part--in this model). In fact, both spin
            ! up and down spectral functions can contribute to the static self-energy. For
            ! this particular pairing model, it is only the spin down that contributes to the
            ! spin up HF diagram (because of the interaction we have chosen)!!!
            !  Here we use the calculated spin-up part in whic is stored in Y_bk_poles but
            ! then we need to use the corresponding spin-down s.p. states in for the
            ! interaction matrix elements (hence the i_1b_bas(i...)+1 indices).
            x1 = x1 + Vpair(i1,i_1b_bas(i3)+1,i2,i_1b_bas(i4)+1,xg) * Y_bk_poles(k,i3)*Y_bk_poles(k,i4)
          end do
         end do
        end do
      end if

      A(i,j) = A(i,j) + x1
     end do
    end do


    if (0 == itr) then
      write(6,*)
      write(6,*)
      write(6,*) ' The Dyson matrix is:'
      write(6,*)
      do i = 1 , Ntot
        write(6,110) i, A(i,1:Ntot)
        if (0 == MOD(i,4)) write(6,*)
      end do
110   format(i4,3(3X,4f7.3))
    end if




    ! LAPACK library (w/ DSYEVD):
    INFO = 0;
    LWopt=0
    LIWopt=0
    call dsyevd('Vectors','Upper',Ntot,A,LDA,W,WORK, LWORK,IWORK,LIWORK,INFO);
    if (0 /= INFO) write(6,*) 'Dyson (DSYEV), Wrong value of IFAIL: INFO= ',INFO
    if (0 == INFO) then
      LWopt = MAX(int(WORK(1) + 0.01) , LWopt)
      LIWopt = MAX(IWORK(1) , LIWopt)
    else
      write(6,*) ' ''DEEGV'' gave INFO = ', INFO
      if (INFO < 0) write(6,*) ' The ',-INFO,'-th aggument in line 408 was illegal:'
      write(6,*) ' Program has been aborted since IERR != 0.'
      write(6,*) ' Aborted at the line 410, Dyson.f!'
      write(6,*) '   --> stop.'
      stop
    end if

    if ((max_itrs == itr) .OR. (0 == itr)) then
      write(6,*)
      write(6,*) '              e^+/-      \sum A        SF       \sum KRS          KSR      Corr. En                Check norm'
    end if

    xA = 0.d0
    xKolt = 0.d0
    xA_run = 0.d0
    xKolt_run = 0.d0
    N_bk_poles = 0
    N_fw_poles = 0
    do i = 1 , Ntot

      xNorm = 0.d0; xSF = 0.d0
      x2 = 0.d0
      do j = 1 , Ntot
        x1 = A(j,i)
        xNorm = xNorm + x1*x1
        if (j <= N_1b_bas) then
          xSF = xSF + x1*x1
!          write(6,*) j,xsp_en(i_1b_bas(j))
          x2  = x2  + x1*x1*(xsp_en(i_1b_bas(j)) + W(i))
        end if
      end do
      xSF   = xSF / xNorm
      xNorm = DSQRT(xNorm)

      if (W(i) < E_Fermi) then
        ! quasiholes
         N_bk_poles = N_bk_poles + 1
         E_bk_poles(N_bk_poles) = W(i)
        SF_bk_poles(N_bk_poles) = xSF
         Y_bk_poles(N_bk_poles,1:N_1b_bas) = A(1:N_1b_bas,i) / xNorm
        xKolt = xKolt + x2 / xNorm / 2.d0
        xA    = xA  + xSF
      else
        ! quasiparticles
         N_fw_poles = N_fw_poles + 1
         E_fw_poles(N_fw_poles) = W(i)
        SF_fw_poles(N_fw_poles) = xSF
         X_fw_poles(N_fw_poles,1:N_1b_bas) = A(1:N_1b_bas,i) / xNorm
      end if

      xKolt_run = xKolt_run + x2 / xNorm / 2.d0
      xA_run    = xA_run  + xSF


      if ((max_itrs == itr) .OR. (0 == itr)) then
         write(6,120) i, W(i),  xA_run, xSF, xKolt_run, xKolt_run*2.d0, xKolt_run*2.d0 - 2.d0 + xg, x2, xNorm
         if (0 == MOD(i,6)) write(6,*)
      end if
120   format(i4,3(4X,4f12.5))

    end do

    !
    !  We have caluclate only the spectral function for particles with
    ! spin up  but the spin spin down would be exactly the same and give
    ! another equal contribution. We multiply by 2.d0 to account for this.
    !
    write(6,121) ' Particle number and Koltun SR :  ', 2.d0*xA, 2.d0*xKolt
121 format(A,2f12.5)

    Dyson = 2.d0*xKolt
  end do ! itr loop

end function Dyson




real*8 function Exact_energy(xdlt,xg)
  implicit none
  real*8, intent(in) :: xdlt, xg

  integer :: Ntot

  ! The stuff needed by the LAPACK eigenvalue pakage:
  ! For 'DSYEV':
  integer, parameter :: LDA = 6
  integer, parameter :: LWORK  = 2+(6+2*LDA)*(LDA)  ! For 'DSYEVD' only...
  integer, parameter :: LIWORK = 3+5*LDA              ! For 'DSYEVD' only...
  real*8  :: A(LDA,LDA)
  real*8  :: W(LDA), WORK(LWORK)
  integer :: IWORK(LIWORK)

  integer :: INFO, LWopt, LIWopt, i


  write(6,*)
  write(6,*) ' --- --- ---'
  write(6,*)
  write(6,*) '  Solving the exact ground state energy by diagonalising the '
  write(6,*) ' following hamiltonian:'
  write(6,*)



  A(:,:) = -xg / 2.d0
  Ntot = 6
  do i = 1 , Ntot
    A(i, Ntot+1-i) = 0.d0
    A(i,i) = 2.d0 * i - xg
    if (i > Ntot/2) A(i,i) = xdlt * 2.d0 * (i-1) - xg
  end do

  do i = 1 , Ntot
  write(6,130) i, A(i,1:Ntot)
  if (0 == MOD(i,3)) write(6,*)
  end do
130 format(i4,3(3X,3f7.3))



  ! LAPACK library (w/ DSYEVD):
  INFO = 0;
  LWopt=0
  LIWopt=0
  Ntot = 5
  call dsyevd('Vectors','Upper',Ntot,A,LDA,W,WORK, LWORK,IWORK,LIWORK,INFO);
  if (0 /= INFO) write(6,*) 'Dyson (DSYEV), Wrong value of IFAIL: INFO= ',INFO
  if (0 == INFO) then
    LWopt = MAX(int(WORK(1) + 0.01) , LWopt)
    LIWopt = MAX(IWORK(1) , LIWopt)
  else
    write(6,*) ' ''DEEGV'' gave INFO = ', INFO
    if (INFO < 0) write(6,*) ' The ',-INFO,'-th aggument in line 408 was illegal:'
    write(6,*) ' Program has been aborted since IERR != 0.'
    write(6,*) ' Aborted at the line 410, Dyson.f!'
    write(6,*) '   --> stop.'
    stop
  end if

  Exact_energy = W(1)


  write(6,*)
  write(6,*)
  write(6,142) ' The exact correlation energy is: ', W(1) - 2.d0 + xg
  write(6,*)
  write(6,*)
  write(6,*) 'Eigenvalues:'
  write(6,140) W(1:Ntot)
  write(6,*)
  write(6,*) 'Eivenvectors:'
  do i = 1 , Ntot
    write(6,141) i, A(i,1:Ntot)
  end do
140 format(4X,5(3X,3f10.4))
141 format(i4,5(3X,3f10.4))
142 format(A,f12.5)

  return
end function Exact_energy



real*8 function SolveCCD(xdlt,xg,dEmbpt2)
  use Bases
  implicit none
  real*8, intent(in)  :: xdlt, xg
  real*8, intent(out) :: dEmbpt2

real*8, allocatable :: V2(:,:), T2(:,:), E0(:,:), Vpp(:,:), Vhh(:,:), TT(:,:)

  integer ::  r, q, itr, ia, ib, ii, ij
  real*8  ::  xv, xe, xe_new

  real*8, external :: Vpair


  write(6,*)
  write(6,*) ' --- --- ---'
  write(6,*)
  write(6,*) '  Solving the CC equation for the correlation energy with only '
  write(6,*) ' the pp/hh ladders. This will also yield the MBPT2 result as the:'
  write(6,*) ' 0-th iteration...'
  write(6,*)


  allocate (  T2(N_2p_bas,N_2h_bas) )
  allocate (  V2(N_2p_bas,N_2h_bas) )
  allocate (  E0(N_2p_bas,N_2h_bas) )

  allocate ( Vpp(N_2p_bas,N_2p_bas) )
  allocate ( Vhh(N_2h_bas,N_2h_bas) )

  dEmbpt2 = 0.d0
  do r = 1 , N_2p_bas
    ia = i_2p_n1(r)
    ib = i_2p_n2(r)
    do q = 1 , N_2h_bas
      ii = i_2h_k1(q)
      ij = i_2h_k2(q)

      xv = Vpair(ia,ib,ii,ij,xg)
      xe = x_2h_en(q) - x_2p_en(r)
      V2(r,q) = xv
      E0(r,q) = xe
      T2(r,q) = xv / xe

      dEmbpt2 = dEmbpt2 + xv*xv/xe
    end do
  end do

  do r = 1 , N_2p_bas
   ia = i_2p_n1(r)
   ib = i_2p_n2(r)
   do q = 1 , N_2p_bas
    ii = i_2p_n1(q)
    ij = i_2p_n2(q)
    Vpp(r,q) = Vpair(ia,ib,ii,ij,xg)
   end do
  end do
  do r = 1 , N_2h_bas
   ia = i_2h_k1(r)
   ib = i_2h_k2(r)
   do q = 1 , N_2h_bas
    ii = i_2h_k1(q)
    ij = i_2h_k2(q)
    Vhh(r,q) = Vpair(ia,ib,ii,ij,xg)
   end do
  end do


  write(6,*) '                     dE_old           dE_new'

  xe_new = dEmbpt2
  itr = 0
  do while (ABS(xe - xe_new) > 1.d-7)

    xe = xe_new
    itr = itr + 1

    TT = V2 + MATMUL(T2, Vhh ) + MATMUL(Vpp, T2)

    if (TsqLadd) TT = TT + MATMUL( MATMUL (T2 , TRANSPOSE( V2) ) , T2 )

    if (FullCCD) then
      !
      ! Not coded yet...
      !
    end if

    T2 = TT

    xe_new = 0.d0
    do r = 1 , N_2p_bas
     do q = 1 , N_2h_bas
      T2(r,q) = T2(r,q) / E0(r,q)
      xe_new = xe_new + V2(r,q)*T2(r,q)
     end do
    end do

    write(6,'(A,I5,2f17.7)') 'Itr n, ',itr, xe, xe_new

  end do

  SolveCCD = xe_new

  deallocate (  T2 )
  deallocate (  V2 )
  deallocate (  E0 )
  deallocate ( Vpp )
  deallocate ( Vhh )


  return
end function SolveCCD




real*8 function Vpair(a,b,c,d,xg)
  use Bases
  implicit none
  integer, intent(in) :: a,b,c,d
  real*8,  intent(in) :: xg
  Vpair = 0.d0
  if (isp_lv(a) /=  isp_lv(b)) return
  if (isp_sp(a) /= -isp_sp(b)) return
  if (isp_lv(c) /=  isp_lv(d)) return
  if (isp_sp(c) /= -isp_sp(d)) return
  Vpair = -xg/2.d0
  if (isp_sp(a) /=  isp_sp(c)) Vpair = -Vpair
  return
end function Vpair
