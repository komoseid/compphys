!---------------------------
! Author: Kristine Onsum Moseid
! Date: 04.12.16
! This program includes different methods to solve the diffusion equations with
! three different FDA schemes in 1D
!
! I have used a module tridiag_solver, which is authored by Lars Frogner.
! j is space in x direction (Jmax)
! k is space in y direction (Kmax)
! i is time (N) (might be confusing but this is how I like it)
!---------------------------------------

program diffusion_solver_oneD
  use tridiag_module

  implicit none

!------------ Declarations --------------------


    integer, parameter           :: Jmax = 101
    integer                      :: N

    real*8,dimension(Jmax)       :: u
    real*8                       :: t_tot
    real*8                       :: dx,dt, alpha

! --------- give values needed to run the solvers ---------------------

    t_tot = 1.d0
    dx = 1.d0/(Jmax-1)
    alpha = 0.4
    dt = dx*dx*alpha
    print*, "Hello!"

    N = int(t_tot/dt)

    print*, "N is: ", N

    u = 0.d0

    u(1)= 0.d0
    u(Jmax) = 1.d0

!---------- uncomment what method you would like to choose from ------
    !call forward_euler(Jmax,N,alpha,u)
    !call crank_nicolson(Jmax,N,alpha,u)
    !call backward_euler(Jmax,N,alpha,u)


end program



!-------------------------------------------------------------/
! ------FORWARD EULER FTCS SCHEME ( EXPLICIT ) ---------------/
!---------- only used in 1D situation ------------------------/
!-------------------------------------------------------------/

  subroutine forward_euler(Jmax, N, alpha, u)
    implicit none

!---------------- Declarations ------------------------
    integer                      :: Jmax ! input maximum iterations in x-direction
    integer                      :: N    ! input maximum iterations in time
    integer                      :: i, j ! i for time, j for space in x
!    integer                      :: out  ! helping value to pick out data to write

    real*8                       :: dt, dx, alpha, beta ! alpha = dt/dx**2
    real*8, dimension(Jmax)      :: u ! fluid velocity of interest
    real*8, dimension(Jmax,2)    :: v ! helping velocity for me to use my regular method
    real*8, dimension(N)         :: t ! time
    integer                      :: inew, iold, isave ! helping value for time
    integer,parameter            :: uunit = 10
    character(len=8)             :: frmt
    real*8, dimension(Jmax)      :: v_1, v_2, v_3, v_4, v_5

!--------------- Initializations -----------------------
    dx = 1.d0/(Jmax -1)
    dt = alpha*dx*dx

    beta = 1 - 2*alpha

    if(dt > (dx*dx/2)) then
      print*, "Warning: the forward euler scheme is now unstable"
      continue
    end if


    print*, "initializing"
    v = 0.d0


    !fill up fake velocity v
    do j=1,Jmax
      v(j,1) = u(j)
      v(j,2) = u(j)
    end do

    iold = 1
    inew = 2
    isave = 0

!--------------- looptime ------------------------------
  print*, "start looping Forward Euler"

    do i=1,N
      do j=2,Jmax-1
        v(j,inew) = beta*v(j,iold) + alpha*(v(j+1,iold)  + v(j-1,iold))
      end do

      if ( i == 200) then
        print*, "saving for time:", i*dt
        v_1 = v(:,inew)
      end if
      if ( i == 800) then
        print*, "saving for time:", i*dt
        v_2 = v(:,inew)
      end if
      if ( i == 3000) then
        print*, "saving for time:", i*dt
        v_3 = v(:,inew)
      end if
      if ( i == 10000) then
        print*, "saving for time:", i*dt
        v_4 = v(:,inew)
      end if
      if ( i == 24998) then
        print*, "saving for time:", i*dt
        v_5 = v(:,inew)
      end if


      isave=iold
      iold=inew
      inew=isave
    end do

!---------- WRITE TO FILE -------------------------------
    open(uunit,file='FTCS_unstable.txt')
    do j=1,Jmax
      write (uunit,fmt='(5f30.14)') v_1(j), v_2(j), v_3(j), v_4(j), v_5(j)
    end do
    close(unit=uunit)

  end subroutine forward_euler


!-------------------------------------------------------------/
!------- BACKWARD EULER FTCS SCHEME ( IMPLICIT ) -------------/
!-- only used in 1D situation - together with tridiag-module -/
!-------------------------------------------------------------/
  subroutine backward_euler(Jmax, N, alpha, u)
    use tridiag_module
    implicit none

!----------------------- Declarations ------------------------------------
    integer                      :: Jmax ! input maximum iterations in x-direction
    integer                      :: N    ! input maximum iterations in time
    integer                      :: i, j ! i for time, j for space in x
    integer                      :: out  ! helping value to pick out data to write

    real*8                       :: dt, dx, alpha ! alpha = dt/dx**2
    real*8, dimension(Jmax)      :: u ! fluid velocity of interest
    real*8, dimension(Jmax)      :: u_600, u_1000, u_1500, u_2400

    real*8, dimension(Jmax)      :: diag ! Array representing the diagonal of the probelm-matrix
    real*8, dimension(Jmax-1)    :: diag_up, diag_low ! Arrays representing the upper and lower "lines" along the diagonal
    real*8, dimension(Jmax)      :: u_1, u_2, u_3, u_4

!----------------------initializations ------------------------------------
    dx = 1.d0/ (Jmax-1)
    dt = alpha*dx*dx

    diag = 1 + (2*alpha)
    diag(1) = 1.d0
    diag(Jmax) = 1.d0

    diag_up = -alpha
    diag_up(1) = 0.d0

    diag_low = -alpha
    diag_low(Jmax-1) = 0.d0


    print*, "about to start looping backward euler"
!------------------- Looptime ---------------------------------------------
    do i=1,N                       ! TIMELOOP
      call solve_tridiag(diag_low, diag, diag_up, u)

      if (i==200) then
        u_1 = u
      else if (i==800) then
        u_2 = u
      else if (i==3000) then
        u_3 = u
      else if (i==24998) then
        u_4 = u
      end if

    end do                         ! END TIME LOOP
    print*, "success!"


!------------------ WRITE TO FILE ---------------
    open(unit=10,file='FTCS_back.txt')
    do j=1,Jmax
      write(unit=10,fmt='(4f20.10)') u_1(j), u_2(j), u_3(j), u_4(j)
    end do
    close(unit=10)
  end subroutine backward_euler


!-------------------------------------------------------------/
!--------------- CRANK-NICOLSON (IMPLICIT ) ------------------/
!-- only used in 1D situation - together with tridiag-module -/
!-------------------------------------------------------------/
  subroutine crank_nicolson(Jmax, N, alpha, u)
    use tridiag_module
    implicit none

!------------------------ Declarations -------------------------
    integer                      :: Jmax ! input maximum iterations in x-direction
    integer                      :: N    ! input maximum iterations in time
    integer                      :: i, j ! i for time, j for space in x
    integer                      :: out  ! helping value to pick out data to write

    real*8                       :: dt, dx, alpha, beta , beta_mark ! alpha = dt/dx**2
    real*8, dimension(Jmax)      :: u
    real*8, dimension(Jmax)      :: u_save, tmp
    real*8, dimension(Jmax)      :: u_1, u_2, u_3, u_4


    real*8, dimension(Jmax)      :: diag ! Array representing the diagonal of the probelm-matrix
    real*8, dimension(Jmax-1)    :: diag_up, diag_low ! Arrays representing the upper and lower "lines" along the diagonal

!--------------------- Initializations ------------------------

    dx = 1.d0/(Jmax-1)
    dt = alpha*dx*dx

    beta = 2*(1 + alpha)/alpha

    u_save= 0.d0

    do j=1,Jmax
      u_save(j) = u(j)
    end do

    diag = beta
    diag(1) = 1.d0
    diag(Jmax) = 1.d0

    diag_up = -1.d0
    diag_up(1) = 0.d0

    diag_low = -1.d0
    diag_low(Jmax-1) = 0.d0

    print*, "about to start looping CN"
!--------------------- Looptime --------------------------------

    do i=1,N                                          ! TIMELOOP
      do j=2,Jmax-1                                     ! SPACELOOP
        u_save(j) = u(j-1) + ((2.d0-2.d0*alpha)/alpha)*u(j) + u(j+1)
      end do                                            ! END SPACELOOP


      call solve_tridiag(diag_low,diag,diag_up,u_save)

      if (i==200) then
        u_1 = u_save
      else if (i==800) then
        u_2 = u_save
      else if (i==3000) then
        u_3 = u_save
      else if (i==24999) then
        u_4 = u_save
      end if


      u=u_save
    end do                                             ! END TIMELOOP

    print*, "success!"

!---------------------------- WRITE TO FILE -------------------
    open(unit=10,file='CN.txt')
    do j=1,Jmax
      write(unit=10,fmt='(4f20.10)') u_1(j), u_2(j), u_3(j), u_4(j)
    end do
    close(unit=10)
  end subroutine crank_nicolson
