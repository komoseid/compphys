!---------------------------
! Author: Kristine Onsum Moseid
! Date: 04.12.16
! This program solve sthe diffusion equations with
! the jacobi-solver in 2D.
!
! j is space in x direction (Jmax)
! k is space in y direction (Kmax)
! i is time (N) (might be confusing but this is how I like it)
!---------------------------------------

program diffusion_twoD
  implicit none
!-------------- DECLARATIONS ----------------------------
  integer, parameter         :: Jmax = 101
  integer                    :: N
  integer                    :: out
  integer                    :: i,j,k


  real*8,dimension(Jmax,Jmax) :: u_new
  real*8                     :: t_tot
  real*8                     :: ds, alpha, dt ! ds is delta (space), dx=dy
  real*8                     :: epsilon
! ------------------- INITIALIZATIONS ---------------------

  u_new= 0.d0
  u_new(:,Jmax) = 1.d0
  t_tot = 3.d0



  alpha = 10
  epsilon = 0.00000001d0
  ds = 1.d0/(Jmax-1)
  dt = alpha*ds*ds
  N = int(t_tot/dt)
  print*, "N is:", N
  print*, "dt is:", dt

  call jacobi(N,Jmax,alpha,u_new,epsilon)



end program diffusion_twoD


subroutine jacobi(N, Jmax, alpha, u_new, epsilon)
  implicit none
!----------- DECLARATIONS ----------------------------------
  integer                       :: Jmax, N
  integer                       :: j,k,i
  integer                       :: max_i ! max iterations
  integer                       :: out
  integer                       :: watchout

  real*8                        :: alpha, dt, dx, beta
  real*8                        :: epsilon
  real*8                        :: t
  real*8                        :: ds

  real*8, dimension(Jmax, Jmax) :: u
  real*8, dimension(Jmax, Jmax) :: u_new
  real*8, dimension(Jmax, Jmax) :: u_save
  real*8, dimension(Jmax,Jmax)  :: u_plot
! -------------- INITIALIZATIONS ---------------------------
  ds = 1.d0/(Jmax-1)
  dt = alpha*ds*ds


  max_i = 10000

  beta = 1 + 4*alpha

  t = 0.d0

!--------------- LOOPTIME ---------------------------------

  print*, "start looping over jacobis method"

  do i=1,N-1                                     ! TIMELOOP
    u = u_new
    watchout = 0

    do j=2,Jmax-1                                ! SPACELOOP X
      do k=2,Jmax-1                              ! SPACELOOP Y

        u_save = u_new
        u_new(j,k) = (u(j,k) + alpha*( u_save(j,k-1) + u_save(j-1,k) + u_save(j+1,k) + u_save(j,k+1)))/beta

        watchout = watchout + 1
        if (maxval(dabs(u_new - u_star)) < epsilon) then
          print*, "Oops, u_new too close to u_star, Exiting.."
          exit
        if (watchout >= max_i) then
          print*, "over 10000 iterations exceeded. Exiting.."
          exit
        end if

      end do                                    ! END SPACELOOP Y
    end do
                                                ! END SPACELOOP X
    if (i==2000) then
      print*, "Only", N-2000, "iterations left"
    end if
    if (i==2000) then !change this after what you are after
      u_plot(:,:) = u_new(:,:)
    end if

  end do                                        ! END TIMELOOP


!---- WRITE TO SPARATE MATRIX FILES ---------

  open(unit=10, file='jacobi_n2000.txt', ACTION="write", STATUS="replace") !change the name of the file as you wish
  do k=1,Jmax
    write(10, '(1000F14.7)')( real(u_plot(j,k)) ,j=1,Jmax) !my matrix writer for making beautiful plots
  end do
  close(unit=10)

  print*, "Success!"
end subroutine jacobi
