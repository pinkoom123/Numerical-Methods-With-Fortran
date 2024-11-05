MODULE physical_pendulum
  IMPLICIT NONE
CONTAINS

  ! Velocity-Verlet Method Subroutine
  SUBROUTINE verlet_method(x0, t0, v0, dt, tf, tlist, vlist, xlist)
    REAL(8), INTENT(IN) :: x0, t0, v0, dt, tf
    REAL(8), DIMENSION(:), INTENT(OUT) :: tlist, xlist, vlist
    REAL(8) :: t, x, at, v
    INTEGER :: i, N

    ! Calculate N based on tf, t0, and dt
    N = CEILING((tf - t0) / dt)

    ! Initialize conditions
    x = x0
    t = t0
    v = v0
    tlist(1) = t
    xlist(1) = x
    vlist(1) = v

    ! Verlet method loop
    DO i = 1, N
      at = -SIN(x)         ! a(t) = -SIN(x) for this problem
      v = v + at * dt / 2.0  
      x = x + v * dt       ! Update x based on the verlet method
      at = -SIN(x)         ! Recalculate acceleration at new position
      v = v + at * dt / 2.0  ! update velocity
      t = t + dt           ! update time
      tlist(i + 1) = t
      xlist(i + 1) = x
      vlist(i + 1) = v
    END DO
  END SUBROUTINE verlet_method

  ! Euler Method Subroutine
  SUBROUTINE euler_method(x0, t0, v0, dt, tf, tlist, vlist, xlist)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x0, dt, t0, v0, tf
    REAL(8), DIMENSION(:), INTENT(OUT) :: tlist, vlist, xlist
    REAL(8) :: t, x, at, v
    INTEGER :: i, N

    ! Calculate N based on tf, t0, and dt
    N = CEILING((tf - t0) / dt)

    ! Initialize conditions
    x = x0
    t = t0
    v = v0
    tlist(1) = t
    xlist(1) = x
    vlist(1) = v

    ! Euler method loop
    DO i = 1, N
      at = -SIN(x)         ! a(t) = -SIN(x) for this problem
      x = x + v * dt       ! Update x based on the Euler method
      v = v + at * dt      ! update velocity
      t = t + dt           ! update time
      tlist(i + 1) = t
      xlist(i + 1) = x
      vlist(i + 1) = v
    END DO
  END SUBROUTINE euler_method

  ! Subroutine to Write Data to File
  SUBROUTINE write_to_file(filename, tlist, xlist, vlist, N)
    CHARACTER(len=*), INTENT(IN) :: filename
    REAL(8), DIMENSION(:), INTENT(IN) :: tlist, xlist, vlist
    INTEGER, INTENT(IN) :: N
    INTEGER :: i, ios

    OPEN(unit=10, IOSTAT=ios, file=filename, status='replace')
    IF (ios /= 0) THEN
      PRINT*, 'Error opening file'
      STOP
    END IF

    ! Write the header
    WRITE(10, '(A)') 'Time(t)      Position x(t)     Velocity v(t)'
  
    ! Write data to file
    DO i = 1, N
      WRITE(10, '(F10.5, F10.5, F10.5)') tlist(i), xlist(i), vlist(i)
    END DO

    CLOSE(10)
  END SUBROUTINE write_to_file

END MODULE physical_pendulum

! Main Program
PROGRAM Differential_solver
  USE physical_pendulum
  IMPLICIT NONE
  REAL(8) :: x0, t0, dt, v0, tf
  REAL(8), DIMENSION(:), ALLOCATABLE :: xlist_euler, xlist_verlet, vlist_euler, vlist_verlet, tlist
  INTEGER :: N

  ! Initialize conditions
  tf = 30.0
  dt = 0.1
  x0 = 1.0
  t0 = 0.0
  v0 = 0.0
  N = CEILING((tf - t0) / dt)

  ! Allocate arrays
  ALLOCATE(tlist(N), xlist_euler(N), vlist_euler(N), xlist_verlet(N), vlist_verlet(N))

  ! Solve using Euler method
  CALL euler_method(x0, t0, v0, dt, tf, tlist, vlist_euler, xlist_euler)

  ! Solve using Verlet method
  CALL verlet_method(x0, t0, v0, dt, tf, tlist, vlist_verlet, xlist_verlet)

  ! Write results to files
  CALL write_to_file('euler.txt', tlist, xlist_euler, vlist_euler, N)
  CALL write_to_file('verlet.txt', tlist, xlist_verlet, vlist_verlet, N)

  ! Deallocate arrays
  DEALLOCATE(tlist, vlist_euler, xlist_euler, vlist_verlet, xlist_verlet)

  PRINT*, 'Data has been written to euler.txt and verlet.txt'
 END PROGRAM Differential_solver
