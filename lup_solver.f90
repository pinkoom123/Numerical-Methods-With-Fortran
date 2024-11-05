module lup_solver

  implicit none

  private

  real , parameter :: eps = epsilon(1.0)

  public :: solve_lup

  contains
  !
  ! Naive implementation of LUP solver. The order of input matrix a is
  ! columnwise (column first)
  !
  !        a11 a21 a31 ... an1 a12 a22 a32 ... an2 ...
  !
  subroutine solve_lup(a,b,x)
    implicit none
    real , dimension(:,:) , intent(in) :: a
    real , dimension(:) , intent(in) :: b
    real , dimension(:) , intent(out) :: x
    real , dimension(:,:) , allocatable :: lu
    integer , dimension(:) , allocatable :: p
    real , dimension(:) , allocatable :: tmp
    integer , dimension(:) , allocatable :: itmp
    integer :: i , j , n
    n = size(b)
    if ( n /= size(a,1) .or. n /= size(a,2) ) then
      write(*,*) 'CANNOT SOLVE: DIMENSION ERROR'
      stop
    end if
    allocate(lu(n,n),p(n),tmp(n),itmp(n))
    lu = a
    x = b
    ! Decompose : Find LU (in compact form) and P
    do j = 1, n-1
      itmp = maxloc(abs(lu(j:n,j)))+j-1
      if ( itmp(1) /= j ) then
        p(j) = itmp(1)
        tmp(:) = lu(j,:)
        lu(j,:) = lu(p(j),:)
        lu(p(j),:) = tmp(:)
      else
        p(j) = j
      end if
      if ( abs(lu(j,j)) < eps ) then
        write(*,*) 'CANNOT SOLVE : SINGULAR MATRIX'
        stop
      end if
      do i = j+1 , n
        lu(i,j) = lu(i,j)/lu(j,j)
        lu(i,j+1:n) = lu(i,j+1:n) - lu(i,j)*lu(j,j+1:n)
      end do
    end do
    if ( abs(lu(n,n)) < eps ) then
      write(*,*) 'CANNOT SOLVE : SINGULAR MATRIX'
      stop
    end if
    ! Permute
    do j = 1, n-1
      if ( p(j) /= j ) then
        tmp(1) = x(j)
        x(j) = x(p(j))
        x(p(j)) = tmp(1)
      end if
    end do
    ! Forward
    do i = 2 , n
      x(i) = x(i) - dot_product(lu(i,1:i-1),x(1:i-1))
    end do
    ! Backward
    x(n) = x(n)/lu(n,n)
    do i = n-1, 1, -1
      x(i) = (x(i) - dot_product(lu(i,i+1:n),x(i+1:n)))/lu(i,i)
    end do
    deallocate(lu,p,tmp,itmp)
  end subroutine solve_lup

end module lup_solver
