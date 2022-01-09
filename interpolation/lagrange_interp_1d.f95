subroutine lagrange_basis_1d ( nd, xd, ni, xi, lb ) 

!*****************************************************************************80
!
!! LAGRANGE_BASIS_1D evaluates a 1D Lagrange basis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = rk ) XD(ND), the interpolation nodes.
!
!    Input, integer ( kind = 4 ) NI, the number of evaluation points.
!
!    Input, real ( kind = rk ) XI(NI), the evaluation points.
!
!    Output, real ( kind = rk ) LB(NI,ND), the value, at the I-th point XI, 
!    of the Jth basis function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = rk ) lb(ni,nd)
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)
  
  do i = 1, ni
    do j = 1, nd
      lb(i,j) = product ( ( xi(i) - xd(1:j-1)  ) / ( xd(j) - xd(1:j-1)  ) ) &
              * product ( ( xi(i) - xd(j+1:nd) ) / ( xd(j) - xd(j+1:nd) ) )
    end do
  end do

  return
end

subroutine lagrange_value_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! LAGRANGE_VALUE_1D evaluates the Lagrange interpolant.
!
!  Discussion:
!
!    The Lagrange interpolant L(ND,XD,YD)(X) is the unique polynomial of
!    degree ND-1 which interpolates the points (XD(I),YD(I)) for I = 1
!    to ND.
!
!    The Lagrange interpolant can be constructed from the Lagrange basis
!    polynomials.  Given ND distinct abscissas, XD(1:ND), the I-th Lagrange 
!    basis polynomial LB(ND,XD,I)(X) is defined as the polynomial of degree 
!    ND - 1 which is 1 at  XD(I) and 0 at the ND - 1 other abscissas.
!
!    Given data values YD at each of the abscissas, the value of the
!    Lagrange interpolant may be written as
!
!      L(ND,XD,YD)(X) = sum ( 1 <= I <= ND ) LB(ND,XD,I)(X) * YD(I)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = rk ) XD(ND), the data points.
!
!    Input, real ( kind = rk ) YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = rk ) XI(NI), the interpolation points.
!
!    Output, real ( kind = rk ) YI(NI), the interpolated values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = rk ) lb(ni,nd)
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) yd(nd)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) yi(ni)

  call lagrange_basis_1d ( nd, xd, ni, xi, lb )

  yi = matmul ( lb, yd )

  return
end


subroutine lagrange_value_1d_cmpl ( nd, xd, yd, ni, xi, yi )
    implicit none
  
    integer, parameter :: rk = kind ( 1.0D+00 )
  
    integer ( kind = 4 ) nd
    integer ( kind = 4 ) ni
  
    real ( kind = rk ) lb(ni,nd)
    real ( kind = rk ) xd(nd)
    complex *16 :: yd(nd)
    real ( kind = rk ) xi(ni)
    complex *16 :: yi(ni)
  
    call lagrange_basis_1d ( nd, xd, ni, xi, lb )
    !write(9402,*) lb(i,j)
    yi = matmul ( lb, yd )
  
    return
  end
  