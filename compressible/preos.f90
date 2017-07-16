subroutine getThermo(T)

  implicit none

  ! 0-based indexing to match python
  double precision, intent(inout) :: T(0:qx-1)

!f2py depend(qx, qy) :: r, u, v, p
!f2py depend(qx, qy, nvar) :: q_l, q_r
!f2py intent(in) :: T
!f2py intent(out) :: a, b, R, dadT, d2adT2

  double precision :: dq(0:nvar-1), q(0:nvar-1)
  double precision :: lvec(0:nvar-1,0:nvar-1), rvec(0:nvar-1,0:nvar-1)

  double precision :: MW = 28.0134d-3, Tc = 126.19, pc = 3.3958d06, rhoc = 313.3, omega = 0.03720
  double precision :: c, R, b, 
  double precision :: a(0:)

  c = 0.37464 + 1.54226*omega - 0.26992*omega**2

  R = 8.314d0/MW

  a = 0.457236*(R*Tc)**2 / pc*(1+c*(1-sqrt(T/Tc)))**2

  b = 0.077796*R*Tc/pc
  G = c*np.sqrt(T/Tc) / (1+c*(1-np.sqrt(T/Tc)))
  dadT = -1./T*a*G
  d2adT2 = 0.457236*R**2. / T/2*c*(1+c)*Tc/pc*np.sqrt(Tc/T)

end subroutine states


subroutine getEnergyfromTandRho(T, rho)

  implicit none

  integer, intent(in) :: idir
  integer, intent(in) :: qx, qy, ng
  integer, intent(in) :: nvar, idens, ixmom, iymom, iener
  integer, intent(in) :: lower_solid, upper_solid
  double precision, intent(in) :: gamma

  ! 0-based indexing to match python 
  double precision, intent(inout) :: U_l(0:qx-1,0:qy-1,0:nvar-1)
  double precision, intent(inout) :: U_r(0:qx-1,0:qy-1,0:nvar-1)
  double precision, intent(  out) :: F(0:qx-1,0:qy-1,0:nvar-1)

!f2py depend(qx, qy, nvar) :: U_l, U_r     
!f2py intent(in) :: U_l, U_r
!f2py intent(out) :: e

  coef = [3.531005280E+00,-1.236609870E-04,-5.029994370E-07,2.435306120E-09, -1.408812350E-12,-1.046976280E+03,2.967474680E+00]
   
  a,b,R,dadT,d2adT2 = getThermo(T)
  h_ideal = T*R*(coef[0] + coef[1]*T/2.0 + coef[2]*T**2/3.0 + coef[3]*T**3/4.0 + coef[4] * T**4/5.0+ coef[5] / T)

  v = 1/rho
  K1 = 1.0/(2.0*sqrt(2.0)*b) * log((v+(1.0-sqrt(2))*b)/(v+(1.0+sqrt(2))*b))
  dep = (a - T*dadT)*K1

  e = h_ideal - R*T + dep

end subroutine getEnergyfromTandRho


subroutine getGamma(V)


  implicit none

  ! 0-based indexing to match python 
  double precision, intent(inout) :: U_l(0:qx-1,0:qy-1,0:nvar-1)
  double precision, intent(inout) :: U_r(0:qx-1,0:qy-1,0:nvar-1)
  double precision, intent(  out) :: gammaS

!f2py depend(qx, qy, nvar) :: U_l, U_r     
!f2py intent(in) :: U_l, U_r
!f2py intent(out) :: F

  double precision, parameter :: smallc = 1.e-10
  double precision, parameter :: smallrho = 1.e-10
  double precision, parameter :: smallp = 1.e-10

  rho = V[0]
  p = V[1]
  T = self.getTfromPandRho(p,rho)
  e = self.getEnergyfromTandRho(T,rho)

  sos = self.getSos(V)

  gammaS = sos**2*rho/p 

end subroutine getGamma
   
subroutine consFlux(idir, gamma, idens, ixmom, iymom, iener, nvar, U_state, F)        
  
  integer, intent(in) :: idir
  double precision, intent(in) :: gamma
  integer, intent(in) :: idens, ixmom, iymom, iener, nvar
  double precision, intent(in) :: U_state(0:nvar-1)
  double precision, intent(out) :: F(0:nvar-1)

  double precision :: p, u, v

  u = U_state(ixmom)/U_state(idens)
  v = U_state(iymom)/U_state(idens)

  p = (U_state(iener) - 0.5d0*U_state(idens)*(u*u + v*v))*(gamma - 1.0d0)

  if (idir == 1) then
     F(idens) = U_state(idens)*u
     F(ixmom) = U_state(ixmom)*u + p
     F(iymom) = U_state(iymom)*u
     F(iener) = (U_state(iener) + p)*u
  else
     F(idens) = U_state(idens)*v
     F(ixmom) = U_state(ixmom)*v 
     F(iymom) = U_state(iymom)*v + p
     F(iener) = (U_state(iener) + p)*v
  endif

end subroutine consFlux
  

subroutine artificial_viscosity(qx, qy, ng, dx, dy, &
                                cvisc, u, v, avisco_x, avisco_y)

  implicit none
  integer, intent(in) :: qx, qy, ng
  double precision, intent(in) :: dx, dy
  double precision, intent(in) :: cvisc

  ! 0-based indexing to match python
  double precision, intent(in) :: u(0:qx-1, 0:qy-1)
  double precision, intent(in) :: v(0:qx-1, 0:qy-1)
  double precision, intent(out) :: avisco_x(0:qx-1, 0:qy-1)
  double precision, intent(out) :: avisco_y(0:qx-1, 0:qy-1)

!f2py depend(qx, qy) :: u, v
!f2py depend(qx, qy) :: avisco_x, avisco_y
!f2py intent(in) :: u, v
!f2py intent(out) :: avisco_x, avisco_y

  ! compute the artifical viscosity.  Here, we compute edge-centered
  ! approximations to the divergence of the velocity.  This follows 
  ! directly Colella & Woodward (1984) Eq. 4.5
  !
  ! data locations:
  !
  !   j+3/2--+---------+---------+---------+
  !          |         |         |         |
  !     j+1  +         |         |         |
  !          |         |         |         |
  !   j+1/2--+---------+---------+---------+
  !          |         |         |         |
  !        j +         X         |         |
  !          |         |         |         |
  !   j-1/2--+---------+----Y----+---------+ 
  !          |         |         |         |
  !      j-1 +         |         |         |
  !          |         |         |         | 
  !   j-3/2--+---------+---------+---------+
  !          |    |    |    |    |    |    | 
  !              i-1        i        i+1   
  !        i-3/2     i-1/2     i+1/2     i+3/2 
  !
  ! X is the location of avisco_x(i,j)
  ! Y is the location of avisco_y(i,j)

  integer :: ilo, ihi, jlo, jhi
  integer :: nx, ny
 
  integer :: i, j

  double precision :: divU_x, divU_y

  nx = qx - 2*ng; ny = qy - 2*ng
  ilo = ng; ihi = ng+nx-1; jlo = ng; jhi = ng+ny-1

  do j = jlo-1, jhi+1
     do i = ilo-1, ihi+1

        ! start by computing the divergence on the x-interface.  The
        ! x-difference is simply the difference of the cell-centered
        ! x-velocities on either side of the x-interface.  For the
        ! y-difference, first average the four cells to the node on
        ! each end of the edge, and then difference these to find the
        ! edge centered y difference.
        divU_x = (u(i,j) - u(i-1,j))/dx + &
             0.25d0*(v(i,j+1) + v(i-1,j+1) - v(i,j-1) - v(i-1,j-1))/dy

        avisco_x(i,j) = cvisc*max(-divU_x*dx, 0.0d0)

        ! now the y-interface value
        divU_y = 0.25d0*(u(i+1,j) + u(i+1,j-1) - u(i-1,j) - u(i-1,j-1))/dx + &
             (v(i,j) - v(i,j-1))/dy

        avisco_y(i,j) = cvisc*max(-divU_y*dy, 0.0d0)
        
     enddo
  enddo

end subroutine artificial_viscosity
