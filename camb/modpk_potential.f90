MODULE potential
  USE modpkparams
  USE internals, ONLY: PI
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pot, getH, getHdot, getEps, dHdphi, dVdphi, d2Vdphi2, getdepsdalpha, powerspectrum, &
       tensorpower, initialphi, norm_u, getHJeps, getHJeta, getHJxi, getHJH, getHJHdot, d3Vdphi3

CONTAINS

  FUNCTION MySech(x)
    REAL*8  :: x,MySech

    IF(ABS(x).GT.40.0) THEN
       MySech=0.0
    ELSE
       MySech=1.0/COSH(x)
    END IF
    RETURN
  END FUNCTION MySech

  FUNCTION pot(phi)
    !
    !     Returns V(phi) given phi
    !
    REAL*8 :: pot,phi,b,c,d,f,meff2,phase
    REAL*8 :: mu, lambda, m_V, finv, p, V0
    REAL*8 :: eps, eta, xi, A1, A2, A3
    REAL*8 :: c_V, d_V, phi_step
    REAL*8 :: fac
    REAL*8 :: C_0, phi_0, q, s, delta, amp1, amp2, ang, ang2
    integer :: ii
    
    ! Parameters for inflection point inflation
    REAL*8 :: lambda_ii, M_ii, delta_ii      

    select case(potential_choice)
    case(1)
       m_V = 10.d0**(0.5d0*vparams(1))
       pot = 0.5*m_V*m_V*phi*phi
    case(2)
       lambda = 10.d0**vparams(1)
       finv = 1.d0/(10.d0**vparams(2))
       pot = lambda**4*(1.d0+cos(finv*phi))
    case(3)
       lambda = 10.d0**vparams(1)
       pot = 0.25*lambda*phi**4
    case(4)
       lambda = 10.d0**vparams(1)
       pot = lambda*phi
    case(5)
       lambda = 10.d0**vparams(1)
       pot = lambda*1.5d0*phi**(2./3.)
    case(6)
       lambda = 10.d0**vparams(1)
       mu = 10.d0**vparams(2)
       pot = lambda**4 - mu*phi**4/4.d0
    case(7)
       eps=vparams(1)
       eta=vparams(2)
       xi=vparams(3)
       A1=sqrt(eps/2.)
       A2=eta/4.
       A3=xi/(12.*sqrt(eps*2))
       pot = 3*H0**2*((1.+A1*phi+A2*phi**2+A3*phi**3)**2-2./3.*(A1+2.*A2*phi+3.*A3*phi**2)**2)
    case(8)
       eps=10.d0**vparams(1)
       eta=vparams(2)
       xi=vparams(3)
       A1=sqrt(eps/2.)
       A2=eta/4.
       A3=xi/(12.*sqrt(eps*2))
       pot = 3*H0**2*((1.+A1*phi+A2*phi**2+A3*phi**3)**2-2./3.*(A1+2.*A2*phi+3.*A3*phi**2)**2)
    case(9)
       mu = 10d0**vparams(1)
       f = 10d0**vparams(2)
       b = vparams(3)
       c = 1d0 - b*f
       pot = mu**3 *((1d0 + phi**2)**.5d0 - b*f * cos(phi/f) -c)
    case(10)
       phase = vparams(4) 
       mu = 10d0**vparams(1)
       f = 10d0**vparams(2)
       b = vparams(3)
       c = - b*f*cos(phase) 
       pot = mu**3 *(phi- b*f * cos(phi/f+phase) -c)
    case(11)
       m_V = 10.d0**(0.5d0*vparams(1))
       c_V = vparams(2)
       d_V = vparams(3)
       phi_step = vparams(4)
       pot=0.5*m_V*m_V*phi*phi*(1.+c_V*TANH((phi-phi_step)/d_V))
    case(12)
       lambda = 10.d0**vparams(1)
       fac=sqrt(2.d0/3.d0)
       pot = lambda**4*(1.d0-exp(-fac*phi))**2.d0
    case(13)
       p=4./3.
       phi_0 = sqrt(2.*p*57.5)
       mu = 10d0**vparams(1)
       f = 10d0**vparams(2)       
       b = vparams(3)
       phase = vparams(4) 
       s = vparams(5)
       q = vparams(6)
       C_0 = vparams(7)
       delta=1d-4
       amp1 = mu**(4.-p)
       amp2 = 0.5d0*b* phi_0**(p-1) * f * p /(s+1.)
       ang = phase + (phi_0/f)*((phi+delta)/phi_0)**(s+1.)
       c = 10.d0
       d = 0.5d0
       ang2 = (phi-c)/d
       pot = amp1 *(phi**p +  amp2*(1.+ tanh(ang2))* cos(ang))
 !MODIFICATION BY LILLIAN   
    case(14)
       !Inflection point inflation
       !phi
       !log priors

       M_ii = vparams(1)
       !Delta_ii = (10.d0**vparams(2))+1.d0 !Delta>1
       Delta_ii = 1.d0 - 10.d0**vparams(2) !Delta<1
       lambda_ii = 10.d0**vparams(3)

       pot = lambda_ii * (&
          + M_ii**2 /2.0d0 * phi**2 &
          - 2.0d0 / 3.0d0 * M_ii * Delta_ii * phi**3 &
          + 1.0d0 / 4.0d0 * phi**4 &
          )
          
    case(15)
       !Inflection point inflation
       !phi
       !uniform Delta prior

       M_ii = vparams(1)
       !Delta_ii = 1.d0 + vparams(2) !Delta>1
       Delta_ii = 1.d0 - vparams(2) !Delta<1
       lambda_ii = 10.d0**vparams(3)

       pot = lambda_ii * (&
          + M_ii**2 /2.0d0 * phi**2 &
          - 2.0d0 / 3.0d0 * M_ii * Delta_ii * phi**3 &
          + 1.0d0 / 4.0d0 * phi**4 &
          )
    
    case(16)
        !Inflection point inflation
        !psi
        !log priors

       M_ii = 10.d0**vparams(1)
       delta_ii = vparams(2)
       lambda_ii = 10.d0**vparams(3)

       pot = lambda_ii/12.0d0 * (&
          M_ii**4 * (1 + 6*delta_ii) + &
          M_ii**3 * 12 * delta_ii * phi + &
          M_ii**1 * 2 * (2 - 3 * delta_ii) * phi**3 + &
          3 * phi**4 &
          )
          
    case(17)
        !Inflection point inflation
        !psi

       M_ii = vparams(1)
       delta_ii = vparams(2)
       lambda_ii = 10.d0**vparams(3)

       pot = lambda_ii/12.0d0 * (&
          M_ii**4 * (1 + 6*delta_ii) + &
          M_ii**3 * 12 * delta_ii * phi + &
          M_ii**1 * 2 * (2 - 3 * delta_ii) * phi**3 + &
          3 * phi**4 &
          )
  !END MODIFICATION BY LILLIAN
    case default
       write(*,*) 'MODPK: Need to set pot(phi) in modpk_potential.f90 for potential_choice =',potential_choice
       STOP
    end select

    RETURN
  END FUNCTION pot


  FUNCTION dVdphi(phi)
    !
    !     Returns dV/dPhi given phi
    !
    REAL*8 :: phi,dVdPhi,dphi,phiplus, b, c, d,f,phase
    REAL*8 :: m_V, lambda, finv, mu, V0, p , meff2, dmeff2dphi
    REAL*8 :: eps, eta, xi, A1, A2, A3
    REAL*8 :: c_V, d_V, phi_step
    REAL*8 :: fac
    REAL*8 :: C_0, phi_0, q, s, delta, amp1, amp2, ang, ang2, pop0
    integer :: ii
    
    ! For inflection point inflation
    REAL*8 :: p1, p2, p3, &
       term0, term1, &
       term2, term3, term4

    ! Parameters for inflection point inflation
    REAL*8 :: lambda_ii, M_ii, delta_ii

    if (vnderivs) then
          phiplus = phi + 0.5*phi*findiffdphi**(1./3.)
          dphi = phiplus - phi
          dVdphi = (pot(phi+dphi)-pot(phi-dphi))/(2.*dphi)
          if (dVdphi.eq.0.d0) then
             write(*,*) 'MODPK: dVdphi=0, possibly signaling a problem with accuracy of numerical derivatives.'
             write(*,*) 'MODPK: Try using vnderivs=F if possible.'
             STOP
          end if
    else
       select case(potential_choice)
       case(1)
          m_V = 10.d0**(0.5d0*vparams(1))
          dVdphi = m_V*m_V*phi
       case(2)
          lambda = 10.d0**vparams(1)
          finv = 1.d0/(10.d0**vparams(2))
          dVdphi = -lambda**4*finv*sin(finv*phi)
       case(3)
          lambda = 10.d0**vparams(1)
          dVdphi = lambda*phi**3
       case(4)
          lambda = 10.d0**vparams(1)
          dVdphi = lambda
       case(5)
          lambda = 10.d0**vparams(1)
          dVdphi = lambda*phi**(-1./3.)
       case(6)
          mu = 10.d0**vparams(2)
          dVdphi = -mu*phi**3
       case(7)
          eps=vparams(1)
          eta=vparams(2)
          xi=vparams(3)
          A1=sqrt(eps/2.)
          A2=eta/4.
          A3=xi/(12.*sqrt(eps*2))
          dVdphi = 3*H0**2*(2*(1.+A1*phi+A2*phi**2+A3*phi**3)*(A1+2*A2*phi+3*A3*phi**2) &
               -4./3*(A1+2.*A2*phi+3.*A3*phi**2)*(2.*A2+6.*A3*phi))
       case(8)
          eps=10.d0**vparams(1)
          eta=vparams(2)
          xi=vparams(3)
          A1=sqrt(eps/2.)
          A2=eta/4.
          A3=xi/(12.*sqrt(eps*2))
          dVdphi = 3*H0**2*(2*(1.+A1*phi+A2*phi**2+A3*phi**3)*(A1+2*A2*phi+3*A3*phi**2) &
               -4./3*(A1+2.*A2*phi+3.*A3*phi**2)*(2.*A2+6.*A3*phi))

       case(9)
          mu = 10d0**vparams(1)
          f = 10d0**vparams(2)
          b = vparams(3)
          dVdphi = mu**3 *(phi*(1d0 + phi**2)**(-.5d0) + b * sin(phi/f))

       case(10)
          phase = vparams(4) 
          mu = 10d0**vparams(1)
          f = 10d0**vparams(2)
          b = vparams(3)
          dVdphi = mu**3 *(1d0 + b * sin(phi/f+phase))

       case(11)
          m_V = 10.d0**(0.5d0*vparams(1))
          c_V = vparams(2)
          d_V = vparams(3)
          phi_step = vparams(4)
          dVdphi=m_V*m_V*(phi*(1.+c_V*TANH((phi-phi_step)/d_V)) &
               &     + 0.5*phi*phi*c_V/d_V*(MySech((phi-phi_step)/d_V)**2.))

       case(12)
          lambda = 10.d0**vparams(1)
          fac=sqrt(2.d0/3.d0)
          dVdphi=2.d0*exp(-fac*phi)*(1.d0 - exp(-(fac*phi)))*fac*lambda**4
                    
       case(13)
          p=4./3.
          phi_0 = sqrt(2.*p*57.5)
          mu = 10d0**vparams(1)
          f = 10d0**vparams(2)       
          b = vparams(3)
          phase = vparams(4) 
          s = vparams(5)
          q = vparams(6)
          C_0 = vparams(7)
          delta=1d-4
          amp1 = mu**(4.-p)
          amp2 = 0.5d0*b* phi_0**(p-1) * f * p /(s+1.)
          ang = phase + (phi_0/f)*((phi+delta)/phi_0)**(s+1.)
          c = 10.d0
          d = 0.5d0
          ang2 = (phi-c)/d
          pop0 = ((phi+delta)/phi_0)
          dVdphi = amp1 *(p*(phi**(p-1.)) +  & 
               & amp2* (cos(ang)/d*MySech(ang2)**2 - (1.+s)/f*(pop0**s)*sin(ang)*(1.d0+tanh(ang2))))
       
       case(14)
          !Inflection point inflation
          !phi
          !log priors

          M_ii = vparams(1)
          !Delta_ii = (10.d0**vparams(2))+1.d0 !Delta>1
          Delta_ii = 1.d0 - 10.d0**vparams(2) !Delta<1
          lambda_ii = 10.d0**vparams(3)


          dVdphi = lambda_ii * ( &
             + M_ii**2 * phi &
             - 2.0d0 * M_ii * Delta_ii * phi**2 &
             + phi**3 &
             )
             
       case(15)
          !Inflection point inflation
          !phi
          !uniform priors

          M_ii = vparams(1)
          !Delta_ii = 1.d0 + vparams(2) !Delta>1
          Delta_ii = 1.d0 - vparams(2) !Delta<1
          lambda_ii = 10.d0**vparams(3)



          dVdphi = lambda_ii * ( &
             + M_ii**2 * phi &
             - 2.0d0 * M_ii * Delta_ii * phi**2 &
             + phi**3 &
             )

       case(16)
          !Inflection point inflation
          !Log priors

           M_ii = 10.d0**vparams(1)
           delta_ii = vparams(2)
           lambda_ii = 10.d0**vparams(3)

          dVdphi = lambda_ii/12.0d0 * (&
             0 + &
             M_ii**3 * 12 * delta_ii + &
             M_ii**1 * 6 * (2 - 3 * delta_ii) * phi**2 + &
             12 * phi**3 &
             )
             
       case(17)
          !Inflection point inflation
          !Uniform priors

          M_ii = vparams(1)
          delta_ii = vparams(2)
          lambda_ii = 10.d0**vparams(3)

          dVdphi = lambda_ii/12.0d0 * (&
             0 + &
             M_ii**3 * 12 * delta_ii + &
             M_ii**1 * 6 * (2 - 3 * delta_ii) * phi**2 + &
             12 * phi**3 &
             )

       case default
          write(*,*) 'MODPK: Need to set dVdphi in modpk_potential.f90 or use numerical derivatives (vnderivs=T)'
          STOP
       end select

    end if

    RETURN
  END FUNCTION dVdphi


  FUNCTION d2Vdphi2(phi)
    !
    !     Returns d^2V/dPhi^2 given phi
    !
    REAL*8 :: phi,d2VdPhi2,dphi,phiplus, b,c,d,f,phase,meff2,d2meff2dphi2,dmeff2dphi
    REAL*8 :: m_V, lambda, finv, mu, V0, p
    REAL*8 :: eps, eta, xi, A1, A2, A3, part1, part2
    REAL*8 :: c_V, d_V, phi_step
    REAL*8 :: fac
    REAL*8 :: C_0, phi_0, q, s, delta, amp1, amp2, ang, pop0, ang2
    integer :: ii
    
    ! For inflection point inflation
    REAL*8 :: p1, p2, p3, &
        term0, term1, &
        term2, term3, term4

    ! Parameters for inflection point inflation
    REAL*8 :: lambda_ii, M_ii, delta_ii

    if (vnderivs) then
       phiplus = phi + 0.2*phi*findiffdphi**(1./4.)
       dphi = phiplus - phi
       d2Vdphi2 = (pot(phi+2.*dphi)+pot(phi-2.*dphi)-2.*pot(phi))/(4.*dphi*dphi)
    else
       select case(potential_choice)
       case(1)
          m_V = 10.d0**(0.5d0*vparams(1))
          d2Vdphi2 = m_V*m_V
       case(2)
          lambda = 10.d0**vparams(1)
          finv = 1.d0/(10.d0**vparams(2))
          d2Vdphi2 = -lambda**4*finv*finv*cos(finv*phi)
       case(3)
          lambda = 10.d0**vparams(1)
          d2Vdphi2 = 3.d0*lambda*phi*phi
       case(4)
          d2Vdphi2 = 0.d0
       case(5)
          lambda = 10.d0**vparams(1)
          d2Vdphi2 = -lambda/3.d0*phi**(-4./3.)
       case(6)
          mu = 10.d0**vparams(2)
          d2Vdphi2 = -3.d0*mu*phi*phi
       case(7)
          eps=vparams(1)
          eta=vparams(2)
          xi=vparams(3)
          A1=sqrt(eps/2.)
          A2=eta/4.
          A3=xi/(12.*sqrt(eps*2))
          part1 = 3*H0**2*(2*(A1+2*A2*phi+3.*A3*phi**2)*(A1+2*A2*phi+3*A3*phi**2)-4./3*(2.*A2+6.*A3*phi)*(2.*A2+6.*A3*phi))
          part2 = 3*H0**2*(2*(1.+A1*phi+A2*phi**2+A3*phi**3)*(2*A2+6*A3*phi)-4./3*(A1+2.*A2*phi+3.*A3*phi**2)*(6.*A3))
          d2Vdphi2 =part1 + part2
       case(8)
          eps=10.d0**vparams(1)
          eta=vparams(2)
          xi=vparams(3)
          A1=sqrt(eps/2.)
          A2=eta/4.
          A3=xi/(12.*sqrt(eps*2))
          part1 = 3*H0**2*(2*(A1+2*A2*phi+3.*A3*phi**2)*(A1+2*A2*phi+3*A3*phi**2)-4./3*(2.*A2+6.*A3*phi)*(2.*A2+6.*A3*phi))
          part2 = 3*H0**2*(2*(1.+A1*phi+A2*phi**2+A3*phi**3)*(2*A2+6*A3*phi)-4./3*(A1+2.*A2*phi+3.*A3*phi**2)*(6.*A3))
          d2Vdphi2 =part1 + part2

       case(9)
          mu = 10d0**vparams(1)
          f = 10d0**vparams(2)
          b = vparams(3)
          d2Vdphi2 = mu**3 *((1d0 + phi**2)**(-.5d0) -  phi**2 *(1d0 + phi**2)**(-1.5d0) + b/f * cos(phi/f))


       case(10)
          phase = vparams(4)
          mu = 10d0**vparams(1)
          f = 10d0**vparams(2)
          b = vparams(3)
          d2Vdphi2 = mu**3 *( b/f * cos(phi/f+phase))

       case(11)
          m_V = 10.d0**(0.5d0*vparams(1))
          c_V = vparams(2)
          d_V = vparams(3)
          phi_step = vparams(4)
          d2Vdphi2=m_V*m_V*((1.+c_V*TANH((phi-phi_step)/d_V)) &
               &     + (MySech((phi-phi_step)/d_V)**2.)*c_V* &
               &     (2./d_V*phi - phi*phi/d_V/d_V*TANH((phi-phi_step)/d_V)))

       case(12)
          lambda = 10.d0**vparams(1)
          fac=sqrt(2.d0/3.d0)
          d2Vdphi2=2.d0*fac*fac*exp(-fac*phi)*lambda**4*(2.d0*exp(-(fac*phi))-1.d0)

       case(13)
          p=4./3.
          phi_0 = sqrt(2.*p*57.5)
          mu = 10d0**vparams(1)
          f = 10d0**vparams(2)       
          b = vparams(3)
          phase = vparams(4) 
          s = vparams(5)
          q = vparams(6)
          C_0 = vparams(7)
          delta=1d-4
          amp1 = mu**(4.-p)
          amp2 = 0.5d0*b* phi_0**(p-1) * f * p /(s+1.)
          ang = phase + (phi_0/f)*((phi+delta)/phi_0)**(s+1.)
          c = 10.d0
          d = 0.5d0
          ang2 = (phi-c)/d
          pop0 = ((phi+delta)/phi_0)
          d2Vdphi2 = amp1 *(p*(p-1)*(phi**(p-2.)) -  amp2 *  &
               & ((2.d0*(1.+s)*pop0**s/d/f*sin(ang) + 2.d0*tanh(ang2)*cos(ang)/d**2*cos(ang))*MySech(ang2)**2 &
               & + ( ((1.+s)/f)**2 * (pop0**(2.*s))*cos(ang) + s*(1.+s)/f/phi_0*pop0**(s-1.)*sin(ang))*(1.+tanh(ang2))))
       
       case(14)
          !Inflection point inflation
          !phi
          !log priors

          M_ii = vparams(1)
          !Delta_ii = (10.d0**vparams(2))+1.d0 !Delta>1
          Delta_ii = 1.d0 - 10.d0**vparams(2) !Delta<1
          lambda_ii = 10.d0**vparams(3)

          d2Vdphi2 = lambda_ii * ( &
             + M_ii**2 &
             - 4.0d0 * M_ii * Delta_ii * phi &
             + 3.0d0 *  phi**2 &
          )
          
       case(15)
          !Inflection point inflation
          !phi
          !uniform priors

          M_ii = vparams(1)
          !Delta_ii = 1.d0 + vparams(2) !Delta>1
          Delta_ii = 1.d0 - vparams(2) !Delta<1
          lambda_ii = 10.d0**vparams(3)

          d2Vdphi2 = lambda_ii * ( &
             + M_ii**2 &
             - 4.0d0 * M_ii * Delta_ii * phi &
             + 3.0d0 *  phi**2 &
          )

       case(16)
          !Inflection point inflation
          !Log priors

          M_ii = 10.d0**vparams(1)
          delta_ii = vparams(2)
          lambda_ii = 10.d0**vparams(3)

          ! The second derivative is initially 0; only need to set the values on
          ! the diagonal.
          d2Vdphi2 = lambda_ii/12.0d0 * ( &
             0 + &
             0 + &
             M_ii**1 * 12 * (2 - 3 * delta_ii) * phi**1 + &
             36 * phi**2 &
             )
             
       case(17)
          !Inflection point inflation
          !Uniform priors

          M_ii = vparams(1)
          delta_ii = vparams(2)
          lambda_ii = 10.d0**vparams(3)

          ! The second derivative is initially 0; only need to set the values on
          ! the diagonal.
          d2Vdphi2 = lambda_ii/12.0d0 * ( &
             0 + &
             0 + &
             M_ii**1 * 12 * (2 - 3 * delta_ii) * phi**1 + &
             36 * phi**2 &
             )

       case default
          write(*,*) 'MODPK: Need to set d2Vdphi2 in modpk_potential.f90 or use numerical derivatives (vnderivs=T)'
          STOP
       end select

    end if

    RETURN
  END FUNCTION d2Vdphi2

  FUNCTION dHdphi(phi)
    !
    ! returns dH/dphi given phi, only used for slow-roll reconstruction
    !
    
    REAL*8 :: phi,dHdphi
    REAL*8 :: eps, eta, xi, A1, A2, A3

    select case(potential_choice)
    case(7)
       eps=vparams(1)
       eta=vparams(2)
       xi=vparams(3)
       A1=sqrt(eps/2.)
       A2=eta/4.
       A3=xi/(12.*sqrt(eps*2))
       dHdphi = H0*(A1 + 2*A2*phi + 3*A3*phi**2)
    case(8)
       eps=10.d0**vparams(1)
       eta=vparams(2)
       xi=vparams(3)
       A1=sqrt(eps/2.)
       A2=eta/4.
       A3=xi/(12.*sqrt(eps*2))
       dHdphi = H0*(A1 + 2*A2*phi + 3*A3*phi**2)
    case default
       write(*,*) 'MODPK: ERROR, dHdphi called for a case not corresponding to slow-roll reconstruction'
       STOP
    end select

    RETURN
  END FUNCTION dHdphi

  FUNCTION d3Vdphi3(phi) result(third_deriv)
    REAL*8 :: phi
    REAL*8 :: third_deriv
    REAL*8 :: m2_V
    REAL*8 :: lambda
    integer :: ii

    REAL*8 :: p_exp

    REAL*8 :: location_phi, step_size, step_slope

    ! For inflection point inflation
    REAL*8 :: p1, p2, p3, &
        term0, term1, &
        term2, term3, term4

    ! Parameters for inflection point inflation
    REAL*8 :: lambda_ii, M_ii, delta_ii

    third_deriv = 0d0

    select case(potential_choice)
    case(14)
       !Inflection point inflation
       !phi
       !log priors

       M_ii = vparams(1)
       !Delta_ii = (10.d0**vparams(2))+1.d0 !Delta>1
       Delta_ii = 1.d0 - 10.d0**vparams(2) !Delta<1
       lambda_ii = 10.d0**vparams(3)

       third_deriv = lambda_ii * ( &
           + 0.0d0 * M_ii**2 &
           - 4.0d0 * M_ii * Delta_ii &
           + 6.0d0 * phi**1 &
           )
    
    case(15)
       !Inflection point inflation
       !phi
       !uniform priors

       M_ii = vparams(1)
       !Delta_ii = 1.d0 + vparams(2) !Delta>1
       Delta_ii = 1.d0 - vparams(2) !Delta<1
       lambda_ii = 10.d0**vparams(3)

       third_deriv = lambda_ii * ( &
           + 0.0d0 * M_ii**2 &
           - 4.0d0 * M_ii * Delta_ii &
           + 6.0d0 * phi**1 &
           )

    case(16)
       !Inflection point inflation
       !Log priors

       M_ii = 10.d0**vparams(1)
       delta_ii = vparams(2)
       lambda_ii = 10.d0**vparams(3)

       ! The third derivative is initially 0; only need to set the values on
       ! the diagonal.
       third_deriv = lambda_ii/12.0d0 * ( &
           0 + &
           0 + &
           M_ii**1 * 12 * (2 - 3 * delta_ii) + &
           72 * phi**1 &
           )
    
    case(17)
       !Inflection point inflation
       !Uniform priors

       M_ii = vparams(1)
       delta_ii = vparams(2)
       lambda_ii = 10.d0**vparams(3)

       ! The third derivative is initially 0; only need to set the values on
       ! the diagonal.
       third_deriv = lambda_ii/12.0d0 * ( &
           0 + &
           0 + &
           M_ii**1 * 12 * (2 - 3 * delta_ii) + &
           72 * phi**1 &
           )

    case default

    write(*,*) "MODECODE: potential_choice =", potential_choice
    write(*,*) "MODPK: Need to set third derivative for this potential choice."

    end select

  END FUNCTION d3Vdphi3

  FUNCTION initialphi(phi0)
    !
    !     Sets initial value of phi (depending on potential, may use
    !     either the user-specified value phi0 or something derived 
    !     from the potential parameters)
    !
    REAL*8 :: initialphi, phi0
    REAL*8 :: phii, Ninit, finv, lambda, mu, phesq, m_V
    REAL*8 :: x1, x2

    Ninit = 70.d0

    select case(potential_choice)
    case(1)
       phii = 2.d0*sqrt(Ninit+0.5d0)
    case(2)
       finv = 1.d0/(10.d0**vparams(2))
       phii = 2.d0/finv*asin(exp(-0.5d0*Ninit*finv*finv)/ &
            sqrt(1.d0+0.5d0*finv*finv))
    case(3)
       phii = sqrt(8.d0*(Ninit+1.d0))
    case(4)
       phii = sqrt(2.d0*Ninit+0.5d0)
    case(5)
       phii = sqrt(4.d0/3.d0*Ninit+2.d0/9.d0)
    case(6)
       lambda = 10.d0**vparams(1)
       mu = 10.d0**vparams(2)
       x1 = lambda**4/mu
       phesq = ((sqrt(2.d0)*x1)**(-4./3.)+1.d0/(4.d0*x1))**(-0.5)
       if (vparams(1)<-3.d0) then
          phii = sqrt(phesq/(2.d0**1.5*Ninit/sqrt(phesq)+1.d0))
       else
          x2 = 4.d0*Ninit + 2.d0*x1/phesq + 0.5d0*phesq
          phii = sqrt(x2)*sqrt(1.d0-sqrt(1.d0-4.d0*x1/x2/x2))
       end if
    case(7)
       phii = 0.
    case(8)
       phii = 0.
    case(9)
       phii = sqrt(2.d0*Ninit+0.5d0)

    case(10)
       phii = sqrt(2.d0*Ninit+0.5d0)

    case(11)
       phii = sqrt(2.d0*Ninit+0.5d0)
       
    case(12)
       phii = 5.8d0

    case(13)
       phii = sqrt(2.d0*Ninit+0.5d0)
    
    ! case(16)
    !     PRINT*, "MODECODE: Setting phi_init algorithmically"
    !     phii = phi0
    !     phii = 2.5e1_dp
    !     ! phii = .1d0 * vparams(2,1)**1.5
    !     ! PRINT*, phi0
    !     ! PRINT*, phii
    ! case(17)
    !     PRINT*, "MODECODE: Setting phi_init algorithmically"
    !     phii = phi0
    !     phii = 2.5e1_dp
    !     ! phii = .1d0 * vparams(2,1)**1.5
    !     ! PRINT*, phi0
    !     ! PRINT*, phii

    case default
       phii = phi0
    end select

    initialphi = phii

    RETURN
  END FUNCTION initialphi


  FUNCTION getEps(phi,dphi)
    !
    !     Returns epsilon given phi and dphi/dalpha
    !     slowroll parameter epsilon_H = 2 M_pl^2 [(dH/dphi)/H]^2
    !
    REAL*8 :: getEps,phi,dphi
    getEps=2.d0*(M_Pl**2)*(((getHdot(phi,dphi)/dphi)/getH(phi,dphi))**2)
    RETURN

  END FUNCTION getEps


  FUNCTION getH(phi,dphi)
    !
    !     Returns H given phi and dphi/dalpha
    !
    REAL*8 :: getH,phi,dphi

    getH=SQRT(pot(phi)/3./M_Pl/M_Pl/ &
         &     (1.0-dphi*dphi/6.0/M_Pl/M_Pl))
    RETURN
  END FUNCTION getH


  FUNCTION getHJH(phi)

    !
    ! returns H given phi, only used for slow-roll reconstruction
    !

    REAL*8 :: getHJH, phi, A1, A2, A3, eps, eta, xi
    
    if (potential_choice .eq. 7) then
       eps=vparams(1)
    else if (potential_choice .eq. 8) then
       eps=10.d0**vparams(1)
    end if
    eta=vparams(2)
    xi=vparams(3)
    
    A1=sqrt(eps/2.)
    A2=eta/4.
    A3=xi/(12.*sqrt(eps*2))

    getHJH = H0*(1.+A1*phi+A2*phi**2+A3*phi**3)
  END FUNCTION getHJH



  FUNCTION getHJHdot(phi,dphi)
    !
    !     Returns dH/dalpha given phi and dphi/dalpha, only used for slow-roll reconstruction
    !
    REAL*8 :: getHJHdot,phi,dphi

    getHJHdot=-dphi*dphi*getHJH(phi)/2./M_Pl/M_Pl
    RETURN
  END FUNCTION getHJHdot


  FUNCTION getHJeps(phi)

    !
    ! returns epsilon given phi, only used for slow-roll reconstruction
    !

    REAL*8 :: getHJeps, phi, A1, A2, A3, eps, eta, xi    
    if (potential_choice .eq. 7) then
       eps=vparams(1)
    else if (potential_choice .eq. 8) then
       eps=10.d0**vparams(1)
    end if
    eta=vparams(2)
    xi=vparams(3)
    
    A1=sqrt(eps/2.)
    A2=eta/4.
    A3=xi/(12.*sqrt(eps*2))

    getHJeps = 2.*((A1+2.*A2*phi+3.*A3*phi**2)/(1.+A1*phi+A2*phi**2+A3*phi**3))**2
  END FUNCTION getHJeps


  FUNCTION getHJeta(phi)

    !
    ! returns eta given phi, only used for slow-roll reconstruction
    !

    REAL*8 :: getHJeta, phi, A1, A2, A3, eps, eta, xi
    
    if (potential_choice .eq. 7) then
       eps=vparams(1)
    else if (potential_choice .eq. 8) then
       eps=10.d0**vparams(1)
    end if
    eta=vparams(2)
    xi=vparams(3)
    
    A1=sqrt(eps/2.)
    A2=eta/4.
    A3=xi/(12.*sqrt(eps*2))

    getHJeta = 2.*((2.*A2+6.*A3*phi)/(1.+A1*phi+A2*phi**2+A3*phi**3))

  END FUNCTION getHJeta

  FUNCTION getHJxi(phi)

    !
    ! returns xi given phi, only used for slow-roll reconstruction
    !

    REAL*8 :: getHJxi, phi, A1, A2, A3, eps, eta, xi
    
    if (potential_choice .eq. 7) then
       eps=vparams(1)
    else if (potential_choice .eq. 8) then
       eps=10.d0**vparams(1)
    end if
    eta=vparams(2)
    xi=vparams(3)
    
    A1=sqrt(eps/2.)
    A2=eta/4.
    A3=xi/(12.*sqrt(eps*2))


    getHJxi = 4.*((A1+2.*A2*phi+3.*A3*phi**2)*6*A3/(1.+A1*phi+A2*phi**2+A3*phi**3)**2)

  END FUNCTION getHJxi
  

  FUNCTION getHdot(phi,dphi)
    !
    !     Returns dH/dalpha given phi and dphi/dalpha
    !
    REAL*8 :: getHdot,phi,dphi

    getHdot=-dphi*dphi*getH(phi,dphi)/2./M_Pl/M_Pl
    RETURN
  END FUNCTION getHdot


  FUNCTION getdepsdalpha(phi,dphi)
    !
    !    Returns depsilon/dalpha given phi and dphi/dalpha
    !    Gets this by differentiating Peiris et al Eq A3 (2003)
    !
    REAL*8 :: getdepsdalpha,phi,dphi,H,dHdalpha,eps
    
    H=getH(phi,dphi)
    dHdalpha=getHdot(phi,dphi)
    eps=getEps(phi,dphi)
    getdepsdalpha=3./H/H*(2.*H*dHdalpha*(1.-eps/3.)-dVdphi(phi)*dphi/3./M_Pl/M_Pl)
    RETURN
  END FUNCTION getdepsdalpha

  
  FUNCTION powerspectrum(u1,u2,z)
    USE internals
    REAL*8 :: powerspectrum,u1,u2
    REAL*8 :: z
    !
    !     Calculates P_R(k) given u1, u2, z
    !
    powerspectrum = (u1*u1/(2.0*k)+ &
         &           u2*u2*k/2.0/(h_ik*a_ik)**2) &
         &           /z**2 &
         &           *(k**3)/2./PI**2
    RETURN
  END FUNCTION powerspectrum

  FUNCTION norm_u(u1,u2)
    USE internals
    REAL*8 :: norm_u,u1,u2
    !
    !     Calculates P_R(k) given u1, u2, z
    !
     norm_u= sqrt((u1*u1/(2.0*k)+ u2*u2*k/2.0/(h_ik*a_ik)**2))
    RETURN
  END FUNCTION norm_u


  FUNCTION tensorpower(v1,v2,a)
    USE internals
    REAL*8 :: tensorpower,v1,v2,a
    !
    !     Calculates P_h(k) given v1, v2
    !
    tensorpower = (v1*v1/(2.0*k)+ &
         &           v2*v2*k/2.0/(h_ik*a_ik)**2) &
         &           /a**2 &
         &           *(k**3)*4./PI**2/(M_Pl**2)
    RETURN
  END FUNCTION tensorpower

END MODULE potential
