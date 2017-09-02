!     Bond stretch hyperelastic model
!     Functional form is the same as original model in ConstitutiveBondStretch.f
!     Interface is different, for use with internal energy based phase field model 
!
!     nprops=4, ninternal=3
!     
      module matdata

      integer, parameter :: nmatprops = 6
      integer, parameter :: ninternal = 3
      end module matdata
****************************************************************************
      subroutine Constitutive(props, F, gfunc, TR, varepsilon0, psi,
     &     Tangents, internal, pnewdt)

      use matdata
      
      implicit none

!     input
      real(8) :: props(nmatprops), F(3,3), gfunc

!     output
      real(8) :: TR(3,3), varepsilon0, psi
      real(8) :: Tangents(3,3,3,3), pnewdt

!     input and output
      real(8) :: internal(ninternal)

!     local variables
      real(8) :: Tangents_dev(3,3,3,3), Tangents_vol(3,3,3,3)
      real(8) :: TR_dev(3,3), TR_vol(3,3), varepsilon_star
      real(8) :: varepsilon_bond, varepsilon_vol, varepsilon
      real(8) :: Finv(3,3), detF, Iden(3,3), Ebdam, k_bulk, k_vol
      real(8) :: Gshear,Kbulk,Ebar,lmd_L,PPFG,PPFGD,PPGG,beta,dbdx
      real(8) :: x,Jm23,psi_ent,dpsidJ,d2psidJ2
      real(8) :: A1,A2,A3,A4,A5,A6,dA5,dlmdbdI1,dy1dI1,dPsiDevdI1,d2bdx2
      real(8) :: Finvtrans(3,3),I1,lmd_b0,lmd_b,lmd_bar
      real(8) :: dA1dI1,dA2dI1,dA3dI1,dA4dI1,dA5dI1,dA6dI1,dlmdb2dI12
      real(8) :: dy12dI12,dpsidevdI12,dI1dF(3,3),d2I1dF2,K2

      integer :: status,i,j,k,l

      call onem(Iden)
      
!     rename material properties
      k_bond = props(1)
      k_vol  = props(2)
      Kbulk  = props(3)
      Gshear = props(4)
      lmd_L  = props(5)
      Ebar   = props(6)

!     rename internal variables
      lmd_b0 = internal(1)
      lmd_b  = internal(2)

!     Deformation quantities
      call mdet(F, detF)
      Jm23 = detF**(-2.d0/3.d0)
      call m3inv(F, Finv)
      Finvtrans = transpose(Finv)
      I1 = Jm23*sum(F*F)
      lmd_bar = dsqrt(I1/3.d0)

      Ebdam = (gfunc+k_bond)*Ebar
      call SolveBondStretch(lmd_bar,Ebdam,Gshear,lmd_L,lmd_b,status)
      if (status .ne. 0) then
         pnewdt = 0.25d0
         return
      end if

      x = lmd_bar/(lmd_L*lmd_b)
      
      call LangInv(x,beta,dbdx,d2bdx2)

!     The undamaged internal energy is needed as the fracture
!     driving force in the phase field balance law.      

!     undamaged internal energy per bond, normalizeb by (kb*theta)
      varepsilon_star = 0.5d0*Ebar/Gshear*(lmd_b - 1.d0)**2

!     undamaged internal energy density due to bond stretch
      varepsilon_bond = varepsilon_star*Gshear*lmd_L**2
!     undamaged volumetric internal energy
      varepsilon_vol = Kbulk/8.d0*(detF**2 + detF**(-2) - 2.d0)

!     total undamaged internal energy
      varepsilon0 = varepsilon_bond + varepsilon_vol

!     entropic deviatoric free energy density
      psi_ent = Gshear*(
     &     x*beta + dlog(2.d0*beta)- beta
     &     - dlog(1.d0 - dexp(-2.d0*beta))
     &     )

!     total free energy density
      psi = (gfunc+k_bond)*varepsilon_bond+(gfunc+k_vol)*varepsilon_vol
     $     + psi_ent

!     The following computations are for the actual (i.e., damaged)
!     Piola stress and tangent operator.
!     For convenience, reduce the bond stiffness and bulk
!     modulus by the damage factor.
      Kbulk = (gfunc+k_vol)*Kbulk
      Ebar  = (gfunc+k_bond)*Ebar

      dpsidJ = 0.25d0*Kbulk*(detF - detF**(-3)) 
      d2psidJ2 = 0.25d0*Kbulk*(1.d0 + 3.d0/detF**4)


!     Stress (including damage)
!     first volumetric part
      TR_vol = dpsidJ*detF*Finvtrans

!     Deviatoric stress
      !!! A5 = d (varepsilon)/ d (lmdb)
      !!! dA5 = d (A5)/ d (lmdb)
      !!! PPFG = d (lmdb*d(varepsilon)/d(lmdb))/d(lmdb)
      !!! Also can be expressed as d (lmdb*A5)/d(lmdb)
      !!! PPFGD = d (PPFG)/d lmdb
      !!! PPGG = d (PPFG*lmdb)/dlmdb
      A5 = Ebar*(lmd_b-1.0)
      dA5 = Ebar 
      PPFG = Ebar*(2.d0*lmd_b-1.d0)
      PPFGD = Ebar*(2.d0)
      PPGG = Ebar*(4.d0*lmd_b-1.0)

      !!! following quantities do not depend on the form of enthalpy.	  
      A1 = Gshear*(beta+x*dbdx)/(6.d0*lmd_L*lmd_bar)
      A2 = PPFG*lmd_b+Gshear*(beta+x*dbdx)*x
      A3 = (PPFG)/(Gshear*(beta+x*dbdx))
      A4 = (2.d0*dexp(-2.d0*beta)/(1.d0-dexp(-2.d0*beta)))
      A6 = Gshear*(beta+x*dbdx+dbdx/beta-dbdx-
     +                          A4*dbdx)
      dlmdbdI1 = A1/A2
      dy1dI1 = A3*dlmdbdI1
      dPsiDevdI1 = (lmd_L**2)*(A5*dlmdbdI1+A6*dy1dI1)
      dI1dF = 2.d0*(detF**(-2.d0/3.d0)*F-I1/3.d0*Finvtrans)
      
      TR_dev = dPsiDevdI1*(2.d0*detF**(-2.d0/3.d0)*F -
     &     2.d0/3.d0*I1*Finvtrans)

      TR = TR_dev + TR_vol
      
!     Tangents

      dA1dI1 = Gshear*(2.d0*dbdx+x*d2bdx2)*dy1dI1/(6.d0*lmd_L*lmd_bar)
     +     - Gshear*(beta+x*dbdx)/((6.d0*lmd_L*lmd_bar)**2)*
     +     (6.d0*lmd_L/(6.d0*lmd_bar))

      dA2dI1 = Gshear*(x*dbdx+beta+2.d0*x*dbdx+x*x*d2bdx2)*dy1dI1
     +     + PPGG*dlmdbdI1
      dA3dI1 = (PPFGD*dlmdbdI1)/(Gshear*(beta+x*dbdx))
     +     - (PPFG)/((Gshear*(beta+x*dbdx))**2)*
     +     (Gshear*(2.d0*dbdx+x*d2bdx2)*dy1dI1)
      dA4dI1 = (-4.d0*dexp(-2.d0*beta)/(1.d0-dexp(-2.d0*beta)))
     +     *dbdx*dy1dI1
     +     - (2.d0*dexp(-2.d0*beta)/((1.d0-dexp(-2.d0*beta))**2))*
     +     (2.d0*dexp(-2.d0*beta)*dbdx)*dy1dI1
      dA5dI1 = dA5*dlmdbdI1
      dA6dI1 = Gshear*(2.d0*dbdx*dy1dI1+x*d2bdx2*dy1dI1
     +     + d2bdx2*dy1dI1/beta-dbdx*dbdx*dy1dI1/(beta*beta)
     +     - d2bdx2*dy1dI1-dA4dI1*dbdx-A4*d2bdx2*dy1dI1)

      dlmdb2dI12 = dA1dI1/A2-(A1/(A2**2))*dA2dI1
      dy12dI12 = A3*dlmdb2dI12+dlmdbdI1*dA3dI1
  
      dpsidevdI12 = (lmd_L**2)*(dA5dI1*dlmdbdI1+A5*dlmdb2dI12 +
     +     dA6dI1*dy1dI1+A6*dy12dI12)
      
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3

                  d2I1dF2 = 2.d0*(-2.d0/3.d0*Jm23*Finv(l,k)*F(i,j)
     &                 + Jm23*Iden(i,k)*Iden(j,l)
     &                 - 1.d0/3.d0*dI1dF(k,l)*Finv(j,i)
     &                 + I1/3.d0*Finv(j,k)*Finv(l,i))
                  
                  Tangents_dev(i,j,k,l) = dPsiDevdI1*d2I1dF2
     &                 + dpsidevdI12*dI1dF(i,j)*dI1dF(k,l)
                  
                  Tangents_vol(i,j,k,l) =
     &                 (d2psidJ2*detF+dpsidJ)*detF*Finv(l,k)*Finv(j,i)
     &                 - dpsidJ*detF*Finv(j,k)*Finv(l,i)
               end do
            end do
         end do
      end do

      Tangents = Tangents_dev + Tangents_vol
      
      internal(2) = lmd_b
      internal(3) = varepsilon_star
      
      end subroutine Constitutive
****************************************************************************
      subroutine SolveBondStretch(lmd_bar,Eb,Gshear,lmd_L,root,status)

      implicit none

!     in
      real(8) :: lmd_bar, Eb, Gshear, lmd_L
!     in, out
      real(8) :: root
!     out
      integer :: status

!     local
      real(8) :: f,df,fl,fh,xl,xh,x1,x2,swap,dx,dxold
      real(8) :: rootMin,rootMax,temp
      integer :: j
      
      integer, parameter :: maxit = 50
      real(8), parameter :: xacc = 1.0d-9, zero = 0.0d0, one = 1.0d0

!     initialize status to normal termination
!     We'll change it if an error occurs
      status = 0
      
      rootMin = max(1.d0, (1.d0+1.d-12)*lmd_bar/lmd_L)
      rootMax = 1.6d0*lmd_bar

!     use old value of bond stretch as initial guess,
!     unless the old value falls outside the new bounds.
!     In that case, take the mean of the bounds
      if ( (root .lt. rootMin) .or. (root .gt. rootMax) ) then
         root = 0.5d0*(rootMin + rootMax)
      end if
      
!     set up bounds on root
      x1 = rootMin
      x2 = rootMax
      call lmdbFunc(x1,FL,DF,Eb,Gshear,lmd_bar,lmd_L)
      call lmdbFunc(x2,FH,DF,Eb,Gshear,lmd_bar,lmd_L)
      if(fl*fh.ge.zero) then
         write(*,*) 'Warning, root not bounded'
!        indicate an error has occurred to caller
         status = 1
         return
      endif

C
C		ORIENT THE SEARCH SO THAT F(XL) < 0.
C
      IF( FL .LT. 0.D0 ) THEN
         XL = X1
         XH = X2
      ELSE
         XH = X1
         XL = X2
         SWAP = FL
         FL = FH
         FH = SWAP
      END IF
C
C		INITIALIZE THE GUESS FOR THE ROOT, THE ''STEP SIZE
C		BEFORE LAST'', AND THE LAST STEP
C
c      ROOT = (set by input argument)
      DXOLD = DABS(X2 - X1)
      DX = DXOLD

      call lmdbFunc(root,F,DF,Eb,Gshear,lmd_bar,lmd_L)
C
C			LOOP OVER ALLOWED ITERATIONS
C
      DO 10 J = 1,MAXIT
C
C			BISECT IF NEWTON OUT OF RANGE, OR NOT DECREASING
C			FAST ENOUGH.
C
         IF( ((ROOT-XH)*DF - F)*((ROOT - XL)*DF -F) .GE. 0.D0
     +        .OR. DABS(2.D0*F) .GT. DABS(DXOLD*DF) ) THEN

            DXOLD = DX
            DX = 0.5D0*(XH-XL)
            ROOT = XL + DX
            IF( XL .EQ. ROOT ) THEN
C
C			CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         ELSE
C
C			NEWTON STEP IS ACCEPTABLE. TAKE IT.
C
            DXOLD = DX
            DX = F/DF
            TEMP = ROOT
            ROOT = ROOT - DX
            IF( TEMP .EQ. ROOT) THEN
C
C			 CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         END IF
C
C		CONVERVEGENCE CRITERION
C
         IF( DABS(DX) .LT. XACC) RETURN

C
C			THE ONE NEW FUNCTION EVALUATION PER ITERATION
C
         call lmdbFunc(root,F,DF,Eb,Gshear,lmd_bar,lmd_L)

C
C		MAINTAIN THE BRACKET ON THE ROOT
C
         IF( F .LT. 0.D0) THEN
            XL = ROOT
            FL = F
         ELSE
            XH = ROOT
            FH = F
         END IF

 10   CONTINUE

      WRITE(*,'(/1X,A)') 'SolveBondStretch EXCEEDING MAXIMUM ITERATIONS'
      status = 1

      end subroutine SolveBondStretch
**********************************************************************
      subroutine lmdbFunc(x,f,df,Eb,G0,lmd_bar,lmd_L)

      ! This subroutine serves as the function we would like to solve for
      ! the bond stretch

      implicit none

      real(8) :: f,df,x,y,dydx,beta,dbdy,dbeta
      real(8) :: Eb,G0,lmd_bar,lmd_L,AAAX,unused

      ! Compute the residual
      !
      y = lmd_bar/(lmd_L*x)
      call LangInv(y,beta,dbdy,unused)
      AAAX = Eb
      !f = AAAX*(dlog(phi)-(G0/Eb)*yyt*beta)
      f = AAAX*(x*(x-1.d0)-(G0/Eb)*y*beta)

      ! Compute the tangent
      !
      dydx = -lmd_bar/(lmd_L*(x**2))
      dbeta = dbdy*dydx
      !dff = AAAX*(1.d0/x
      !+          -(G0/Eb)*dydx*beta
      !+          -(G0/Eb)*y*dbeta)
      df = AAAX*(2.d0*x-1.d0
     +          -(G0/Eb)*dydx*beta
     +          -(G0/Eb)*y*dbeta)

      return
      end subroutine lmdbFunc
****************************************************************************
c     Inverse Langevin function
c      
      subroutine LangInv(y, f, df, ddf)
      implicit none

!     in
      real(8) :: y
!     out
      real(8) :: f, df, ddf
      
      y = min(y, 1.d0 - 1.d-12)
      f = y*(3.d0-y**2)/(1.d0-y**2)
      df = (3.d0-y**2+4.d0*y**2/(1.d0-y**2))/(1.d0-y**2)
      ddf = (6.d0*y)/(y**2 - 1.d0) - (8.d0*y**3)/(y**2 - 1.d0)**2 +
     &     (8.d0*y**3*(y**2 - 3.d0))/(y**2 - 1.d0)**3 -
     &     (6.d0*y*(y**2 - 3.d0))/(y**2 - 1.d0)**2
      
      end subroutine LangInv
****************************************************************************
      subroutine InitializeInternalVariables(props,internals)

      use matdata
      
      implicit none

      integer :: stat
      real(8) :: internals(ninternal), props(nmatprops), Kbulk
      real(8) :: lmd_b0, lmd_b, varepsilon_b, Ebar, Gshear, lmd_L

      if (ninternal .ne. 3) then
         write(*,*) 'Incorrect number of internal variables.'
         write(*,*) 'Check number of state variables',
     &        ' assigned in user element'
         call xit
      end if

      Kbulk  = props(1)
      Gshear = props(2)
      lmd_L  = props(3)
      Ebar   = props(4)
      
c     initial bond stretch
      call SolveBondStretch(1.d0,Ebar,Gshear,lmd_L,lmd_b0,stat)
      if (stat .ne. 0) then
         write(*,*) 'Initialization of bond stretch failed'
         call xit
      end if

c     current bond stretch
      lmd_b = lmd_b0

c     bond stretch energy
      varepsilon_b = 0.d0

      internals(1) = lmd_b0
      internals(2) = lmd_b
      internals(3) = varepsilon_b
      
      end subroutine InitializeInternalVariables
****************************************************************************
