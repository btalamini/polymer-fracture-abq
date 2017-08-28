#define STAGGERED 1
#define TEST_MODE 0
!     set ONLY_DAMAGE_IN_TENSION to 1 for the asymmetric tension/compresion model
!     set STAGGERED to 1 to use with the staggered solution approach.
!     (Use with "*Solution technique, type=separated" in input file)
!
!     set TEST_MODE to 1 to remove all of the includes of ABA_PARAM.INC
!     This is useful for sompiling the UEL for debugging.
************************************************************************
#include "ElementRoutines.for"
#include "UtilityRoutines.for"
#include "ConstitutiveBondStretchEnthalpic.f"

      
!     
!     User element for the threshold phase field fracture model. The model
!     is also time dependent (Allen-Cahn type). The temperature degree of
!     freedom is used for the phase field.
!     The underlying constiutive law is the bond stretch 8 chain model.
! 
!     2D plane stress only.
!
!     Nodal unknowns:
!     1 U1
!     2 U2
!     11 NT11 (for phase field)
!
!     Supported elements:
!     4 node bilinear quad
!     8 node quadratic serendipity quad
!     9 node biquadratic quad
!     3 node linear triangle
!     6 node quadratic triangle      
!      
!     The damage field and the displacement have the same interpolation.
!
!     You MUST use Abaqus 6.10-1 or higher for this to work.
!     
!     Brandon Talamini, 2015-2016, Implemented in Abaqus 6.14
!
***********************************************************************
!
!     User element statement in the input file (set ? values as needed):
!     2D case:
!     *User Element,Nodes=?,Type=U1,Iproperties=2,Properties=8,Coordinates=2,Variables=?
!     1,2,11
!
!     State Variables
!     --------------------------------------------------------------
!     Global SDV's (for use in UVARM to make contour plots)
!       1) phi, phase field (stored to control dummy element deletion)
!       2) Cauchy Stress T(1,1)
!       3) Cauchy Stress T(2,2)
!       4) Cauchy Stress T(3,3)
!       5) Cauchy Stress T(1,2)
!       6) Cauchy Stress T(2,3)
!       7) Cauchy Stress T(3,1)
!       8) pressure, -tr(T)/3
!       9) lambda_bar = sqrt(I1/3), where I1 is J**(-2/3) trace(C)
!      10) J, determinant of deformation gradient
!      11) lambda33, out of plane stretch
!      12) vareps0, undamaged internal energy density
!      13) Actual strain energy density (including damage)
!      14) psi_tau, Free energy density
!      15) k*vareps0, regularization energy
!      16) lambda_b, bond_stretch
!
!     Local SDV's (for use internal to the UEL)
!      svars(1+j)  = vareps0, undamaged internal energy, (maximum over time)
!      svars(2+j)  = lam3,   out of plane stretch
!      svars(3+j)  = vd,     total energy dissipated through rate-dependent fracture term
!      + number of internal variables of constitutive model
!
!     In the input file, set
!     * 'ngSdv'= number of global SDV's (16)
!     * 'nlSdv'= number of local SDV's (3 + #internals)
!     * 'User output variables'= ngSdv      
!     * 'variables'=(nlSdv*nIntPt)
!
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     vareps_c  = props(1)  critical value of the internal energy
!     len    = props(2)  regualrization length
!     beta   = props(3)  kinetic modulus for phase field evolution
!     k      = props(4)  Rest stiffness, ratio of remaining stiffness at failure
!     Kbulk  = props(5)  Bulk modulus
!     Gshear = props(6)  Young's modulus
!     lmd_L  = props(7)  Locking stretch
!     Ebar   = props(8)  Segment stretching modulus
!     --------------------------------------------------------------
!     Integer properties
!     nlSdv  = jprops(1) 
!     ngSdv  = jprops(2)
***********************************************************************

      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  nInt
      !   Total number of integration points per elemnet. The same
      !   number for the real and dummy elements.  You must set that
      !   parameter value here.
      !
      ! nDim
      !  Spatial dimension. You must set that parameter here.
      !
      !  elemPt
      !   Element pointer
      !
      !  elOffset
      !   Offset between user element ID and dummy element ID
      !
      !  debug_wait
      !   Set to true to force UEL into infinite loop on its
      !   first call, so that we can catch the program in the 
      !   debugger.      

      integer elemPt,numUvarm,elCount,err
      integer, parameter :: numElem = 4
      integer, parameter :: nInt = 1
      integer, parameter :: nDim = 2
      integer, parameter :: elOffset = 100000
      logical, save :: debug_wait = .false.

      real(8), allocatable :: globalSdv(:,:,:)

      
      end module global

***********************************************************************
***********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.
     
      use global

#if !TEST_MODE
      include 'ABA_PARAM.INC'
#endif
c      IMPLICIT NONE
     
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

c     match up the dummy visualization element with
c     the corresponding user element
      elCount = noel - elOffset

c     copy the globally stored output variables
c     into the user variables for this element and
c     quad point.
      if(allocated(globalSdv)) then
         uvar(1:nuvarm) = globalSdv(1:nuvarm,npt,elCount)
      end if

      return
      end subroutine uvarm

****************************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)

      USE global
      USE matdata, ONLY: nmatprops, ninternal
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
      
      real*8 u(nNode,nDim),du(nNode,nDim),uOld(nNode,nDim),v(nNode,nDim)
      real*8 phi(nNode),dphi(nNode), phiOld(nNode)

      integer i,j,k,l,intpt,a11,b11,jj,ai
      integer a,b,a12,b12,stat,ngSdv
      integer nlSdv,lenJobName,lenOutDir
      integer ndof_node,row,col

      real(8) :: Ru(nDim*nNode),Rd(nNode),Kdd(nNode,nNode)
      real(8) :: Kud(nDim*nNode,nNode),Kuu(nDim*nNode,nDim*nNode)
      real(8) :: xi(nInt,nDim),w(nInt),sh(nNode),dshxi(nNode,nDim)
      real(8) :: detMapJ,dsh(nNode,nDim),dvR
      real(8) :: X(mcrd),u_tau(nDim),phi_tau,phi_t,F_tau(3,3)
      real(8) :: grad_phi(nDim),phi_dot,detF,gfunc
      real(8) :: vareps_c,len,beta,kappa,lambda_b
      real(8) :: TR_tau(3,3),T_tau(3,3),Tangents(3,3,3,3)
      real(8) :: xip(nDim),varpi,psi_tau,I1,lambda_bar
      real(8) :: d2psi_dF_dphi(3,3),H_tau,strain_energy_density
      real(8) :: body(nDim),Iden(3,3),rest_stiffness
      real(8) :: vareps0,vareps0_t,vareps0_tau,P33,lam33,press,C33
      real(8) :: matprops(nmatprops),dlam33,dlam33old,lam33L,lam33H
      real(8) :: P33L,P33H,detF2D,limitA,limitB,swap,res0,vd
      real(8) :: lam33Min,lam33Max
      real(8), dimension(:), allocatable :: internals, internals_old

      real(8), parameter :: zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,
     &     Pi=4.d0*atan(1.d0),three=3.d0,third=1.d0/3.d0

      character*256 jobName,outDir,fileName

      ! infinite loop to hold the program while debugger is attached
      ! set debug=.false. to use code normally
      do while (debug_wait)
         continue
      end do

c     ------------------------------------------------
c     Check input data
      
      ! Check lflags to make sure correct abaqus
      ! options have been chosen, otherwise abort
      Call KCheckLflags(lflags)

      ! Check number of supplied real number properties
      if (nprops .ne. (4+nmatprops)) then
         write(*,*) 'Incorrect number of real number properties',
     &        ' for phase field element'
         call xit
      end if

      ! Check number of supplied integer properties
      if (njprop .ne. 2) then
         write(*,*) 'Incorrect number of integer properties',
     &        ' for phase field element'
         call xit
      end if

c     End check of input data
c     ------------------------------------------------
      
      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.zero) return
      
      ! phase field properties
      vareps_c = props(1)
      len      = props(2)
      beta     = props(3)
      rest_stiffness = props(4)
      kappa = two*vareps_c*len**2

      ! mechanical material properties
      matprops = props(5:nprops)


      ! handle steady state case
      !
      ! The easiest way to handle this is
      ! to set beta=0
      if (lflags(1) .eq. 71) beta = 0.d0  

      ! Get element integer parameters
      !
      nlSdv  = jprops(1) ! number of local sdv's per integ point
      ngSdv  = jprops(2) ! number of global sdv's per integ point
      ndof_node = ndofel / nnode ! dofs per node (assumes equal interpolations for all fields)
      !ninternal = nlSdv - 3     ! internal variables for constitutive law,
                                ! sdv(1) is psi0_t
                                ! sdv(2) is lam33
                                ! sdv(3) is vd
      
      if (nsvars .ne. (ninternal+3)*nInt) then
         write(*,*) 'Incorrect number of solution dependent',
     &        ' variables allocated.'
         write(*,*) 'check the VARIABLES option with the *USER ELEMENT',
     &        ' card.'
         call xit
      end if
      
      allocate(internals(ninternal))
      allocate(internals_old(ninternal))

      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
         !
         ! allocate memory for the globalSdv's
         !
         ! numElem needs to be set in the MODULE
         ! nInt needs to be set in the UEL
         !
         stat=0
         allocate(globalSdv(ngSdv,nInt,numElem),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt=',nInt
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv
         write(*,*) '-------------------------------------------------'
      endif

      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rd  = zero

      Kuu = zero
      Kud = zero
      Kdd = zero
      
      Energy = zero          ! initialize strain energy

      ! no body forces for now
      body = zero
      

      ! Extract nodal unknowns from abaqus array
      ! displacements, velocities, order parameter
      ! (velocity from backward difference)
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
            v(i,j) = du(i,j)/dtime
         enddo
         k = k + 1
         phi(i) = Uall(k)
         dphi(i) = DUall(k,1)
         phiOld(i) = phi(i) - dphi(i)
      enddo

      ! Identity tensor
      !
      call onem(Iden)

      !
      ! Obtain integration point local coordinates and weights
      !
      call QuadratureRule(nDim, nNode, nInt, xi, w)

      
      !----------------------------------------------------------------
      ! Begin the loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nInt

         if ((kstep .eq. 1) .and. (kinc .le. 1)) then
            ! initialize
            vareps0_t = 0.d0
            lam33 = one
            vd = 0
            svars(1+jj) = vareps0_t
            svars(2+jj) = lam33
            svars(3+jj) = vd
            if (ninternal .gt. 0) then
               call InitializeInternalVariables(matprops,internals)
c              svars(4+jj:4+jj+ninternal-1) = internals
               do i=1,ninternal
                  svars(3+jj+i) = internals(i)
               end do
            end if
         else
            ! Obtain state variables from previous increment
            !
            vareps0_t  = svars(1+jj)
            lam33   = svars(2+jj)
            vd      = svars(3+jj)
            if (ninternal .gt. 0) then
c              internals = svars(4+jj:4+jj+ninternal-1)
               do i=1,ninternal
                  internals(i) = svars(3+jj+i)
               end do
            end if
         end if

         ! Obtain shape functions and their gradients
         !
         call CalcShape(nDim, nNode, nInt, intpt, xi, sh, dshxi)
         
!        Map shape function derivatives from local to
!        global reference coordinate system
         call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)

         ! abort calculation if element is not well formed
         if(stat.eq.0) then
            write(*,*) 'Invalid element with negative Jacobian 
     &detected. Aborting. Element number ',jelem 
            call Xit
         endif

         ! discrete volume element
         dvR = detMapJ*w(intpt)

         ! interpolate position and displacement
         X=zero
         u_tau=zero
         do a=1,nNode
            do i=1,nDim
               X(i) = X(i) + sh(a)*coords(i,a)
               u_tau(i) = u_tau(i) + sh(a)*u(a,i)
            end do
         end do


         ! Interpolate order parameter and its gradient
         phi_tau = zero
         phi_t = zero
         phi_dot = zero
         grad_phi = zero
         do k=1,nNode
            phi_tau = phi_tau + phi(k)*sh(k)
            phi_t   = phi_t + phiOld(k)*sh(k)
            do i=1,nDim
               grad_phi(i) = grad_phi(i) + phi(k)*dsh(k,i)
            enddo
         enddo
         phi_dot = (phi_tau - phi_t)/dtime

         
         ! Interpolate deformation gradient
         F_tau = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
               enddo
            enddo
         enddo
         F_tau(3,3) = lam33
         
         ! degradation function
         gfunc = (1.0d0 - phi_tau)**2 + rest_stiffness

         internals_old = internals
         
         ! Constitutive relations for stress
         ! Iterate to find plane stress state
         detF2D = F_tau(1,1)*F_tau(2,2)-F_tau(1,2)*F_tau(2,1)
         limitA = 1.0d0/detF2D
         limitB = 1.d0
         lam33Min = 0.9d0*min(limitA,limitB)
         lam33Max = 1.1d0*max(limitA,limitB)

         if (lam33Min .lt. zero) lam33Min = 1.d-9
         
         internals = internals_old
         F_tau(3,3) = lam33Min
         call Constitutive(matprops,F_tau,gfunc,TR_tau,vareps0_tau,
     &        strain_energy_density,Tangents,internals,pnewdt)
         P33L = TR_tau(3,3)

         internals = internals_old
         F_tau(3,3) = lam33Max
         call Constitutive(matprops,F_tau,gfunc,TR_tau,vareps0_tau,
     &        strain_energy_density,Tangents,internals,pnewdt)
         P33H = TR_tau(3,3)

         if (P33L*P33H .gt. 0) then
            write(*,*) 'F33 not bounded for plane stress'
            write(*,*) 'element=',jelem
            write(*,*) 'step=',kstep,' inc=',kinc
            write(*,*) 'P33L=',P33L,' P33H=',P33H
            write(*,*) 'lam33Min=',lam33Min,' lam33Max=',
     &           lam33Max
            pnewdt = 0.25d0
            return
         end if

         if (P33L .lt. 0.d0) then
            lam33L = lam33Min
            lam33H = lam33Max
         else
            lam33H = lam33Max
            lam33L = lam33Min
            swap = P33L
            P33L = P33H
            P33H = swap
         end if
         
         internals = internals_old
         lam33 = 0.5d0*(lam33L+lam33H)
         F_tau(3,3) = lam33
         call Constitutive(matprops,F_tau,gfunc,TR_tau,vareps0_tau,
     &        strain_energy_density,Tangents,internals,pnewdt)

         C33 = Tangents(3,3,3,3)
         P33 = TR_tau(3,3)
         dlam33old = dabs(lam33H-lam33L)
         dlam33 = dlam33old
         res0 = dabs(P33)
         i = 1
         do while (dabs(P33) .gt. 1.d-8*res0)

            if (i .gt. 20) then
               write(*,*) 'WARNING:'
               write(*,*) 'Plane stress condition could not be enforced'
               write(*,*) 'Current rel error norm=',dabs(P33/res0)
               write(*,*) 'Element ',jelem,', GP ',intpt
               write(*,*) 'step ',kstep,' inc ',kinc
               pnewdt = 0.25d0
               exit
            end if

!           bisect if Newton out of range or decreasing slowly
            if ( ((lam33-lam33H)*C33 - P33)*((lam33-lam33L)*C33-P33)
     &           .GE. 0.D0
     +           .OR. DABS(2.D0*P33) .GT. DABS(dlam33old*C33) ) THEN
               
               dlam33old = dlam33
               dlam33 = 0.5D0*(lam33H-lam33L)
               lam33 = lam33L + dlam33
               
            else
c              NEWTON STEP IS ACCEPTABLE. TAKE IT.
c              
               dlam33old = dlam33
               dlam33 = P33/C33
               lam33 = lam33 - dlam33
            end if
            
            F_tau(3,3) = lam33

c           compute new residual
            internals = internals_old            
            call Constitutive(matprops,F_tau,gfunc,TR_tau,vareps0_tau,
     &           strain_energy_density,Tangents,internals,pnewdt)

            C33 = Tangents(3,3,3,3)
            P33 = TR_tau(3,3)

c           maintain brackets on root
            if ( P33 .LT. 0.D0) then
               lam33L = lam33
               P33L = P33
            else
               lam33H = lam33
               P33H = P33
            end if

            if (dabs(dlam33) .lt. 1.d-9) exit
            i = i + 1
         end do

         ! static condensation of mat tangent op
         do i=1,2
            do j=1,2
               do k=1,2
                  do l=1,2
                     Tangents(i,j,k,l) = Tangents(i,j,k,l) -
     &                    Tangents(i,j,3,3)*Tangents(3,3,k,l)/
     &                    Tangents(3,3,3,3)
                  end do
               end do
            end do
         end do

         ! Compute Jacobian determinant now that
         ! lambda33 is known
         call mdet(F_tau, detF)

c         d2psi_dF_dphi = -2.0d0*(1.0d0 - phi_tau) *
c     &        (Pstar_dev+Pstar_vol)


!        total free energy density
         psi_tau = strain_energy_density
     &        + 0.5d0*kappa*DOT_PRODUCT(grad_phi, grad_phi)
     &        + 2.0*vareps_c*phi_tau

!        viscous fracture dissipation density
         vd = vd + beta*phi_dot**2*dtime
         
!        update internal variable driving fracture
         vareps0_tau=MAX(vareps0_t, vareps0_tau)
#if STAGGERED
!        Numerical approximation:
!        use strain energy at time t instead of current time.
!        Decouples phase field balance law from momentum balance
         vareps0 = vareps0_t
#else
         vareps0 = vareps0_tau
#endif
         
!        Add a small value to the energy so that stiffness matrix
!        is not singular when there is no deformation
         vareps0 = vareps0 + 1e-8*vareps_c

!        generalized forces for phase field
         H_tau = max(vareps0 - vareps_c, zero)
         varpi = -two*(1.d0 - phi_tau)*H_tau + two*vareps_c*phi_tau
         xip = kappa*grad_phi

!        Output quantities
         T_tau = matmul(TR_tau, transpose(F_tau))/detF
         call tracem(T_tau, press)
         press = -press/3.d0
         I1 = sum(F_tau*F_tau)/detF**(two/three)
         lambda_bar = dsqrt(I1/3.d0)
      
         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         svars(1+jj) = vareps0_tau
         svars(2+jj) = lam33
         svars(3+jj) = vd
         do i=1,ninternal
            svars(3+jj+i) = internals(i)
         end do
c         if (ninternal .gt. 0) then
c            svars(4+jj:4+jj+ninternal-1) = internals
c         end if
         
         jj = jj + nlSdv ! setup for the next intPt


         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         globalSdv(1,intPt,jelem)  = phi_tau
         globalSdv(2,intPt,jelem)  = T_tau(1,1)
         globalSdv(3,intPt,jelem)  = T_tau(2,2)
         globalSdv(4,intPt,jelem)  = T_tau(3,3)
         globalSdv(5,intPt,jelem)  = T_tau(1,2)
         globalSdv(6,intPt,jelem)  = T_tau(2,3)
         globalSdv(7,intPt,jelem)  = T_tau(3,1)
         globalSdv(8,intPt,jelem)  = press
         globalSdv(9,intPt,jelem)  = lambda_bar
         globalSdv(10,intPt,jelem) = detF
         globalSdv(11,intPt,jelem) = lam33
         globalSdv(12,intPt,jelem) = gfunc*vareps0_tau
         globalSdv(13,intPt,jelem) = strain_energy_density
         globalSdv(14,intPt,jelem) = psi_tau
         globalSdv(15,intPt,jelem) = rest_stiffness*vareps0_tau
         lambda_b = internals(2)
         globalSdv(16,intPt,jelem) = lambda_b
         
         !  energies
         ! Strain energy
         Energy(2) = Energy(2) + strain_energy_density*dvR
         ! fracture energy - store in plastic dissipation variable
         Energy(4) = Energy(4) +
     &        (two*vareps_c*phi_tau +
     &         0.5d0*kappa*dot_product(grad_phi,grad_phi)
     &        )*dvR
         ! rate dependent fracture dissipation - store in viscous dissipation variable
         Energy(5) = Energy(5) + vd*dvR

         ! Compute/update the displacement residual vector
         !
         do a=1,nNode
            do i=1,nDim
               do j=1,nDim
                  Ru(nDim*(a-1)+i) = Ru(nDim*(a-1)+i) 
     &                 -TR_tau(i,j)*dsh(a,j)*dvR
     &                 +body(i)*sh(a)*dvR
               end do
            end do
         end do
         

         ! Compute phase field residual
         !
         do a=1,nNode
            Rd(a) = Rd(a) - (varpi+beta*phi_dot)*sh(a)*dvR
            do i=1,nDim
               Rd(a) = Rd(a) - xip(i)*dsh(a,i)*dvR
            end do
         end do

         ! Compute/update the displacement tangent matrix
         !
         do a=1,nNode
            do i=1,nDim
               row=nDim*(a-1)+i
               do b=1,nNode
                  do k=1,nDim
                     col=nDim*(b-1)+k
                     do j=1,nDim
                        do l=1,nDim
                           Kuu(row,col) = Kuu(row,col) +
     &                          tangents(i,j,k,l)*dsh(b,l)*dsh(a,j)*dvR
                        end do
                     end do
                  end do
               end do
            end do
         end do


!        Set off-diagonal terms of stiffness.
!        These terms are ignored by Abaqus if
!        the separated (staggered) solution procedure
!        is chosen in the input file.
         do a=1,nNode
            do i=1,nDim
               ai=nDim*(a-1)+i
               do b=1,nNode
                  do j=1,nDim
                     Kud(ai,b)=Kud(ai,b)
     &                    + dsh(a,j)*sh(b)*d2psi_dF_dphi(i,j)*dvR
                  end do
               end do
            end do
         end do

         do a=1,nNode
            do b=1,nNode
               Kdd(a,b) = Kdd(a,b) +
     &              (2.d0*(H_tau + vareps_c) + beta/dtime) *
     &              sh(a)*sh(b)*dvR
               do i=1,nDim
                  Kdd(a,b)=Kdd(a,b) + kappa*dsh(a,i)*dsh(b,i)*dvR
               end do
            end do
         end do
         
      enddo
      !
      ! End the loop over integration points
      !----------------------------------------------------------------

c     release memory from internal variable arrays
      deallocate(internals)
      deallocate(internals_old)

      !----------------------------------------------------------------
      ! Add in any stabilization stiffness and forces to what you 
      !  already have for the element, this could be hourglass, etc..
      !
      !
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.  
      !
      ! Return Abaqus the right hand side vector
      !
      rhs(:,1) = zero
      do i=1,nNode
         A11 = (ndof_node)*(i-1)
         A12 = nDim*(i-1)
         !
         ! displacement
         !
!         rhs(A11:A11+nDim,1) = Ru(A12:A12+nDim)
!        rhs(A11+1,1) = Ru(A12+1)
         do k=1,nDim
            rhs(A11+k,1) = Ru(A12+k)
         end do
         !
         ! damage
         !
         rhs(A11+nDim+1,1) = Rd(i)
         !
      enddo
      !
      ! Return Abaqus the tangent matrix
      !
      amatrx = zero
      do i=1,nNode
         do j=1,nNode
            A11 = (ndof_node)*(i-1)
            A12 = nDim*(i-1)
            B11 = (ndof_node)*(j-1)
            B12 = nDim*(j-1)
            !
            ! displacement
            !
            do k=1,nDim
               do l=1,nDim
                  amatrx(A11+l,B11+k) = Kuu(A12+l,B12+k)
               end do
            end do
            !
            ! displacement - order parameter coupling
            do k=1,nDim
               amatrx(A11+k,B11+nDim+1) = Kud(A12+k,j)
               amatrx(A11+nDim+1,B11+k) = Kud(B12+k,i)
            end do
            !amatrx(A11,B11+2)   = Kud(A12,j)
            !amatrx(A11+1,B11+2) = Kud(A12+1,j)            
            !amatrx(A11+2,B11)   = Kud(nDim*(j-1)+1,i)
            !amatrx(A11+2,B11+1) = Kud(nDim*(j-1)+1+1,i)
            !
            ! order parameter
            !
            amatrx(A11+nDim+1,B11+nDim+1) = Kdd(i,j)
         enddo
      enddo
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------

      return
      end subroutine uel
!**************************************************************************

!     check the lflags array to make sure the UEL is being used with
!     appropriate options.
!      
      SUBROUTINE KCheckLflags(lflags)

      implicit none
      integer :: lflags
      dimension lflags(*)

      ! Check procedure type
      if ( (lflags(1).ne.71) .and. (lflags(1).ne.72) .and. 
     &     (lflags(1).ne.73) ) then
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and check the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         call xit
      endif

      ! Make sure Abaqus knows you are doing a large
      !  deformation problem
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a finite displacement analysis'
         write(*,*) 'go in and set nlgeom=no'
         call xit
      endif

      ! Ensure that we are performing a general step,
      ! not a linear purturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear purturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear purturbation step'
         call xit         
      endif

      return
      end subroutine KCheckLflags
****************************************************************************
!     Dummy material
!     Returns zero stress and Jacobians.
!     Material point is deleted when the phase field
!     reaches a user set value.
!
!     The material point deletion is purely for visualization.
!     The UEL continues to exist and to contribute force and
!     stiffness.
!     
      subroutine UMAT(stress, statev, ddsdde, sse, spd, scd,
     &     rpl, ddsddt, drplde, drpldt,
     &     stran, dstran, time, dtime, temp, dtemp, predef, dpred,
     &     cmname, ndi, nshr, ntens, nstatv, props, nprops, coords,
     &     drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer,
     &     kspt, jstep, kinc)

      use global

#if !TEST_MODE
      include 'ABA_PARAM.INC'
#endif
c      implicit none
      
      character(len=80) cmname

!     Formal arguments to function
      
      real(8) :: stress(ntens), statev(nstatv), ddsdde(ntens, ntens)
      real(8) :: ddsddt(ntens), drplde(ntens), drpldt, stran(ntens)
      real(8) :: dstran(ntens), time(2), dtime, predef(1), dpred(1)
      real(8) :: props(nprops), coords(3), drot(3,3), dfgrd0(3,3)
      real(8) :: dfgrd1(3,3), sse, spd, scd, rpl, temp, dtemp, pnewdt
      real(8) :: celent

      integer :: ndi, nshr, ntens, nstatv, nprops, noel, npt, layer
      integer :: kspt, jstep(4), kinc

!     local variables
      real(8) :: phi, phi_fail


!     Get phase field damage variable at this location
!     from global array.
!     The phase field is in slot 1
      phi = 0
      if (allocated(globalSdv)) then
         phi = globalSdv(1, npt, noel - elOffset)
      end if

!     The only material property is the value
!     of the phase field at which to set the
!     element deletion flag.      
      phi_fail = props(1)

!     Fail this material point if the damage variable
!     exceeds the threshold.
      if (phi .ge. phi_fail) then
         statev(1) = 0
      end if

      stress = 0
      ddsdde = 0
 
      return
      end subroutine UMAT
****************************************************************************
****************************************************************************
!     
!     Compressible Gent model impelemented with UHYPER
!     Not used in UEL, but useful for verification and testing
!      
      subroutine UHYPER(BI1, BI2, AJ, U, UI1, UI2, UI3, temp, noel,
     &     cmname, incmpflag, numstatev, statev, numfieldv, fieldv,
     &     fieldvinc, numprops, props)

#if !TEST_MODE
      include 'ABA_PARAM.INC'
#endif

      character(len=80) :: cmname
      dimension U(2), UI1(3), UI2(6), UI3(6), statev(*), fieldv(*),
     &     fieldvinc(*), props(*)

      real(8) :: Gshear, Kbulk, Im, Gbar, logJ

!     give meaningful names to mat properties
      Gshear   = props(1)
      Kbulk    = props(2)
      Im       = props(3)

      
!     Constitutive reponse functions
      U(1) = -0.5d0*Gshear*Im*dlog(1.d0 - (BI1 - 3.d0)/Im) + 
     &     0.5d0*Kbulk*(AJ - 1.d0)**2                        ! free energy density in reference config 
      U(2) = -0.5d0*Gshear*Im*dlog(1.d0 - (BI1 - 3.d0)/Im)   ! deviatoric part of strain energy, not really needed

      Gbar = Gshear*Im/(Im - (BI1 - 3.d0))
      
      UI1 = 0
      UI1(1) = 0.5d0*Gbar        ! dpsiR/dI1
      UI1(3) = Kbulk*(AJ - 1.d0) ! dpsiR/dJ

      UI2 = 0
      UI2(1) = 0.5d0*Gbar/(Im - (BI1 - 3.d0)) ! d^2 psiR/d I1^2
      UI2(3) = Kbulk                          ! d^2 psiR/d J^2

      UI3 = 0
c      UI3(6) = Kbulk/AJ**3*(2.d0*logJ - 3.d0) ! d^3 psiR/d J^3

      end subroutine UHYPER
****************************************************************************
      subroutine kprtmat(A,m,n,unit)

      implicit none
      
      integer :: m,n,k,l,unit
      real(8) :: A(m,n) 

      do k=1,m
         write(unit,*) (A(k,l), l=1,n)
      enddo
      
      return
      end
****************************************************************************
      subroutine QuadratureRule(nDim, nNode, nInt, xi, w)

      implicit none
      
      integer :: nDim, nNode, nInt, nIntPt
      real(8) :: xi(nInt, nDim), w(nInt)
      
         
      if((nNode.eq.4).or.(nNode.eq.8).or.(nNode.eq.9)) then
!        
!        gauss integration for a rectangular element
!
         select case (nInt)
         case (1)
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         case (4)
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         case (9)
            call xint2D9pt(xi,w,nIntPt) ! 9-pt integration, nInt=9 above
         case DEFAULT
            write(*,*) 'Error in UEL'
            write(*,*) 'Illegal integration rule requested, nInt=',
     &           nInt
            call Xit
         end select
      elseif((nNode.eq.3).or.(nNode.eq.6)) then
!
!        gauss integration for a triangular element
!
         select case (nInt)
         case (1)
            call xintTri1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         case (3)
            call xintTri3pt(xi,w,nIntPt) ! 3-pt integration, nInt=3 above
         case (4)
            call xintTri4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         case (6)
            call xintTri6pt(xi,w,nIntPt) ! 6-pt integration, nInt=6 above
         case DEFAULT
            write(*,*) 'Error in UEL:'
            write(*,*) 'Illegal integration rule requested, nInt=',
     &           nInt
            call Xit
         end select
      else
         write(*,*) 'Unknown element specification:'
         write(*,*) 'nNode=',nNode
         call xit
      endif


      end subroutine QuadratureRule
****************************************************************************
      subroutine CalcShape(nDim, nNode, nIntPt, intpt, xi, sh, dshxi)

      implicit none

      integer :: nDim, nNode, nIntPt, intpt
      real(8) :: xi(nIntPt, nDim), sh(nNode), dshxi(nNode, nDim)

      
      select case (nNode)
      case (3)
         call calcShapeTriLinear(nIntPt,xi,intpt,sh,dshxi)
      case (4)
         call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
      case (6)
         call calcShapeTriQuad(nIntPt,xi,intpt,sh,dshxi)
      case (8)
         call calcShape2DQuad(nIntPt,xi,intpt,sh,dshxi)
      case (9)
         call calcShape2DQuadFull(nIntPt,xi,intpt,sh,dshxi)
      case default
         write(*,*)
     &        'Illegal value for number of nodes per element'
         write(*,*) 'nNode=',nNode
         call Xit
      end select
         
      end subroutine CalcShape
****************************************************************************
