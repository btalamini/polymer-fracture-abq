program TestNewton
  ! this program compares the derivative in the constitutive update
  ! with finite difference estimates

  implicit none

  integer, parameter :: nint=4, nlsdv=1, ngsdv=6

  real(8) Eyoung, anu, Gc, len
  
  integer, parameter :: ndofel=12, nrhs=1, nsvars=nint*nlsdv, nprops=4, mcrd=2, nnode=4, &
       jtype=1, kstep=1, kinc=1, jelem=1, ndload=1, npredef=1, mlvarx=1, mdload=1, &
       njprop=3, nparams=1

  real(8) :: rhs(mlvarx,ndofel), amatrx(ndofel,ndofel), svars(nsvars), energy(8), &
       props(nprops), coords(mcrd,nnode), UAll(ndofel), DuAll(mlvarx,ndofel), Vel(ndofel), &
       Accn(ndofel), time(2), dtime, params(nparams), adlmag(mdload,1), &
       predef(2,npredef,nnode), ddlmag(mdload,1), pnewdt, period

  integer :: jdltyp(mdload,1), lflags(5), jprops(njprop)

  real(8) :: rp(ndofel), rm(ndofel), Kh(ndofel,ndofel), Kexact(ndofel,ndofel), du, &
       max_error, error

  integer i,j,imax,jmax

  Eyoung = 218d6
  anu    = 0.3d0
  Gc     = 20d3
  len    = 4d-2

  ! set up material properties
  props(1) = Eyoung
  props(2) = anu
  props(3) = Gc
  props(4) = len

  jprops(1) = nlsdv
  jprops(2) = ngsdv
  jprops(3) = 1

  ! set coordinates for element to be unit square
  coords(1,1) = 0.0d0
  coords(2,1) = 0.0d0

  coords(1,2) = 1.0d0
  coords(2,2) = 0.0d0

  coords(1,3) = 1.0d0
  coords(2,3) = 1.0d0

  coords(1,4) = 0.0d0
  coords(2,4) = 1.0d0

  ! initialize finite element arrays
  rhs = 0
  amatrx = 0
  svars = 0

  UAll  = 0
  DuAll = 0
  Vel   = 0
  Accn  = 0

  ! time not used in this material
  time(1) = 1.0d0
  time(2) = 1.0d0
  dtime = 1.0d0

  ! Abaqus solution parameter flags
  lflags(1) = 1
  lflags(2) = 0
  lflags(3) = 1
  lflags(4) = 0
  lflags(5) = 0


  !
  ! test analytical derivative of root equation
  !
  
  !Uall(4) = 5d-3
  !Uall(7) = 5d-3
  !Uall(8) = -anu/(1.0d0-anu)*Uall(7)
  !Uall(11) = -anu/(1.0d0-anu)*Uall(7)

  !Uall(3) = 9.829d-2
  !Uall(6) = 9.829d-2
  !Uall(9) = 9.829d-2
  !Uall(12) = 9.829d-2

  Kexact=0
  
  call UEL(rhs,Kexact,svars,energy,ndofel,nrhs,nsvars, &
       props,nprops,coords,mcrd,nnode,Uall,Duall,  &
       Vel,Accn,jtype,time,dtime,kstep,kinc,jelem, &
       params,ndload,jdltyp,adlmag,predef,npredef, &
       lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop, &
       period)
  
  rp=0
  rm=0
  du = 1.d-8

  do i=1,ndofel
     Uall(i) = Uall(i) + du
     
     call UEL(rp,amatrx,svars,energy,ndofel,nrhs,nsvars, &
          props,nprops,coords,mcrd,nnode,Uall,Duall,  &
          Vel,Accn,jtype,time,dtime,kstep,kinc,jelem, &
          params,ndload,jdltyp,adlmag,predef,npredef, &
          lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop, &
          period)
     
     Uall(i) = Uall(i) - 2.d0*du
     
     call UEL(rm,amatrx,svars,energy,ndofel,nrhs,nsvars, &
          props,nprops,coords,mcrd,nnode,Uall,Duall,  &
          Vel,Accn,jtype,time,dtime,kstep,kinc,jelem, &
          params,ndload,jdltyp,adlmag,predef,npredef, &
          lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop, &
          period)
     
     Uall(i) = Uall(i) + du
     
     do j=1,ndofel
        Kh(i,j)=-(rp(j)-rm(j))/(2.0d0*du)
     end do
     
  end do

  ! find entry with worst error
  max_error=0
  imax = 0
  jmax = 0
  do i=1,ndofel
     do j=1,ndofel
        error = DABS(Kexact(i,j)-Kh(i,j))
        if (error .gt. max_error) then
           max_error = error
           imax = i
           jmax = j
        end if
     end do
  end do
  
  
  write(*,*) 'Max error: ',max_error
  if (max_error .gt. 0.0d0) then
     write(*,*) 'Entry (', imax, ',', jmax, ')'
     write(*,*) 'Kh(', imax, ',', jmax, ') = ', Kh(imax,jmax)
     write(*,*) 'Kexact(', imax, ',', jmax, ') = ', Kexact(imax,jmax)
  end if

end program TestNewton


subroutine Xit
  write(*,*) 'subroutine Xit() called'
end subroutine Xit
