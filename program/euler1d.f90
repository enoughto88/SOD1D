!***********************************************************************
!*****                 1-D Euler Equation Solver                   *****
!*****         Developed by Rei Kawashima (Univ. of Tokyo)         *****
!*****                  Last Update 12.16.2016                     *****
!***********************************************************************

subroutine Euler1D(NXMAX,DX,DT,GAMMA,rho,vel,pre,cfl)

   implicit none
   integer                            ,intent(in)    :: NXMAX           !Number of cells
   double precision                   ,intent(in)    :: DX              !Grid sizing
   double precision                   ,intent(in)    :: DT              !Time step interval
   double precision                   ,intent(in)    :: GAMMA           !Specific heat ratio
   double precision,dimension(1:NXMAX),intent(inout) :: rho             !Density
   double precision,dimension(1:NXMAX),intent(inout) :: vel             !Velocity
   double precision,dimension(1:NXMAX),intent(inout) :: pre             !Pressure
   double precision                   ,intent(out)   :: cfl             !CFL number
   double precision,dimension(-1:NXMAX+2,3)          :: uvar            !Primitive variables inclusing ghost cells
   double precision,dimension(1:NXMAX,3)             :: qvar            !Conservative variables
   double precision,dimension(1:NXMAX+1,3,2)         :: umuscl          !MUSCL interpolated primitive variables, defined at cell boundary
   double precision,dimension(1:NXMAX+1,3)           :: flux            !Numerical flux, defined at cell boundary
   double precision,dimension(1:NXMAX,3)             :: rhs             !Right hand side
   double precision,dimension(1:NXMAX,3)             :: delta           !Delta of conservative variables
   integer :: i,k
   double precision :: MaxCharSpeed


   !Primitive variable vector
   do i = 1,NXMAX
      uvar(i,1)  = rho(i)
      uvar(i,2)  = vel(i)
      uvar(i,3)  = pre(i)
   enddo

   !Boundary condition, SOD shock-tube problem, 1st-order
   uvar(0   ,1)    = 2.0d0*1.0d0          -uvar(1    ,1)
   uvar(0   ,2)    = 2.0d0*0.0d0          -uvar(1    ,2)
   uvar(0   ,3)    = 2.0d0*1.0d0          -uvar(1    ,3)
   uvar(NXMAX+1,1) = 2.0d0*0.125d0        -uvar(NXMAX,1)
   uvar(NXMAX+1,2) = 2.0d0*0.0d0          -uvar(NXMAX,2)
   uvar(NXMAX+1,3) = 2.0d0*0.1d0          -uvar(NXMAX,3)
   uvar(-1  ,1)    = 2.0d0*uvar(0   ,1)   -uvar(1    ,1)
   uvar(-1  ,2)    = 2.0d0*uvar(0   ,2)   -uvar(1    ,2)
   uvar(-1  ,3)    = 2.0d0*uvar(0   ,3)   -uvar(1    ,3)
   uvar(NXMAX+2,1) = 2.0d0*uvar(NXMAX+1,1)-uvar(NXMAX,1)
   uvar(NXMAX+2,2) = 2.0d0*uvar(NXMAX+1,2)-uvar(NXMAX,2)
   uvar(NXMAX+2,3) = 2.0d0*uvar(NXMAX+1,3)-uvar(NXMAX,3)

   !Conservative variable vector
   do i = 1,NXMAX
      qvar(i,1) = uvar(i,1)
      qvar(i,2) = uvar(i,1)*uvar(i,2)
      qvar(i,3) = uvar(i,3)/(GAMMA-1.0d0)+0.5d0*uvar(i,1)*uvar(i,2)*uvar(i,2)
   enddo

   !MUSCL interoplated primitive variable vector
   call MUSCL1D(NXMAX,3,uvar,umuscl,2)  !Second order

   !Flux calculated by AUSM scheme
   call AUSM(NXMAX,GAMMA,umuscl,flux)

   !CFL number
   MaxCharSpeed = 0.0d0
   do i = 1,NXMAX
      MaxCharSpeed = dmax1(MaxCharSpeed,dabs(vel(i))+dsqrt(GAMMA*pre(i)/rho(i)))
      cfl = MaxCharSpeed*DT/DX
   enddo

   !Right-hand side calculation
   do i = 1,NXMAX
      do k = 1,3
         rhs(i,k) =-DT/DX*(flux(i+1,k)-flux(i,k))
      enddo
   enddo

   !Explicit method
   do i = 1,NXMAX
      do k = 1,3
         delta(i,k) = rhs(i,k)
      enddo
   enddo

   !Variable update
   do i = 1,NXMAX
      do k = 1,3
         qvar(i,k) = qvar(i,k)+delta(i,k)
      enddo
   enddo

   !Calculation of elemental quantities
   do i = 1,NXMAX
      rho(i) = qvar(i,1)
      vel(i) = qvar(i,2)/qvar(i,1)
      pre(i) = (qvar(i,3)-0.5d0*qvar(i,2)*qvar(i,2)/qvar(i,1))*(GAMMA-1.0d0)
   enddo


   return
endsubroutine


!***********************************************************************
!*****                        AUSM Scheme                          *****
!*****        Ref.: M.-S. Liou, J. of Comp. Phys. 129 (1996)       *****
!*****                  Last Update 2016.12.08                     *****
!***********************************************************************

subroutine AUSM(NXMAX,GAMMA,umuscl,flux)

   implicit none
   integer                                  ,intent(in)  :: NXMAX       !Number of cells
   double precision                         ,intent(in)  :: GAMMA       !Specific heat ratio
   double precision,dimension(1:NXMAX+1,3,2),intent(in)  :: umuscl      !MUSCL interpolated primitive variables, defined at cell boundary
   double precision,dimension(1:NXMAX+1,3)  ,intent(out) :: flux        !Numerical flux, defined at cell boundary
   double precision,dimension(1:NXMAX+1,3,2)             :: pmuscl      !
   double precision,dimension(1:NXMAX+1,2)               :: aborder     !
   double precision,dimension(1:NXMAX+1,2)               :: Mach        !
   double precision,dimension(1:NXMAX+1)                 :: mborder     !
   double precision,dimension(1:NXMAX+1)                 :: Mplus,Mminus!
   double precision,dimension(1:NXMAX+1,2)               :: machpm      !Mach number, defined at cell boundary
   double precision,dimension(1:NXMAX+1,2)               :: prespm      !Pressure coef., defined at cell boundary
   integer :: i,k

   !Calculation of speed of sound and Mach number
   do i = 1,NXMAX+1
      aborder(i,1)= dsqrt(GAMMA*umuscl(i,3,1)/umuscl(i,1,1))
      aborder(i,2)= dsqrt(GAMMA*umuscl(i,3,2)/umuscl(i,1,2))
      Mach(i,1)   = umuscl(i,2,1)/aborder(i,1)
      Mach(i,2)   = umuscl(i,2,2)/aborder(i,2)
   enddo

   !Calculation of m^{+-}_{i+1/2}
   do i = 1,NXMAX+1
      if(dabs(Mach(i,1)).ge.1.0d0) then
         Mplus(i)  = 0.5d0*(Mach(i,1)+dabs(Mach(i,1)))
      else
         Mplus(i)  = 0.25d0*(Mach(i,1)+1.0d0)**2.0d0
      endif
      if(dabs(Mach(i,2)).ge.1.0d0) then
         Mminus(i) = 0.5d0*(Mach(i,2)-dabs(Mach(i,2)))
      else
         Mminus(i) =-0.25d0*(Mach(i,2)-1.0d0)**2.0d0
      endif
      mborder(i)  = Mplus(i)+Mminus(i)
      machpm(i,1) = 0.5d0*(mborder(i)+dabs(mborder(i)))
      machpm(i,2) = 0.5d0*(mborder(i)-dabs(mborder(i)))
   enddo

   !Calculation of pressure at cell boundaries
   do i = 1,NXMAX+1
      if(dabs(Mach(i,1)).ge.1.0d0) then
         prespm(i,1) = 0.5d0*(1.0d0+dsign(1.0d0,Mach(i,1)))
      else
         prespm(i,1) = 0.25d0*(Mach(i,1)+1.0d0)**2.0d0*(2.0d0-Mach(i,1))
      endif
      if(dabs(Mach(i,2)).ge.1.0d0) then
         prespm(i,2) = 0.5d0*(1.0d0-dsign(1.0d0,Mach(i,2)))
      else
         prespm(i,2) = 0.25d0*(Mach(i,2)-1.0d0)**2.0d0*(2.0d0+Mach(i,2))
      endif
   enddo

   !MUSCL interoplated flux vector
   do i = 1,NXMAX+1
      do k = 1,2
         pmuscl(i,1,k) = umuscl(i,1,k)
         pmuscl(i,2,k) = umuscl(i,1,k)*umuscl(i,2,k)
         pmuscl(i,3,k) = GAMMA/(GAMMA-1.0d0)*umuscl(i,3,k)+0.5d0*umuscl(i,1,k)*umuscl(i,2,k)**2.0d0
      enddo
   enddo

   !Flux
   do i = 1,NXMAX+1
      flux(i,1) = aborder(i,1)*machpm(i,1)*pmuscl(i,1,1)&
                 +aborder(i,2)*machpm(i,2)*pmuscl(i,1,2)
      flux(i,2) = aborder(i,1)*machpm(i,1)*pmuscl(i,2,1)&
                 +aborder(i,2)*machpm(i,2)*pmuscl(i,2,2)&
                 +             prespm(i,1)*umuscl(i,3,1)&
                 +             prespm(i,2)*umuscl(i,3,2)
      flux(i,3) = aborder(i,1)*machpm(i,1)*pmuscl(i,3,1)&
                 +aborder(i,2)*machpm(i,2)*pmuscl(i,3,2)
   enddo


   return
endsubroutine



!***********************************************************************
!*****                     TVD-MUSCL Method                        *****
!*****         Developed by Rei Kawashima (Univ. of Tokyo)         *****
!*****                  Last Update 2016.12.18                     *****
!***********************************************************************

subroutine MUSCL1D(NXMAX,NVAR,uvar,umuscl,order)

   implicit none
   integer                                     ,intent(in)   :: NXMAX   !Number of cells
   integer                                     ,intent(in)   :: NVAR    !Number of variables
   double precision,dimension(-1:NXMAX+2,NVAR) ,intent(in)   :: uvar    !Primitive variables inclusing ghost cells
   double precision,dimension(1:NXMAX+1,NVAR,2),intent(out)  :: umuscl  !MUSCL interpolated primitive variables, defined at cell boundary
   integer                                     ,intent(in)   :: order   !Order of accuracy
   integer          :: i,k
   double precision :: limitdm,limitdp,dp,dm,bb,kk,aa

   ! o: cell center, |: cell boundary, *: interpolated value
   !                |       |       |       |       |
   !                |   o   |   o  *|*  o   |   o   |
   !                |       |       |       |       |
   !Cell index         i-2     i-1      i      i+1
   !Boundary index         i-1      i      i+1


   !MUSCL interpolation using minmod limiter (k = 1/3)
   if(order.eq.1) then
      aa = 0.0d0
      kk = 0.0d0
      bb = 0.0d0
   elseif(order.eq.2) then
      aa = 1.0d0
      kk =-1.0d0
      bb = 1.0d0
   elseif(order.eq.3) then
      aa = 1.0d0
      kk = 1.0d0/3.0d0
      bb = (3.0d0-kk)/(1.0d0-kk)
   endif
   do i = 1,NXMAX+1
      do k = 1,NVAR
         dp = uvar(i+1,k)-uvar(i  ,k)
         dm = uvar(i-1,k)-uvar(i-2,k)
         !Minmod limiter function is used
         limitdp = dsign(1.0d0,dp)*dmax1(0.0d0,dmin1(dabs(dp),dsign(1.0d0,dp)*bb*dm))
         limitdm = dsign(1.0d0,dm)*dmax1(0.0d0,dmin1(dabs(dm),dsign(1.0d0,dm)*bb*dp))
         umuscl(i,k,1) = uvar(i-1,k)+aa*0.25d0*((1.0d0-kk)*limitdm+(1.0d0+kk)*limitdp)
         umuscl(i,k,2) = uvar(i  ,k)-aa*0.25d0*((1.0d0-kk)*limitdp+(1.0d0+kk)*limitdm)
      enddo
   enddo


   return
endsubroutine






