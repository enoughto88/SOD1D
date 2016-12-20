!***********************************************************************
!*****                  1-D Euler equation solver                  *****
!*****         Developed by Rei Kawashima (Univ. of Tokyo)         *****
!*****                   Last Update 12.14.2016                    *****
!***********************************************************************

program main
   use parameters_mod
   implicit none
   character(len=6)  :: cname
   character(len=30) :: dname
   integer :: nx
   integer :: it                                                        !Time step
   integer :: nstp = 0                                                  !Sampling step
   double precision,dimension(1:NXMAX)      :: rho                      !Density
   double precision,dimension(1:NXMAX)      :: vel                      !Velocity
   double precision,dimension(1:NXMAX)      :: pre                      !Pressure
   double precision :: cfl                                              !CFL number

   !Program start processes
   write(*,*) 'Program start...'

   !Initial condition setting, SOD shock-tube problem
   do nx = 1,NXMAX
      if(nx.le.int(dble(NXMAX)/2.0d0)) then
         rho(nx) = 1.0d0
         vel(nx) = 0.0d0
         pre(nx) = 1.0d0
      else
         rho(nx) = 0.125d0
         vel(nx) = 0.0d0
         pre(nx) = 0.1d0
      endif
   enddo

   MARCH: do it = 1,ITMAX
      if(mod(it,ISMP).eq.0) nstp = nstp+1

      !Electron fluid solver
      call Euler1D(NXMAX,DX,DT,GAMMA,rho,vel,pre,cfl)

      !èOutput
      if(mod(it,ISMP).eq.0) then
         !Information for monitoring
         write(*,*) '***********************************************'
         write(*,'(A10,I10,A)')   '          ',nstp,'-th Step'
         write(*,'(A20,E12.3,A)') '       CFL number =',cfl             ,' [s]'

         !Distribution output
         write(cname,'(I6.6)') nstp
         dname = '/output/distribution/'
         open(unit=20,file=trim(TOPDIR)//trim(dname)//'distribution.'//&
            cname//'.dat',form='formatted',status='unknown')
            do nx = 1,NXMAX
               write(20,'(3E15.5)') &
                  rho(nx),vel(nx),pre(nx)
            enddo
         close(20)

      endif
   enddo MARCH

   !Program end processes
   write(*,*) 'Program end...'


   stop
end program main




