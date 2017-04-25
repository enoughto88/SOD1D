module parameters_mod

   implicit none
   !Program parameter
   character(len=80)             :: TOPDIR= '/home/kawashima/SOD1D'     !Top directory in the Linux Workstation
   integer,parameter             :: ITMAX = 200                         !Maximum number of time steps
   integer,parameter             :: ISMP  = 10                          !Sampling interval
   double precision,parameter    :: DT    = 1.0d-3                      !Time step inverval
   !Grid parameter
   integer,parameter             :: NXMAX = 100                         !Number of cells
   double precision,parameter    :: XL    = 1.2d0                       !Length of calculation field in x-direction
   double precision,parameter    :: DX    = XL/dble(NXMAX)              !Length of cell in x-direction
   !Physical parameter
   double precision,parameter    :: GAMMA = 1.4d0                       !Specific heat ratio

end module parameters_mod


