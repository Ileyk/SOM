module preliminaries
use miscellaneous
contains

! -----------------------------------------------------------------------------------
subroutine cmpt_params
include 'def.f90'

! Compute initial total amount of angular momentum
Jtot_ini = 1.d0 + (1.d0/5.d0)*om0_*R_**2.d0*(1.d0+q_)
! Minimum total ang. mom. to be able to synchronize
Jc   = ((4.d0/5.d0)/(3.d0/5.d0)**0.75d0)*R_**0.5d0*(1.d0+q_)**0.25d0
! Omega synch. when Jtot=Jc (=> om1=om2), to pick up 1st guess for om1 and om2
OmC = ((3.d0/5.d0)*(1.d0+q_)*R_**2.d0)**(-3.d0/4.d0)

end
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine prm_chcks
include 'def.f90'
double precision :: om1=0.d0, om2=0.d0

open(9,file=prm_file)

! 1.
! Is the initial total amount of angular momentum above the minimum value
! necessary to make synchronization possible?
if (Jtot_ini>Jc) write(9,*) 'Jtot initial > J min for synch => synchronization possible', Jtot_ini, Jc
if (Jtot_ini<Jc) write(9,*) 'Jtot initial < J min for synch => doomed to merge', Jtot_ini, Jc

! 2.
! If synch possible, compute the 2 omegas possible
if (Jtot_ini>Jc) then
  ! Given Jtot_ini>Jc, which roots om1 & om2? w/ om1 < om2
  call get_roots(Jtot_ini,om1,om2) ! test : for R=0.9 & q=5, Jtot_ini~Jc & omC~0.76375d0
endif

! 3.
! Given the initial conditions, where are we heading towards? Synch or merge?
if (Jtot_ini>Jc) then
  write(9,*) 'om1 =', om1
  write(9,*) 'om2 =', om2
  if (1.d0<om1) write(9,*) 'Orb. sep. shrinking to synchronization : Om.orb.0 < om1 < om2'
  if (om1<1.d0 .and. 1.d0<om2) write(9,*) 'Orb. sep. increasing to synchronization: om1 < Om.orb.0 < om2'
  if (om2<1.d0) write(9,*) 'Orb. sep. shrinking until merger: om1 < om2 < Om.orb.0'
else
  write(9,*) 'Doomed to merge'
endif

! 4.
! If the star had an important convective envelope (which it should not have),
! would dynamical tides arise?
if (1.d0<2.d0*om0_) then
  write(9,*) 'BEWARE : if the env. was conv., dyn. tides would be significant.'
endif

close(9)

end subroutine prm_chcks
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
subroutine get_roots(Jtot,om1,om2)
include 'def.f90'
double precision, intent(in) :: Jtot
double precision, intent(out) :: om1, om2
double precision :: Jsynch

double precision :: ep, p, dp, omX
integer :: i=0, j, maxi

maxi=50

ep=1.0d-10

do j=1,2

  ! First guess
  if (j==1) omX=OmC/2.d0 ! find om1 (<OmC)
  if (j==2) omX=OmC*2.d0 ! find om2 (>OmC)

LoopNR:  do i=1,maxi

   p  = (1.d0/5.d0)*(1.d0+q_)*R_**2.d0*omX + omX**(-1.d0/3.d0) - Jtot

   dp = (1.d0/5.d0)*(1.d0+q_)*R_**2.d0     - (1.d0/3.d0)*omX**(-4.d0/3.d0)

   if (omX-p/dp>0.d0) then ! always the case beyond omC, right?
     omX=omX-p/dp
   else ! to prevent an iteration of om1 to go below 0.
     omX=omX/2.d0
   endif
   p  = (1.d0/5.d0)*(1.d0+q_)*R_**2.d0*omX + omX**(-1.d0/3.d0) - Jtot

   ! Convergence condition
   if (dabs(p)<ep) then
      if (j==1) om1=OmX ! find om1 (<OmC)
      if (j==2) om2=OmX ! find om2 (>OmC)
      exit LoopNR
   endif

enddo LoopNR

if (i>maxi) then
   print*, i, j, Jtot, dabs(p)
   print*, 'Issue 2'
   stop
endif

enddo

end subroutine get_roots
! -----------------------------------------------------------------------------------

end module preliminaries
