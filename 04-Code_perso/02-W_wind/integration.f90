module integration
use miscellaneous
contains

! -----------------------------------------------------------------------------------
subroutine getdt(X,qdt)
include 'def.f90'
double precision, intent(in) :: X
double precision, intent(out) :: qdt
double precision :: Xdot

call derivs(X,Xdot)

qdt = CFL * (X/dabs(Xdot))

print*, qdt, X, Xdot

if (ISNAN(qdt)) then
   print*, 'qdt_cfl', qdt!, dx(1), dw(1)!v, a
   stop
endif

! double precision :: dx(1:3), dw(1:nw), v, a
! double precision :: factor, r2
! call derivs(x,w,dx,dw,with_compact_object)
!
! ! r = dsqrt(x(1)**two+x(2)**two+x(3)**two)
! v = dsqrt(w(1)**two+w(2)**two+w(3)**two)
! a = dsqrt(dw(1)**two+dw(2)**two+dw(3)**two)

! f is a function designed such as the time step drops when one gets close
! from the "OK sphere" of radius rout_ around the compact object.
! It enables us to get arrival radii similar for the streamlines selected
! and post-processed afterwhile. Otherwise, there is a "step effect"
! with discrete shifts due to the fact that some streamlines have additional time steps.
! f goes from 1 to 1E-2 for r2/rout (r2 the distance to the CO) going from
! a (and beyond) to 1.
! The risk : it slows down the non selected streamlines which pass around...
! r2 = dsqrt(x(1)**two+(x(2)+eqpar(a_))**two+x(3)**two)
! ! Linear profile
! !f = ( (1.d0-1.d-2)/(eqpar(a_)/eqpar(rout_)-one) ) * (r2/eqpar(rout_)) + 1.d-2
! ! Exponential profile such as f stays ~ @ 1 except for r/rout below 1.5 where it exponentially drops
! factor = 1.d-2 + (1.d0-1.d-2) * ( one - dexp(- (r2/eqpar(rout_)-one) / (1.5d0-one)) )
!print*, f
! The use of eqpar(roche2_) rather than eqpar(a_) is to control the radial shifts @ arrival
! (@ fixed courantpar). Indeed, otherwise, some streamlines have an/few additional time step(s)
! compared to other which results in shifts in arrival theta and phi
! (can be seen with surf_pts :
! p '/Users/ielm/Desktop/NS_SgXB_ballistic_v5/output/pts' u 3:2 ----> arrival pts from
! p '/Users/ielm/Desktop/NS_SgXB_ballistic_v5/output/pts' u 12:11 --> those departure ones)
! qdt_cfl = min(courantpar*      dabs(eqpar(NStoL1_)/v),&
!                      courantpar*dsqrt(dabs(eqpar(NStoL1_)/a)),&
!                      courantpar*      dabs((eqpar(NStoL1_)**two)/(v*a))**(1.d0/3.d0))
! ! print*, qdt_cfl!, 'sf', w, dw
!
! if (ISNAN(qdt_cfl)) then
!    print*, 'qdt_cfl', qdt_cfl!, dx(1), dw(1)!v, a
!    stop
! endif

! if (qdt_cfl<(one/eqpar(Om_))/1.d8) call crash('Error#06 : time step too small') !X!

! print*, courantpar*      dabs(eqpar(a_)/v), courantpar*dsqrt(dabs(eqpar(a_)/a)), &
! courantpar*      dabs((eqpar(a_)**two)/(v*a))**(1.d0/3.d0)

! print*, qdt_cfl, courantpar * &
!   ((eqpar(a_)**two)/dsqrt((w(1)**two+w(2)**two+w(3)**two)*(dw(1)**two+dw(2)**two+dw(3)**two)))**(2.d0/3.d0)

! qdt = dtparam * &
!   ((eqpar(a_)**two)/dsqrt((w(1)**two+w(2)**two+w(3)**two)*(dw(1)**two+dw(2)**two+dw(3)**two)))**(2.d0/3.d0)

end subroutine getdt
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
subroutine derivs(X,Xdot)
include 'def.f90'
double precision, intent(in) :: X
double precision, intent(out) :: Xdot
double precision :: tau_eff

tau_eff = tau_ * (2.d0/3.d0) * (1.d0/2.d0) * q_

Xdot = - (1.d0/5.d0) * (R_**8.d0/tau_eff) * X**(-12.d0) * ( X**(-3.d0) - om0_ + (5.d0/(R_**2.d0*(1.d0+q_)))*(X-1.d0) )

! double precision :: xtmp(1:3)
! double precision :: gl, grav1, grav2, cent, cor, r1, r2, r
! integer :: i
! double precision :: r1x, r1y, r1z, r2x, r2y, r2z, rx, ry, rz, &
!                     grav1x, grav1y, grav1z, grav2x, grav2y, grav2z, centx, centy, centz, corx, cory, corz, glx, gly, glz
! double precision :: rad_andreas_tmp
!
! ! x_dot = w : velocities
! do i=1,3
!    dx(i)=w(i)
! enddo
!
! ! This force is computed in the frame centered on the star
! !call compute_glines_from_betalaw(x,gl)
! ! print*, x(1)-xprev(1)
! ! print*, w(1)-wprev(1)
!
! ! call compute_glines_realistic(x,xprev,w,wprev,gl)
! ! print*, x(1)**two+x(2)**two+x(3)**two, gl
! ! print*, 'd',gl * (x(2)/dsqrt(x(1)**two+x(2)**two+x(3)**two)), - eqpar(GM1mod_) * ( x(2) / ( x(1)**two+x(2)**two+x(3)**two)**1.5d0 )
!
! ! Beware! The following forces are computed in the frame centered on the center of mass (using xtmp)! In agreement w/ notes
! xtmp=x
! ! xtmp(1)=xtmp(1)-eqpar(CM_)
! xtmp(2)=xtmp(2)+eqpar(CM_)
!
! ! acceleration_x
! ! r1 = dsqrt((xtmp(1)+(one/eqpar(egg_))*(one      /(one+shpe(q_))))**two+xtmp(2)**two+xtmp(3)**two)
! ! r  = dsqrt(xtmp(1)**two                                                +xtmp(2)**two+zero)       ! projected distance
! ! r2 = dsqrt((xtmp(1)-(one/eqpar(egg_))*(shpe(q_)/(one+shpe(q_))))**two+xtmp(2)**two+xtmp(3)**two)
!
! r1x = xtmp(1) !+(one/eqpar(egg_))*(one/(one+shpe(q_)))
! r1y = xtmp(2) - eqpar(CM_) ! - eqpar(a_) * (one     /(one+shpe(q_)))
! r1z = xtmp(3)
! r2x = xtmp(1) !-(one/eqpar(egg_))*(shpe(q_)/(one+shpe(q_)))
! r2y = xtmp(2) + (one/(shpe(fll_)*eqpar(egg_))) * (shpe(q_)/(one+shpe(q_))) ! + eqpar(a_) * (shpe(q_)/(one+shpe(q_)))
! r2z = xtmp(3)
! rx  = xtmp(1)
! ry  = xtmp(2)
! rz  = zero
!
! r1 = dsqrt(r1x**two+r1y**two+r1z**two)
! r2 = dsqrt(r2x**two+r2y**two+r2z**two)
! r  = dsqrt(rx **two+ry **two+rz **two)
! !print*, gl, (shpe(q_)*(one-shpe(gmEdd_))) / r1**2.d0
!
! if (law .EQ. 'Andr') then
!   ! Interpolate rad_andreas based on the distance to the star
!   call interp_rad_andreas(r1,rad_andreas_tmp)
! else
!   rad_andreas_tmp=(shpe(bet_)*(one-one/r1)**shpe(bet_))/(r1**(one+shpe(bet_))*(r1-one)**(one-shpe(bet_)))
! endif
! ! print*, r1, r2, 'h'
! ! print*, x(1), x(2), x(3), eqpar(CM_)
! ! grav1 = (shpe(q_)/r1**3.d0) * rad_andreas_tmp ! - (shpe(q_)*(one-shpe(gmEdd_))) / r1**3.d0
! ! grav1x = grav1 * r1x
! ! grav1y = grav1 * r1y
! ! grav1z = grav1 * r1z
! ! grav2 = - one                           / r2**3.d0
! ! grav2x = grav2 * r2x
! ! grav2y = grav2 * r2y
! ! grav2z = grav2 * r2z
! ! cent = (one+shpe(q_))/eqpar(a_)**3.d0 ! *eqpar(egg_)**3.d0
! ! centx = cent * rx
! ! centy = cent * ry
! ! centz = zero
! ! cor = - two * dsqrt(one+shpe(q_)) / eqpar(a_)**1.5d0 ! * eqpar(egg_)**1.5d0
! ! corx = cor * (-w(2))
! ! cory = cor * ( w(1))
! ! corz = zero
!
! grav1 = rad_andreas_tmp / r1 ! - (shpe(q_)*(one-shpe(gmEdd_))) / r1**3.d0
! grav1x = grav1 * r1x
! grav1y = grav1 * r1y
! grav1z = grav1 * r1z
! grav2 = - (shpe(eta_)**(-2.d0)/(shpe(fll_)*eqpar(egg_)*(one+shpe(q_)))) / r2**3.d0
! grav2x = grav2 * r2x
! grav2y = grav2 * r2y
! grav2z = grav2 * r2z
! cent = ((shpe(fll_)*eqpar(egg_))/shpe(eta_))**2.d0 ! (one+shpe(q_))/eqpar(a_)**3.d0 ! *eqpar(egg_)**3.d0
! centx = cent * rx
! centy = cent * ry
! centz = zero
! cor = - two * ((shpe(fll_)*eqpar(egg_))/shpe(eta_)) ! dsqrt(one+shpe(q_)) / eqpar(a_)**1.5d0 ! * eqpar(egg_)**1.5d0
! corx = cor * (-w(2))
! cory = cor * ( w(1))
! corz = zero
!
! ! glx = gl * r1x
! ! gly = gl * r1y
! ! glz = gl * r1z
!
! ! Full case : binary system
! if (with_compact_object) then
! dw(1) = grav1x + grav2x + centx + corx ! + glx
! dw(2) = grav1y + grav2y + centy + cory ! + gly
! dw(3) = grav1z + grav2z + centz + corz ! + glz
! else
! dw(1) = grav1x + centx + corx ! + glx
! dw(2) = grav1y + centy + cory ! + gly
! dw(3) = grav1z + centz + corz ! + glz
! endif
! ! Special case : isolated star
! ! dw(1) = grav1x + glx
! ! dw(2) = grav1y + gly
! ! dw(3) = grav1z + glz

end subroutine derivs
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
subroutine rk4(X,qdt,Xdot)
include 'def.f90' ! To access nw
double precision, intent(inout) :: X
double precision, intent(in) :: qdt, Xdot
double precision :: Xt, dXt, dXm
double precision:: qdt2, qdt6
integer :: i

qdt2=half*qdt
qdt6=qdt/6.d0

Xt=X+qdt2*Xdot
! do i=1,3
!    xt(i)=x(i)+qdt2*dx(i)
!    wt(i)=w(i)+qdt2*dw(i)
! enddo
call derivs(Xt,dXt)
Xt=X+qdt2*dXt
! do i=1,3
!    xt(i)=x(i)+qdt2*dxt(i)
!    wt(i)=w(i)+qdt2*dwt(i)
! enddo
call derivs(Xt,dXm)
Xt=X+qdt*dXm
dXm=dXt+dXm
! do i=1,3
!    xt(i)=x(i)+qdt*dxm(i)
!    wt(i)=w(i)+qdt*dwm(i)
!    dxm(i)=dxt(i)+dxm(i)
!    dwm(i)=dwt(i)+dwm(i)
! enddo
call derivs(Xt,dXt)
X=X+qdt6*(Xdot+dXt+2.d0*dXm)
! call derivs(xt,wt,dxt,dwt,with_compact_object)
! do i=1,3
!    x(i)=x(i)+qdt6*(dx(i)+dxt(i)+two*dxm(i))
!    w(i)=w(i)+qdt6*(dw(i)+dwt(i)+two*dwm(i))
! enddo
! print*, w(1)

end subroutine rk4
! -----------------------------------------------------------------------------------

end module integration
