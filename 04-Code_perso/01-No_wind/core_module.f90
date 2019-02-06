module core_module
contains

! ! -----------------------------------------------------------------------------------
subroutine cr_cmpt
use miscellaneous
use integration
include 'def.f90'

double precision :: X
double precision :: Omrb, om, Jtot ! orbital and stellar spin pulsations
double precision :: Omrb_old, om_old
double precision :: qt, qdt ! time and time step
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call followup("First step of RK4 - - - - - - - - - - - - - - - - - - - - - - - -")

call cpu_time(chrono_1)

! Initialization of orbital and spin pulsations
Omrb_old = 1.d0
om_old   = om0_
X        = 1.d0

qt=0.d0  ! initialize the time

open(9,file=out_file)

do while (qt<tmax)

  ! Compute Omorb and om from X
  Omrb = X**(-3.d0)
  om   = om0_ - (5.d0/(R_**2.d0*(1.d0+q_))) * (X-1.d0)

  Jtot = X + ((R_**2.d0*(1.d0+q_))/5.d0)*om

  write(9,*) qt/(2.*dpi), qdt, X, Omrb, om, Jtot

  call advance(X,qt)

  call check_mrgr(X) ! simply when a < R for now, but RLOF rather later on

  call check_Jtot_cons(Jtot)

enddo

close(9)

call chrono(chrono_1,chrono_1_mess)


!          call check_exit(x,qt,time_to_stop,flag)
!
!          ! call get_hull(x,11,i0max,j0max,r_hull,t_hull,p_hull,hull_mask)
!
!           ! write(16,*) x(1), x(2), x(3)
!
!           if (time_to_stop .and. flag) then
!             if(x(1)>0.) then
!               ph_tmp = datan(x(2)/x(1))
!             else
!               ph_tmp = dpi + datan(x(2)/x(1))
!             endif
!             iphi = (ph_tmp-(-dpi/two)) / (2.*dpi/16.) + 1
!             ! print*, x(1),x(2),ph_tmp,iphi
!             if (x(3)>max_(1,iphi)) then
!               max_(1,iphi)=x(3)
!               max_(2,iphi)=dsqrt(x(1)**two+(x(2)+eqpar(CM_))**two)
!               ! print*, pos(ind,1), pos(ind,2)
!               ! print*, max_(1), max_(2)
!               ! print*, i0, j0
!             endif
!           endif
!           ! if ( i0==i0max ) then
!
!          ! Main addition (apart from different normalization) for ULX paper October 2018 :
!          ! the main measure, fraction of Mdot accreted, is made here
!          ! if (time_to_stop .and. flag) then
!          !   print*, ind," accepted"
!          !   surf_star_accreted=surf_star_accreted+delta_ph*(dcos(pos(ind,1)-delta_th/two)-dcos(pos(ind,1)+delta_th/two))
!          ! endif
!
!       enddo
!
!       ! if (max_(1)>true_max(1) .and. max_(1)<max_(2)) then
!       !   true_max(1)=max_(1)
!       !   true_max(2)=max_(2)
!       !   ! print*, true_max(1), true_max(2)
!       ! endif
!
!       ! If not accepted, do not save the initial pos
!       ! (which will be overwrite at next (i0,j0)... except if it is the last! so we set to zero by hand)
!       ! ind identifies the accepted ones
!       if (flag .EQV. .false.) then
!          pos(ind,1)=zero
!          pos(ind,2)=zero
!          ind=ind-1
!       endif
!
!    enddo
! enddo
!
! do iphi=1,16
!   print*, iphi, max_(1,iphi), max_(2,iphi), max_(1,iphi)/max_(2,iphi)
! enddo
! do iphi=1,16/4
!   max_(1,minloc(max_(1,:)))=1.d99
! enddo
! print*, 'median (in phi) aspect ratio : ', max_(1,minloc(max_(1,:)))/max_(2,minloc(max_(1,:))), &
!   one/dsin(datan(max_(1,minloc(max_(1,:)))/max_(2,minloc(max_(1,:))))), &
!   'for eta=', shpe(eta_), 'and beta=', shpe(bet_), '@ r=', r_end,'stellar radii'
!
! ! to plot the contour of the simulation space
! ! do i0=1,300
! !   ! write(16,*) eqpar(a_)+eqpar(rout_)*dcos((float(i0)/300.d0)*two*dpi), eqpar(rout_)*dsin((float(i0)/300.d0)*two*dpi), one
! ! aa=eqpar(NStoL1_)*dsin((float(i0)/300.d0)*two*dpi)
! ! bb=-one/(shpe(fll_)*eqpar(egg_))+eqpar(NStoL1_)*dcos((float(i0)/300.d0)*two*dpi)
! ! write(16,*)aa,bb,one
! ! enddo
!
! close(16)
!
! write(4,*) ind, 'in ', i0max*j0max, 'streamlines accepted'
! inquire(file=info_file, exist=exist)
! if (exist) then
!   open(29, file=info_file, status="old", position="append", action="write")
! else
!   open(29, file=info_file, status="new", action="write")
!   write(29,*) ' q = M_star / M_compact | beta | eta = v_inf / aOmega | fll || bof | OK | yes | yes+ | bof '! header
! end if
! write(29,'(4F12.7)',advance='no') shpe(q_), shpe(bet_), shpe(eta_), shpe(fll_)
! write(29,'(F12.7)',advance='no') ( dble(2*ind)/dble((2*i0max)*j0max) ) * ( om_angle / (4.d0*dpi) )
! write(29,'(F12.7)',advance='no') (2.*surf_star_accreted/(4.d0*dpi))
! print*, 100.*(2.*surf_star_accreted/(4.d0*dpi)), '%'
! print*, 100.d0 * ( dble(2*ind)/dble((2*i0max)*j0max) ) * ( om_angle / (4.d0*dpi) ) , &
!            '% : 1st estimation of the fraction of streamlines entering rout'
! ! check whether the # of accepted streamline will be ~ to the # of entry cells
! ! if (ind*4<Nth_rmp*Nph_rmp*5) call crash('Too few streamlines accepted')
! ! if (dble(ind)/dble(i0max*j0max)<0.3d0) call crash('Too few streamlines accepted')
! ! 1st estimation of the fraction of the stellar wind outflow mass rate entering within eqpar(Rout_) around CO
! ! Factor two since we only consider the up part so the # of accepted lines within 4pi would be twice as large
! ! write(4,*) 100.d0 * ( dble(2*ind) / (((4.d0*dpi)/(two*om_angle))*dble(i0max*j0max)) ) , &
! !            '% : 1st estimation of the fraction of streamlines entering rout'
! write(4,*) 100.d0 * ( dble(2*ind)/dble((2*i0max)*j0max) ) * ( om_angle / (4.d0*dpi) ) , &
!            '% : 1st estimation of the fraction of streamlines entering rout'
! stop
! ! Sanity check to make sure that working on a restricted area
! ! with th_min /=0 and phi not all around the star does not prevent
! ! good streamlines to be taken into account (although other non adjacent good
! ! streamlines could still exist while not being detected by this sanity check,
! ! but @ least, it is a necesary condition that this sanity check is passed)
! ! if ( minval(pos(:,1),MASK=(pos(:,1)/=zero)) == th_min + delta_th * (dble(1    )-half) .or. &
! !      minval(pos(:,2),MASK=(pos(:,2)/=zero)) == ph_min + delta_ph * (dble(1    )-half) .or. &
! !      maxval(pos(:,2),MASK=(pos(:,2)/=zero)) == ph_min + delta_ph * (dble(j0max)-half)) then
! ! print*, minval(pos(:,1),MASK=(pos(:,1)/=zero)) == th_min + delta_th * (dble(1    )-half)
! ! print*, minval(pos(:,2),MASK=(pos(:,2)/=zero)) == ph_min + delta_ph * (dble(1    )-half)
! ! print*, maxval(pos(:,2),MASK=(pos(:,2)/=zero)) == ph_min + delta_ph * (dble(j0max)-half)
! ! call crash("Error#20 : the initial patch is too restrictive. Enlarge it w/ th_min, ph_min & ph_max.")
! ! endif
! ! Put in pos_old the initial position pos of the pts which worked
! allocate(pos_old(ind,2))
! j0=1
! do i0=1,i0max*j0max
!    if (pos(i0,1)/=zero .and. pos(i0,2)/=zero) then
!       pos_old(j0,:)=pos(i0,:)
!       j0=j0+1
!    endif
! enddo
!
! deallocate(pos)
! id=0
! ! Refinement
! delta_th=half*delta_th
! delta_ph=half*delta_ph
!
! ! we do one refinement loop only
! call chrono(chrono_1,chrono_1_mess)
!
! open(16,file=draw_good_in_plane)
!
! call cpu_time(chrono_2)
!
! ! print*, 'voui',dlog(20.*dpi*delta_ph*2./dsqrt((dble(2*ind)/dble((2*i0max)*j0max))*om_angle))/dlog(two)
! ! print*, delta_ph*2., dsqrt((dble(2*ind)/dble((2*i0max)*j0max))*om_angle/dpi)
! ! kmax=1+int(dlog(20.*dpi*delta_ph*2./dsqrt((dble(2*ind)/dble((2*i0max)*j0max))*om_angle))/dlog(two))
! kmax=1
! print*, 'kmax = ', kmax
! ! - - - -
! do k=1,kmax
!   call cpu_time(chrono_3)
!   call followup("Following steps of RK4 - - - - - - - - - - - - - - - - - - - - - - -")
!   if (k==kmax) then
!     allocate(arriv_pos(ind*4,3),arriv_vit(ind*4,3))!,pos_ini_tmp(ind*4,3))
!     arriv_pos=zero
!     arriv_vit=zero
!   endif
!   ! allocate(arriv_pos(ind*4,3),arriv_vit(ind*4,3),length(ind*4))
!   ! arriv_pos=zero
!   ! arriv_vit=zero
!   ! length=zero
!   allocate(pos(ind*4,2)) ! Maximum # of future accepted after refinement
!   pos=zero
!
!   surf_star_accreted=zero
!
!   ! For each successful patch...
!   do i0=1,ind
!     ! ... build the 4 neighbors...
!     do j0=1,4
!        ! % of advancement
!        if (int(dble(4*ind)/100.d0)>zero) then
!        if (mod(j0+(i0-1)*4,int(dble(4*ind)/10.d0))==0) then
!           write(mess,'(I4.4)'), (j0+(i0-1)*4)/int(dble(4*ind)/100.d0)
!           mess=trim(mess)//'%'
!           call followup(trim(mess))
!        endif
!        endif
!
!        flag=.true.
!        ! Initialisation
!        x(1) = dist_andreas(1)*1.00001d0 ! to avoid issue at interpolation ! eqpar(Rstr_) !*eqpar(ratioRadiiIni_) !R! ! r
!        w(1) = zero
!        w(2) = zero
!        w(3) = zero
!
!        select case(j0)
!        case(1)
!           x(2)=pos_old(i0,1)-delta_th/two
!           x(3)=pos_old(i0,2)-delta_ph/two
!        case(2)
!           x(2)=pos_old(i0,1)+delta_th/two
!           x(3)=pos_old(i0,2)-delta_ph/two
!        case(3)
!           x(2)=pos_old(i0,1)-delta_th/two
!           x(3)=pos_old(i0,2)+delta_ph/two
!        case(4)
!           x(2)=pos_old(i0,1)+delta_th/two
!           x(3)=pos_old(i0,2)+delta_ph/two
!        end select
!
!        ! For next iteration k, save the launching position of the patch (canceled if not accepted)
!        id=id+1
!        pos(id,1)=x(2)
!        pos(id,2)=x(3)
!
!        ! Strict repeat from initialization above ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!        call sph_to_cart(x,w)
!
!        qt=0.d0  ! initialize the time
!
!        ! Integration
!        time_to_stop=.false.
!
!        do while (time_to_stop .EQV. .false.)
!
!           xtmp=x
!           call advance(x,w,qt,.true.)
!           ! length(id)=length(id)+dsqrt((x(1)-xtmp(1))**2.+(x(2)-xtmp(2))**2.+(x(3)-xtmp(3))**2.)
!           call check_exit(x,qt,time_to_stop,flag)
!
!           ! Save only the trajectories closest from the equatorial plane
!           ! (and above the equatorial plane)
!           ! for use in draw_me_a_sheep
!           if ( pos_old(i0,1)-dpi/two<zero .and. &
!               dabs(pos_old(i0,1)-dpi/two)==minval(dabs(pos_old(:,1)-dpi/two)) .and. &
!               (j0==2 .or. j0==4) ) then
!             ! write(16,*) x(1), x(2), x(3)
!           endif
!
!          if (k==kmax .and. time_to_stop .and. flag) then
!            surf_star_accreted=surf_star_accreted+delta_ph*(dcos(pos(id,1)-delta_th/two)-dcos(pos(id,1)+delta_th/two))
!          endif
!
!        enddo
!
!        ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!        ! If not accepted, do not save the initial pos
!        ! (which will be overwrite at next (i0,j0)... except if it is the last! so we set to zero by hand)
!        ! ind identifies the accepted ones
!        if (flag .EQV. .false.) then
!           pos(id,1)=zero
!           pos(id,2)=zero
!           ! length(id)=zero
!           id=id-1
!        endif
!        ! Save the final good points and velocities for the last iteration
!        if (k==kmax .and. (flag .EQV. .true.)) then
!           arriv_pos(id,:)=x(1:3)
!           arriv_vit(id,:)=w(1:3)
!           ! x(2)=x(2)+eqpar(a_)
!           ! call cart_to_sph_1pt(x,w)
!           ! write(16,*) x(1), x(2), x(3)
!        endif
!     enddo ! j0 the neighbors
!   enddo ! i0 the accepted
!   ! write(4,*) 'curvilinear length ranges between', minval(length,MASK=(length/=0.)), 'and', maxval(length)
!
!   ! write(4,*) 100.d0*(dble(id)/dble(4.d0*ind)),"% accepted after refinement"
!   write(4,*) id,"accepted for",Nth_rmp*Nph_rmp,"entry points"
!
!   deallocate(pos_old)
!   allocate(pos_old(id,2))
!   if (count(pos(:,1)>smalldble)/=id) call crash("Error#09 : id /= 4*ind. WTF?!")
!   do i0=1,id
!        pos_old(i0,1)=pos(i0,1)
!        pos_old(i0,2)=pos(i0,2)
!   enddo
!   deallocate(pos)
!   ind=id
!   id=0
!   delta_th=half*delta_th
!   delta_ph=half*delta_ph
!   call chrono(chrono_3,chrono_3_mess)
!
! enddo
!
! close(16)
!
! ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
! ! Here, we use the data from file#16 (alias draw_tmp which is erased @ the end of the processe)
! ! to draw all the selected streamlines closest to the equatorial plane
! ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
! write(29,'(F12.7)',advance='no') (2.*surf_star_accreted/(4.d0*dpi))
! print*, 100.*(2.*surf_star_accreted/(4.d0*dpi)), '%'
!
! ! Draw the config + some streamlines
! call draw_me_a_sheep_v2
! ! Remove the draw_tmp file
! call system('rm '//trim(folder)//'draw_all_in_plane')
! ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
!
! close(1)

call chrono(chrono_2,chrono_2_mess)

! ! write(4,*) 'Final # of pts arrived, to be processed : ', ind
! ! Final # of pts arrived
! if (ind==0) call crash('Error#10 : no good points detected. Rise the resolution.')
! allocate(pos(ind,3),pos_ini(ind,3))
! ! pos_ini in spherical
! pos_ini(1:ind,1)=one ! eqpar(Rstr_)
! pos_ini(1:ind,2:3)=pos_old(1:ind,1:2)
! allocate(pos_old_2(ind,nw))
!
! if (count(arriv_pos(:,1)/=zero)/=size(pos_old_2,DIM=1) .or. &
!     count(arriv_pos(:,1)/=zero)/=size(pos    ,DIM=1)) then
!    write(3,*), 'Error#11 : inconsistency. Should be = ind', count(arriv_pos(:,1)/=zero), size(pos,DIM=1), size(pos_old_2,DIM=1)
!    stop
! endif
!
! ! print*, eqpar(a_)
! ! arriv_pos(:,2)=arriv_pos(:,2)+eqpar(a_)
! ! print*, minval(arriv_pos(:,1),MASK=(arriv_pos(:,3)/=zero)), maxval(arriv_pos(:,1),MASK=(arriv_pos(:,3)/=zero))
! ! print*, minval(arriv_pos(:,2),MASK=(arriv_pos(:,3)/=zero)), maxval(arriv_pos(:,2),MASK=(arriv_pos(:,3)/=zero))
! ! print*, minval(arriv_pos(:,3),MASK=(arriv_pos(:,3)/=zero)), maxval(arriv_pos(:,3),MASK=(arriv_pos(:,3)/=zero))
! ! call cart_to_sph(arriv_pos,arriv_vit,ind)
! ! print*, minval(arriv_pos(:,1),MASK=(arriv_pos(:,3)/=zero)), maxval(arriv_pos(:,1),MASK=(arriv_pos(:,3)/=zero))
! ! if (dabs(maxval(arriv_pos(:,1))-minval(arriv_pos(:,1),MASK=(arriv_pos(:,1)/=zero)))/&
! !                            dabs(minval(arriv_pos(:,1),MASK=(arriv_pos(:,1)/=zero)))>1.d-5) then
! !    print*, dabs(maxval(arriv_pos(:,1))-minval(arriv_pos(:,1),MASK=(arriv_pos(:,1)/=zero)))/&
! !                               dabs(minval(arriv_pos(:,1),MASK=(arriv_pos(:,1)/=zero)))
! !    print*, 'Error#4 : why such a discrepancy in arrival radius? Should be ~R0'
! !    stop
! ! endif
!
! ! pos and pos_old_2 play the role of final position and velocity
! ! pos and pos_old_2 in cartesian while
! pos(1:ind,:)=arriv_pos(1:ind,:)
! pos_old_2(1:ind,1:3)=arriv_vit(1:ind,1:3)
! ! Also notice that, since they are indexed w/ ind, pos_ini, pos and pos_old_2
! ! all three corresponds to the same trajectory
! deallocate(arriv_pos,arriv_vit)!,pos_ini_tmp)
! pos_old_2(:,4:nw)=zero
! call cpu_time(chrono_4)
! first_call=.false.
! call process(pos,pos_ini,pos_old_2,ind,first_call) ! length
! ! To modify by adding an argument to prevent the writing of pts and the call to remap
! call chrono(chrono_4,chrono_4_mess)

! deallocate(pos,pos_ini,pos_old_2)
!- - - -

! close(2)

end subroutine cr_cmpt
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
subroutine advance(X,qt)
use integration
include 'def.f90'
double precision, intent(inout) :: X
double precision, intent(inout) :: qt
double precision :: qdt, Xdot

! call getdt(x,w,qdt,with_compact_object)
call getdt(X,qdt)

call derivs(X,Xdot)

call rk4(X,qdt,Xdot)

qt=qt+qdt

end subroutine advance
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
subroutine check_mrgr(X,mrgr)
use miscellaneous
include 'def.f90'
double precision, intent(in) :: X
logical, intent(out) :: mrgr

if (X**2.d0<R_) then
  mrgr=.true.
  ! call crash('Error#01 : the system merged')
endif

end subroutine check_mrgr
! ! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine check_Jtot_cons(Jtot)
use miscellaneous
include 'def.f90'
double precision, intent(in) :: Jtot
double precision :: err

err=dabs(Jtot-Jtot_ini)/Jtot_ini
if (err>accrcy) then
  print*, 'Total angular momentum not conserved', Jtot_ini, Jtot, err,'>',accrcy
  call crash('Error#02 : inaccurate integration, total angular momentum not conserved')
endif

end subroutine check_Jtot_cons
! -----------------------------------------------------------------------------------

end module core_module
