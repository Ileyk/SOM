program tides
! use init
use core_module
use IO
use miscellaneous
use preliminaries
include 'def.f90'

call cpu_time(chrono_0)
! Compilation - - - - - - -
! gfortran miscellaneous.f90 mod_kracken.f90 draw.f90 IO.f90 cart_and_sph.f90 hull.f90 rad_acc_mod.f90 parameters.f90 init.f90 integration.f90 postprocess.f90 core_module.f90 tides.f90 -o tides
! gfortran miscellaneous.f90 mod_kracken.f90 IO.f90 preliminaries.f90 integration.f90 core_module.f90 tides.f90 -o tides

call give_filenames_v2

call read_arguments()
call read_par_files()

print*, 'Preliminary checks - - - - - - -'

! Compute quantities from 4 D.O.F.
! (eg initial tot. ang. mom., min. tot. ang. mom.
! to enable synchronization and corresponding omega)
call cmpt_params
! Check Darwin
call prm_chcks

print*, 'Preliminary checks - - - - - - - -'

print*, 'Core compute starts - - - - - - - - - -'

call cr_cmpt

print*, 'Terminated - - - - - - - - - - - -'

call chrono(chrono_0,chrono_0_mess)

end program tides
