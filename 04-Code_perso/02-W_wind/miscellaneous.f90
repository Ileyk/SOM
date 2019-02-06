module miscellaneous
contains

! -----------------------------------------------------------------------------------
subroutine crash(message)
character(LEN=*), intent(in) :: message
write(3,*), message
stop
end subroutine crash
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
subroutine followup(message)
include "def.f90"
character(LEN=*), intent(in) :: message
open(4,file=log_file,position="append")
write(4,*), message
close(4)
end subroutine followup
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
! Computes the time ellapsed from start
subroutine chrono(start,start_mess)
include 'def.f90'
double precision, intent(in) :: start
character(len=8), intent(in) :: start_mess ! to have the name of the chrono
double precision :: finish
character(len=200) :: duration

call cpu_time(finish)
write(duration,'(F7.1)') finish-start
duration='It took '//trim(duration)//' sec since '//trim(start_mess)
call followup(duration)

end subroutine chrono
! -----------------------------------------------------------------------------------


end module miscellaneous
