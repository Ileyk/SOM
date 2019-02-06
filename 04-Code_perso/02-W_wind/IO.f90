module IO
  implicit none
  public

contains

! -----------------------------------------------------------------------------------
!> Read the command line arguments passed to amrvac
subroutine read_arguments()
  use mod_kracken
  include 'def.f90'

  integer                          :: len, ier, n
  integer                          :: ibegin
  integer                          :: iterm

  ! Specify the options and their default values
  call kracken('cmd','-i amrvac.par -if ' // undefined // &
       ' -slice 0 -collapse 0 --help .false. -convert .false.')

  ! Get the par file(s)
  call retrev('cmd_i', par_file, len, ier)

  ! Split the input files, in case multiple were given
  call get_fields_string(par_file, " ,'"""//char(9))

end subroutine read_arguments
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
!> Read in the user-supplied parameter-file
subroutine read_par_files()
  use miscellaneous
  include 'def.f90'

  integer :: unitpar=9
  logical            :: file_exists
  double precision   :: R, tau, q, om0

  namelist /physical/ R, tau, q, om0

  namelist /integration/ tmax, CFL, accrcy

  open(3,file=err_file)

  R     = bigdble
  tau   = bigdble
  q     = bigdble
  om0   = bigdble

  tmax   = bigdble
  CFL    = bigdble
  accrcy = bigdble

  print *, "Reading " // trim(par_file)

  ! Check whether the file exists
  inquire(file=trim(par_file), exist=file_exists)

  if (.not. file_exists) &
    call crash("The parameter file " // trim(par_file) // " does not exist")

  open(unitpar, file=trim(par_file), status='old')
  ! Try to read in the namelists. They can be absent or in a different
  ! order, since we rewind before each read.
  rewind(unitpar)
  read(unitpar, physical, end=101)

101    rewind(unitpar)
       read(unitpar, integration, end=102)

102    close(unitpar)

! physical - - -
if (R>bigdble/two .or. tau>bigdble/two .or. q>bigdble/two .or. om0>bigdble/two) &
  call crash("Missing physical parameters")
R_   = R
tau_ = tau
q_   = q
om0_ = om0

! integration - - -
if (tmax>bigdble/two) call crash("Missing max integration time")
if (CFL>bigdble/two) call crash("Missing integration parameters")
if (accrcy>bigdble/two) call crash("Missing accuracy needed on L conservation")

end subroutine read_par_files
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine give_filenames_v2
include 'def.f90'
character(LEN=10) :: f_name, q_name, M2_name !, a_name, b_name!, af_name

folder='./out/'

! folders
log_file=trim(folder)//'log'
out_file=trim(folder)//'out'
err_file=trim(folder)//'err'
prm_file=trim(folder)//'prm'

! Since lines are appended in log_file each time (see in miscellaneous.f90),
! we need to erase the previous one first
call system("rm -f "//log_file)
! et tant qu a faire...
call system("rm -f "//err_file)

end subroutine give_filenames_v2
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
!> Routine to find entries in a string
subroutine get_fields_string(line, delims)
  !> The line from which we want to read
  character(len=*), intent(inout) :: line
  !> A string with delimiters. For example delims = " ,'"""//char(9)
  character(len=*), intent(in)    :: delims

  integer :: ixs_start
  integer :: ixs_end
  integer :: ix, ix_prev

  ix_prev = 0

  ! Find the starting point of the next entry (a non-delimiter value)
  ix = verify(line(ix_prev+1:), delims)

  ixs_start = ix_prev + ix ! This is the absolute position in 'line'

  ! Get the end point of the current entry (next delimiter index minus one)
  ix = scan(line(ixs_start+1:), delims) - 1

  if (ix == -1) then              ! If there is no last delimiter,
    ixs_end = len(line) ! the end of the line is the endpoint
  else
    ixs_end = ixs_start + ix
  end if

  line = line(ixs_start:ixs_end)

end subroutine get_fields_string
! -----------------------------------------------------------------------------------

end module IO
