implicit none

double precision, parameter :: smalldble=1.d-10, bigdble=1.d99
integer         , parameter :: intbigdble=100000000
double precision, parameter :: one=1.d0, dpi=dacos(-1.d0), zero=0.d0, half=0.5d0, two=2.d0
character(len=*), parameter :: undefined = 'undefined'

! CGS
double precision, parameter :: kboltz = 1.3806488D-16
double precision, parameter :: stefan = 5.670373d-5
double precision, parameter :: mnucl = 1.6733D-24
double precision, parameter :: msun = 1.989D+33
double precision, parameter :: rsolar = 7.0D+10
double precision, parameter :: secinday= 8.64d4
double precision, parameter :: secinyear= 3.1557600d7
double precision, parameter :: ggrav  = 6.67384D-8
double precision, parameter :: clight = 2.99792458d10
double precision, parameter :: Lsun = 3.846d23
double precision, parameter :: AU = 1.49597871d13
double precision, parameter :: Jsun = 1.1d47 ! ang. mom. for 25d
double precision, parameter :: sigmaE = 6.652458734d-25 ! = sigma_thomson_scattering_on_free_electrons
double precision, parameter :: kpE = 3.98d-1 ! 3.1d-2! = sigma_thomson / mass_proton.
! For kpE, we used the formula (65) of KUDRITZKI+89 w/ I_HE=2 & N_He/N_H = 0.2

! 3 physical D.O.F. + IC
DOUBLE PRECISION :: R_, tau_, q_, om0_
! numerical parameters
double precision :: tmax, CFL, accrcy
! Quantities computed from 4 D.O.F.
double precision :: Jtot_ini, Jc, OmC

double precision :: chrono_0, chrono_1, chrono_2, chrono_3, chrono_4, chrono_5
COMMON chrono_0, chrono_1, chrono_2, chrono_3, chrono_4, chrono_5
character(len=8), parameter :: chrono_0_mess='chrono_0', chrono_1_mess='chrono_1', &
chrono_2_mess='chrono_2', chrono_3_mess='chrono_3', chrono_4_mess='chrono_4', chrono_5_mess='chrono_5'

COMMON R_, tau_, q_, om0_
COMMON tmax, CFL, accrcy
COMMON Jtot_ini, Jc, OmC
character(len=500) :: folder
COMMON folder
character(len=400) :: mess
COMMON mess
character(len=400) :: log_file, out_file, err_file, par_file, prm_file
COMMON log_file, out_file, err_file, par_file, prm_file
logical :: mrgr=.false.
COMMON mrgr
