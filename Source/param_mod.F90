
module param

  use amrex_fort_module, only : rt => amrex_real

  ! Phys bcs as defined in incflo 
   integer, parameter :: undef_cell =   0 ! undefined
   integer, parameter :: pinf_      =  10 ! pressure inflow cell
   integer, parameter :: pout_      =  11 ! pressure outflow cell
   integer, parameter :: minf_      =  20 ! mass flux inflow cell
   integer, parameter :: nsw_       = 100 ! wall with no-slip b.c.
   integer, parameter :: fsw_       = 101 ! wall with free-slip
   integer, parameter :: psw_       = 102 ! wall with partial-slip b.c.
   integer, parameter :: cycl_      =  50 ! cyclic b.c.
   integer, parameter :: cycp_      =  51 ! cyclic b.c. with pressure drop

! Common constants
   real(rt), parameter :: zero = 0.0d0
   real(rt), parameter :: half = 0.5d0
   real(rt), parameter :: one  = 1.0d0
   real(rt), parameter :: two  = 2.0d0
   real(rt), parameter :: my_huge  = 1.0d20
   
 end module param
