#include <AMReX_winstd.H>

#include <NavierStokes.H>
#include <NS_error_F.H>

using std::string;

void
NavierStokes::error_setup()
{
    // The lines below define routines to be called to tag cells for error
    // estimation -- the arguments of each "add" call are:
    //   1. Name of variable (state variable or derived quantity) which will be
    //      passed into the Fortran subroutine.
    //   2. Number of ghost cells each array needs in each call to the Fortran
    //      subroutine.
    //   3. Type of Fortran subroutine -- this determines the argument list of
    //      the Fortran subroutine.  These types are pre-defined.
    //   4. Name of Fortran subroutine.

    //
    // **************  DEFINE ERROR ESTIMATION QUANTITIES  *************
    //
    err_list.add("mag_vort", 0, ErrorRec::Special,
                  BL_FORT_PROC_CALL(FORT_MVERROR,fort_mverror));
}
