
#include <NavierStokes.H>
#include <NS_error_F.H>

using std::string;

using namespace amrex;

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
    if (do_density_ref)   {
        err_list.add("density",  1, ErrorRec::Special,FORT_DENERROR);
	amrex::Print() << "Refining on DENSITY" << std::endl;
    }
    if (do_tracer_ref)    {
        err_list.add("tracer",   1, ErrorRec::Special,FORT_ADVERROR);
	amrex::Print() << "Refining on TRACER" << std::endl;
    }
    if (do_tracer2_ref)    {
	err_list.add("tracer2",   1, ErrorRec::Special,FORT_ADV2ERROR);
	amrex::Print() << "Refining on TRACER2" << std::endl;
    }
    if (do_vorticity_ref) {
        err_list.add("mag_vort", 0, ErrorRec::Special,FORT_MVERROR);
	amrex::Print() << "Refining on MAG_VORT" << std::endl;
    }
    if (do_temp_ref) {
        err_list.add("temp", 1, ErrorRec::Special, FORT_TEMPERROR);
	amrex::Print() << "Refining on TEMP and/or GRAD T" << std::endl;
    }
}
