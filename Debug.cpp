//==========================================================
// This file contains functions for debugging iamr
// in octave
//==========================================================

#include <stdio.h>
#include <strstream.h>
#include <Misc.H>
#include <CoordSys.H>
#include <Geometry.H>
#include <BoxDomain.H>
#include <RunStats.H>
#include <ParmParse.H>
#include <ErrorList.H>
#include <Godunov.H>
#include <NavierStokes.H>
#include <MultiGrid.H>
#include <ArrayLim.H>

#include <netcdfIO.H>

//------------- a hack for sending all of the IO to one place
//------------- this eliminates having to type in paths in the debugger

extern aString PREFIXSTR;


//==========================================================
// debug functions for the NavierStokes object
//==========================================================

extern "C" void dump_aofs( NavierStokes *tmp  )
{
    mfab2cdf( PREFIXSTR.c_str(), tmp->aofs );
}

extern "C" void dump_umac( NavierStokes *tmp )
{
#if (BL_SPACEDIM == 2 )
    medge2cdf( PREFIXSTR.c_str(), &tmp->u_mac[0], &tmp->u_mac[1] );
#else
    medge2cdf( PREFIXSTR.c_str(), &tmp->u_mac[0], &tmp->u_mac[1], &tmp->u_mac[2] );
#endif
}


#if (BL_SPACEDIM == 2 )
extern "C" void dump_flux( FARRAYBOX *uflx, FARRAYBOX *vflx )
#else
extern "C" void dump_flux( FARRAYBOX *uflx, FARRAYBOX *vflx, FARRAYBOX *wflx )
#endif
{
#if (BL_SPACEDIM == 2 )
    edge2cdf( PREFIXSTR.c_str(), uflx, vflx );
#else
    edge2cdf( PREFIXSTR.c_str(), uflx, vflx, wflx );
#endif
}


// this function dumps the multifabs for testing out the nodal solver
extern "C" void dump_nodal( MultiFab *Unew, MultiFab *Pnew )
{
    mfab2cdf( "test/Unew", Unew );
    mfab2cdf( "test/Pnew", Pnew );
}    


extern "C" void dump_state( NavierStokes *tmp, int new_test )
{
    MultiFab &Snew = tmp->get_new_data(State_Type);
    MultiFab &Sold = tmp->get_old_data(State_Type);
    if (new_test)
        mfab2cdf( PREFIXSTR.c_str(), &Snew );
    else
        mfab2cdf( PREFIXSTR.c_str(), &Sold );
}


extern "C" void dump_press( NavierStokes *tmp, int new_test )
{
    MultiFab &Pnew = tmp->get_new_data(Press_Type);
    MultiFab &Pold = tmp->get_old_data(Press_Type);
    if (new_test)
        mfab2cdf( PREFIXSTR.c_str(), &Pnew );
    else
        mfab2cdf( PREFIXSTR.c_str(), &Pold );
}

extern "C" void dump_fab( FARRAYBOX *fab )
{
    fab2cdf( PREFIXSTR.c_str(), fab );
}

extern "C" void dump_mfab( MultiFab *mfab )
{
    mfab2cdf( PREFIXSTR.c_str(), mfab );
}



// This function dumps the state and pressure variablies
extern "C" void dump_all( NavierStokes *tmp )
{
    MultiFab &Snew = tmp->get_new_data(State_Type);
    MultiFab &Sold = tmp->get_old_data(State_Type);
    MultiFab &Pnew = tmp->get_new_data(Press_Type);
    MultiFab &Pold = tmp->get_old_data(Press_Type);

    mfab2cdf( "Snew", &Snew );
    mfab2cdf( "Sold", &Sold );
    mfab2cdf( "Pnew", &Pnew );
    mfab2cdf( "Pold", &Pold );
}

// This function reads the state and pressure variables
extern "C" void read_all( NavierStokes *tmp )
{
    MultiFab &Snew = tmp->get_new_data(State_Type);
    MultiFab &Sold = tmp->get_old_data(State_Type);
    MultiFab &Pnew = tmp->get_new_data(Press_Type);
    MultiFab &Pold = tmp->get_old_data(Press_Type);

    cdf2mfab( "Snew", &Snew );
    cdf2mfab( "Sold", &Sold );
    cdf2mfab( "Pnew", &Pnew );
    cdf2mfab( "Pold", &Pold );
}


//==========================================================
// debug functions for the Godunov object
//==========================================================

extern "C" void dump_uad( Godunov *tmp )
{
#if (BL_SPACEDIM == 2)
    edge2cdf( PREFIXSTR.c_str(), &tmp->uad, &tmp->vad );
#else
    edge2cdf( PREFIXSTR.c_str(), &tmp->uad, &tmp->vad, &tmp->wad );
#endif
}


extern "C" void dump_work( Godunov *tmp )
{
    fab2cdf( PREFIXSTR.c_str(), &tmp->work );
}


extern "C" void dump_1d( Godunov *tmp )
{
    int nitems;
    int sz = tmp->scr_size;
    FILE *dfile = fopen( "test/test_scr.dat", "wb" );

    nitems = sz;
    fwrite( &nitems, sizeof(int), 1, dfile );
    nitems = 3*(BL_SPACEDIM)*sz;
    fwrite( &nitems, sizeof(int), 1, dfile );

    nitems = fwrite( tmp->stxlo,  sizeof(REAL), sz, dfile );  assert(nitems == sz );
    nitems = fwrite( tmp->stxhi,  sizeof(REAL), sz, dfile );  assert(nitems == sz );
    nitems = fwrite( tmp->slxscr, sizeof(REAL), sz, dfile );  assert(nitems == sz );

    nitems = fwrite( tmp->stylo,  sizeof(REAL), sz, dfile );  assert(nitems == sz );
    nitems = fwrite( tmp->styhi,  sizeof(REAL), sz, dfile );  assert(nitems == sz );
    nitems = fwrite( tmp->slyscr, sizeof(REAL), sz, dfile );  assert(nitems == sz );
#if (BL_SPACEDIM == 3)
    nitems = fwrite( tmp->stzlo,  sizeof(REAL), sz, dfile );  assert(nitems == sz );
    nitems = fwrite( tmp->stzhi,  sizeof(REAL), sz, dfile );  assert(nitems == sz );
    nitems = fwrite( tmp->slzscr, sizeof(REAL), sz, dfile );  assert(nitems == sz );
#endif
    fclose(dfile);
}


