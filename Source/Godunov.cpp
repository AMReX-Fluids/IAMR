
//
// Godunov is the object which calculates advective terms for iamr.
//

#include <AMReX_LO_BCTYPES.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <Godunov.H>
#include <iamr_godunov.H>
#include <GODUNOV_F.H>

#include <algorithm>

#define GEOM_GROW 1
#define XVEL 0
#define YVEL 1
#define ZVEL 2

using namespace amrex;

namespace
{
    bool initialized = false;

    const int use_unlimited_slopes_DEF = 0;
}
//
// Set default values for these in Initialize()!!!
//
static int hyp_grow;

int  Godunov::verbose;
int  Godunov::ppm_type;
int  Godunov::slope_order;
int  Godunov::corner_couple;
int  Godunov::use_forces_in_trans;

void
Godunov::Initialize ()
{
    if (initialized) return;
    //
    // Set defaults here!!!
    //
    hyp_grow = 3;

    Godunov::verbose             = 0;
    Godunov::ppm_type            = 0;
    Godunov::slope_order         = 4;
    Godunov::corner_couple       = 1;
    Godunov::use_forces_in_trans = 0;

    int use_unlimited_slopes = use_unlimited_slopes_DEF;

    ParmParse pp("godunov");

    pp.query("v",                    verbose);
    pp.query("ppm_type",             ppm_type);
    pp.query("slope_order",          slope_order);
    pp.query("corner_couple",        corner_couple);
    pp.query("use_forces_in_trans",  use_forces_in_trans);
    pp.query("use_unlimited_slopes", use_unlimited_slopes);

    if (ppm_type == 2)
    {
	hyp_grow = 4;
    }

#ifdef AMREX_USE_EB

//     //nghost from incflo.

    hyp_grow = 4;
#endif

#if (BL_SPACEDIM==2)
    BL_ASSERT(slope_order==1 || slope_order==2 || slope_order==4);
#else
    BL_ASSERT(slope_order==1 || slope_order==4);
#endif

    set_params(slope_order, use_unlimited_slopes);

    amrex::ExecOnFinalize(Godunov::Finalize);

    initialized = true;
}

void
Godunov::Finalize ()
{
    initialized = false;
}

//
// Construct the Godunov Object.
//
Godunov::Godunov (int max_size)
{
  Initialize();
}

Godunov::~Godunov ()
{
    ;
}



int
Godunov::hypgrow ()
{
    return hyp_grow;
}
