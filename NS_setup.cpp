//BL_COPYRIGHT_NOTICE

//
// $Id: NS_setup.cpp,v 1.36 2000-04-12 19:08:43 sstanley Exp $
//

#include <NavierStokes.H>
#include <RegType.H>
#include <ParmParse.H>
#include <ErrorList.H>
#include <PROB_F.H>
#include <DERIVE_F.H>
#include <Misc.H>
#include <FArrayBox.H>
#include <CoordSys.H>

static Box the_same_box (const Box& b)    { return b;           }
static Box grow_box_by_one (const Box& b) { return ::grow(b,1); }

//
// Components are  Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
//
static int norm_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD, EXT_DIR, EXT_DIR
};

static int tang_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, HOEXTRAP, EXT_DIR
};

static int scalar_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, FOEXTRAP
};

static int press_bc[] =
{
    INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, FOEXTRAP
};

static int temp_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, EXT_DIR, EXT_DIR
};

static int divu_bc[] =
{
    INT_DIR, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static int dsdt_bc[] =
{
    INT_DIR, EXT_DIR, EXT_DIR, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};


static
void
set_x_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,norm_vel_bc[lo_bc[0]]);
    bc.setHi(0,norm_vel_bc[hi_bc[0]]);
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_y_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
    bc.setLo(1,norm_vel_bc[lo_bc[1]]);
    bc.setHi(1,norm_vel_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

#if (BL_SPACEDIM == 3)
static
void
set_z_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
    bc.setLo(2,norm_vel_bc[lo_bc[2]]);
    bc.setHi(2,norm_vel_bc[hi_bc[2]]);
}
#endif

static
void
set_scalar_bc (BCRec&       bc,
               const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i,scalar_bc[lo_bc[i]]);
        bc.setHi(i,scalar_bc[hi_bc[i]]);
    }
}

static
void
set_temp_bc (BCRec&       bc,
             const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i,temp_bc[lo_bc[i]]);
        bc.setHi(i,temp_bc[hi_bc[i]]);
    }
}

static
void
set_pressure_bc (BCRec&       bc,
                 const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i,press_bc[lo_bc[i]]);
        bc.setHi(i,press_bc[hi_bc[i]]);
    }
}

static 
void 
set_divu_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
        bc.setLo(i,divu_bc[lo_bc[i]]);
        bc.setHi(i,divu_bc[hi_bc[i]]);
    }
}

static 
void 
set_dsdt_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
        bc.setLo(i,dsdt_bc[lo_bc[i]]);
        bc.setHi(i,dsdt_bc[hi_bc[i]]);
    }
}


void
NavierStokes::variableSetUp ()
{
    BL_ASSERT(desc_lst.length() == 0);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        phys_bc.setLo(dir,SlipWall);
        phys_bc.setHi(dir,SlipWall);
    }

    read_params();

    BCRec bc;
    //
    // Set number of state variables.
    //
    NUM_STATE = Density + 1;
    int Trac = NUM_STATE++;
    if (do_temp) NUM_STATE++;
    NUM_SCALARS = NUM_STATE - Density;

    //
    // **************  DEFINE VELOCITY VARIABLES  ********************
    //
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,1,NUM_STATE,
                           &cell_cons_interp);
    set_x_vel_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Xvel,"x_velocity",bc,FORT_XVELFILL);
    set_y_vel_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Yvel,"y_velocity",bc,FORT_YVELFILL);
#if (BL_SPACEDIM == 3)
    set_z_vel_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Zvel,"z_velocity",bc,FORT_ZVELFILL);
#endif
    //
    // **************  DEFINE SCALAR VARIABLES  ********************
    //
    set_scalar_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Density,"density",bc,FORT_DENFILL);
    desc_lst.setComponent(State_Type,Trac,"tracer",bc,FORT_ADVFILL);
    //
    // **************  DEFINE TEMPERATURE  ********************
    //
    if (do_temp)
    {
        set_temp_bc(bc,phys_bc);
        desc_lst.setComponent(State_Type,Temp,"temp",bc,FORT_TEMPFILL);
    }

    is_conservative.resize(NUM_STATE);
    is_diffusive.resize(NUM_STATE);
    for (int i = 0; i < NUM_STATE; i++)
    {
        is_conservative[i] = false;
        is_diffusive[i] = false;
        if (visc_coef[i] > 0.0)
            is_diffusive[i] = true;
    }
    is_conservative[Density] = true;
    if (do_temp) is_conservative[Temp] = false;
    is_conservative[Trac] = false;
    if (is_diffusive[Density])
    {
        BoxLib::Error("Density cannot diffuse, bad visc_coef");
    }
    //
    // ---- pressure
    //
    desc_lst.addDescriptor(Press_Type,IndexType::TheNodeType(),
                           StateDescriptor::Interval,1,1,
                           &node_bilinear_interp);
    set_pressure_bc(bc,phys_bc);
    desc_lst.setComponent(Press_Type,Pressure,"pressure",bc,FORT_PRESFILL);

    if (do_temp)
    {
	// stick Divu_Type on the end of the descriptor list
	Divu_Type = desc_lst.length();
	int nGrowDivu = 1;
	desc_lst.addDescriptor(Divu_Type,IndexType::TheCellType(),
                               StateDescriptor::Point,nGrowDivu,1,
			       &cell_cons_interp);
	set_divu_bc(bc,phys_bc);
	desc_lst.setComponent(Divu_Type,Divu,"divu",bc,FORT_DIVUFILL);
	
	// stick Dsdt_Type on the end of the descriptor list
	Dsdt_Type = desc_lst.length();
	int nGrowDsdt = 0;
	desc_lst.addDescriptor(Dsdt_Type,IndexType::TheCellType(),
                               StateDescriptor::Point,nGrowDsdt,1,
			       &cell_cons_interp);
	set_dsdt_bc(bc,phys_bc);
	desc_lst.setComponent(Dsdt_Type,Dsdt,"dsdt",bc,FORT_DSDTFILL);
    }
    //
    // **************  DEFINE DERIVED QUANTITIES ********************
    //
    // kinetic energy
    //
    derive_lst.add("energy",IndexType::TheCellType(),1,FORT_DERKENG,the_same_box);
    derive_lst.addComponent("energy",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("energy",desc_lst,State_Type,Xvel,BL_SPACEDIM);

    derive_lst.add("mag_vel",IndexType::TheCellType(),1,FORT_DERMVEL,the_same_box);
    derive_lst.addComponent("mag_vel",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // magnitude of vorticity
    //
    derive_lst.add("mag_vort",IndexType::TheCellType(),1,FORT_DERMGVORT,grow_box_by_one);
    derive_lst.addComponent("mag_vort",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // divergence of velocity field
    //
    derive_lst.add("diveru",IndexType::TheCellType(),1,FORT_DERMGDIVU,grow_box_by_one);
    derive_lst.addComponent("diveru",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // average pressure
    //
    derive_lst.add("avg_pressure",IndexType::TheCellType(),1,FORT_DERAVGPRES,
                   the_same_box);
    derive_lst.addComponent("avg_pressure",desc_lst,Press_Type,Pressure,1);
    //
    // pressure gradient in X direction
    //
    derive_lst.add("gradpx",IndexType::TheCellType(),1,FORT_DERGRDPX,the_same_box);
    derive_lst.addComponent("gradpx",desc_lst,Press_Type,Pressure,1);
    //
    // pressure gradient in Y direction
    //
    derive_lst.add("gradpy",IndexType::TheCellType(),1,FORT_DERGRDPY,the_same_box);
    derive_lst.addComponent("gradpy",desc_lst,Press_Type,Pressure,1);
    //
    // magnitude of pressure gradient 
    //
    derive_lst.add("gradp",IndexType::TheCellType(),1,FORT_DERGRDP,the_same_box);
    derive_lst.addComponent("gradp",desc_lst,Press_Type,Pressure,1);

#if (BL_SPACEDIM == 3)
    //
    // pressure gradient in Z direction
    //
    derive_lst.add("gradpz",IndexType::TheCellType(),1,FORT_DERGRDPZ,the_same_box);
    derive_lst.addComponent("gradpz",desc_lst,Press_Type,Pressure,1);
#endif
    //
    // **************  DEFINE ERROR ESTIMATION QUANTITIES  *************
    //
    err_list.add("tracer",1,ErrorRec::Special,FORT_ADVERROR);
    err_list.add("mag_vort",0,ErrorRec::Special,FORT_MVERROR);
}

void
NavierStokes::sum_integrated_quantities ()
{
    int finest_level = parent->finestLevel();
    Real time        = state[State_Type].curTime();
    Real mass        = 0.0;
    Real trac        = 0.0;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        NavierStokes& ns_level = getLevel(lev);
        mass += ns_level.volWgtSum("density",time);
        trac += ns_level.volWgtSum("tracer",time);
    }

    if (ParallelDescriptor::IOProcessor())
    {
        int old_prec = cout.precision(12);
        cout << '\n';
        cout << "TIME= " << time << " MASS= " << mass << '\n';
        cout << "TIME= " << time << " TRAC= " << trac << '\n';
        cout.precision(old_prec);
    }
}

