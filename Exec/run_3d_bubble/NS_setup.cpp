

#include <NavierStokes.H>
#include <NS_BC.H>
#include <RegType.H>
#include <AMReX_ErrorList.H>
#include <PROB_NS_F.H>
#include <DERIVE_F.H>
#include <AMReX_FArrayBox.H>

using namespace amrex;

static Box the_same_box (const Box& b)    { return b;                 }
static Box grow_box_by_one (const Box& b) { return amrex::grow(b,1); }

// NOTE: the int arrays norm_vel_bc, tang_vel_bc, scalar_bc, temp_bc, press_bc, divu_bc, dsdt_bc 
//                      are now all defined in NS_BC.H in iamrlib

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

typedef StateDescriptor::BndryFunc BndryFunc;

void
NavierStokes::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        phys_bc.setLo(dir,SlipWall);
        phys_bc.setHi(dir,SlipWall);
    }

    Initialize();

    BCRec bc;
    //
    // Set number of state variables.
    //
    NUM_STATE = Density + 1;
    int Trac = NUM_STATE++;
    NUM_SCALARS = NUM_STATE - Density;

    if (do_scalar_update_in_order) {
	// Need to check numbers and values of scalar update
	// Idea is to specify the scalar index (counting Density as zero)
	int maxComp=NUM_SCALARS-1;
	for (int iComp=0; iComp<maxComp; iComp++) {
	    if ((scalarUpdateOrder[iComp]>maxComp)||(scalarUpdateOrder[iComp]<1))
		amrex::Abort("Scalar Update Order out of bounds");
	    for (int jComp=iComp+1; jComp<maxComp; jComp++)
		if (scalarUpdateOrder[iComp]==scalarUpdateOrder[jComp])
		    amrex::Abort("Scalar Update Order values not unique");
	}
    }

    //
    // **************  DEFINE VELOCITY VARIABLES  ********************
    //
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,1,NUM_STATE,
                           &cell_cons_interp);
    set_x_vel_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Xvel,"x_velocity",bc,BndryFunc(FORT_XVELFILL));

    set_y_vel_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Yvel,"y_velocity",bc,BndryFunc(FORT_YVELFILL));

#if (BL_SPACEDIM == 3)
    set_z_vel_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Zvel,"z_velocity",bc,BndryFunc(FORT_ZVELFILL));
#endif
    //
    // **************  DEFINE SCALAR VARIABLES  ********************
    //
    set_scalar_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Density,"density",bc,BndryFunc(FORT_DENFILL));

    set_scalar_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Trac,"tracer",bc,BndryFunc(FORT_ADVFILL));

    is_diffusive.resize(NUM_STATE);
    advectionType.resize(NUM_STATE);
    diffusionType.resize(NUM_STATE);
    for (int i = 0; i < NUM_STATE; i++)
    {
        advectionType[i] = NonConservative;
        diffusionType[i] = RhoInverse_Laplacian_S;
        is_diffusive[i] = false;
        if (visc_coef[i] > 0.0)
            is_diffusive[i] = true;
    }

    if (do_mom_diff == 1)
      for (int d = 0; d < BL_SPACEDIM; d++)
        advectionType[Xvel+d] = Conservative;

    advectionType[Density] = Conservative;

    advectionType[Trac] = NonConservative;
    diffusionType[Trac] = Laplacian_S; if (do_cons_trac) {
      advectionType[Trac] = Conservative;
      diffusionType[Trac] = Laplacian_SoverRho;
      if (ParallelDescriptor::IOProcessor())
	std::cout << "Using conservative advection update for tracer." << std::endl;
    }

    if (is_diffusive[Density])
    {
        amrex::Error("Density cannot diffuse, bad visc_coef");
    }
    //
    // ---- pressure
    //
    desc_lst.addDescriptor(Press_Type,IndexType::TheNodeType(),
                           StateDescriptor::Interval,1,1,
                           &node_bilinear_interp);

    set_pressure_bc(bc,phys_bc);
    desc_lst.setComponent(Press_Type,Pressure,"pressure",bc,BndryFunc(FORT_PRESFILL));

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
#if (BL_SPACEDIM == 3)
    //
    //  vorticity vector field
    //
    derive_lst.add("vort_x",IndexType::TheCellType(),1,FORT_DERVORTX,grow_box_by_one);
    derive_lst.addComponent("vort_x",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    derive_lst.add("vort_y",IndexType::TheCellType(),1,FORT_DERVORTY,grow_box_by_one);
    derive_lst.addComponent("vort_y",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    derive_lst.add("vort_z",IndexType::TheCellType(),1,FORT_DERVORTZ,grow_box_by_one);
    derive_lst.addComponent("vort_z",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    derive_lst.add("DMagVort",IndexType::TheCellType(),1,FORT_DERDMAG,grow_box_by_one);
    derive_lst.addComponent("DMagVort",desc_lst,State_Type,Xvel,BL_SPACEDIM);
#endif
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

#ifdef AMREX_PARTICLES
    //
    // The particle count at this level.
    //
    derive_lst.add("particle_count",IndexType::TheCellType(),1,
                   FORT_DERNULL,the_same_box);
    derive_lst.addComponent("particle_count",desc_lst,State_Type,Density,1);
    //
    // The total # of particles at our level or above.
    //
    derive_lst.add("total_particle_count",IndexType::TheCellType(),1,
                   FORT_DERNULL,the_same_box);
    derive_lst.addComponent("total_particle_count",desc_lst,State_Type,Density,1);
#endif

    //
    // **************  DEFINE ERROR ESTIMATION QUANTITIES  *************
    //
    error_setup();
}
