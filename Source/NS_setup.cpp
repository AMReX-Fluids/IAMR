

#include <NavierStokes.H>
#include <NS_BC.H>
#include <RegType.H>
#include <AMReX_ErrorList.H>
#include <PROB_NS_F.H>
#include <DERIVE_F.H>
#include <NS_error_F.H>
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
    int Trac2;
    if (do_trac2)
	Trac2 = NUM_STATE++;
    if (do_temp) NUM_STATE++;
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
    // FIXME? - cribbed from CNS
#ifdef AMREX_USE_EB
    bool state_data_extrap = false;
    bool store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
    			   StateDescriptor::Point,NUM_GROW,NUM_STATE,
    			   &eb_cell_cons_interp,state_data_extrap,store_in_checkpoint);
#else
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,1,NUM_STATE,
                           &cell_cons_interp);
#endif
    
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

    if (do_trac2)
    {
       set_scalar_bc(bc,phys_bc);
       desc_lst.setComponent(State_Type,Trac2,"tracer2",bc,BndryFunc(FORT_ADV2FILL));
    }
    //
    // **************  DEFINE TEMPERATURE  ********************
    //
    if (do_temp)
    {
        set_temp_bc(bc,phys_bc);
        desc_lst.setComponent(State_Type,Temp,"temp",bc,BndryFunc(FORT_TEMPFILL));
    }

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
    if (do_temp) advectionType[Temp] = NonConservative;

    advectionType[Trac] = NonConservative;
    diffusionType[Trac] = Laplacian_S; if (do_cons_trac) {
      advectionType[Trac] = Conservative;
      diffusionType[Trac] = Laplacian_SoverRho;
      amrex::Print() << "Using conservative advection update for tracer.\n";
    }

    if (do_trac2) {
	advectionType[Trac2] = NonConservative;
	diffusionType[Trac2] = Laplacian_S;
	if (do_cons_trac2) {
	  advectionType[Trac2] = Conservative;
	  diffusionType[Trac2] = Laplacian_SoverRho;
	  amrex::Print() << "Using conservative advection update for tracer2.\n";
	}
    }

    if (is_diffusive[Density])
    {
        amrex::Error("Density cannot diffuse, bad visc_coef");
    }
    //
    // ---- pressure
    //
#if 1
    desc_lst.addDescriptor(Press_Type,IndexType::TheNodeType(),
                           StateDescriptor::Interval,1,1,
                           &node_bilinear_interp);

    set_pressure_bc(bc,phys_bc);
    desc_lst.setComponent(Press_Type,Pressure,"pressure",bc,BndryFunc(FORT_PRESFILL));
#else
    desc_lst.addDescriptor(Press_Type,IndexType::TheNodeType(),
                           StateDescriptor::Point,1,1,
                           &node_bilinear_interp,true);

    set_pressure_bc(bc,phys_bc);
    desc_lst.setComponent(Press_Type,Pressure,"pressure",bc,BndryFunc(FORT_PRESFILL));
    //
    // ---- time derivative of pressure
    //
    Dpdt_Type = desc_lst.length();
    desc_lst.addDescriptor(Dpdt_Type,IndexType::TheNodeType(),
                           StateDescriptor::Interval,1,1,
                           &node_bilinear_interp);
    set_pressure_bc(bc,phys_bc);
    desc_lst.setComponent(Dpdt_Type,Dpdt,"dpdt",bc,BndryFunc(FORT_PRESFILL));
#endif

    if (do_temp)
    {
	// stick Divu_Type on the end of the descriptor list
	Divu_Type = desc_lst.size();
	int nGrowDivu = 1;
	desc_lst.addDescriptor(Divu_Type,IndexType::TheCellType(),
                               StateDescriptor::Point,nGrowDivu,1,
			       &cell_cons_interp);
	set_divu_bc(bc,phys_bc);
	desc_lst.setComponent(Divu_Type,Divu,"divu",bc,BndryFunc(FORT_DIVUFILL));
	
	// stick Dsdt_Type on the end of the descriptor list
	Dsdt_Type = desc_lst.size();
	int nGrowDsdt = 0;
	desc_lst.addDescriptor(Dsdt_Type,IndexType::TheCellType(),
                               StateDescriptor::Point,nGrowDsdt,1,
			       &cell_cons_interp);
	set_dsdt_bc(bc,phys_bc);
	desc_lst.setComponent(Dsdt_Type,Dsdt,"dsdt",bc,BndryFunc(FORT_DSDTFILL));
    }
    //
    // **************  DEFINE DERIVED QUANTITIES ********************
    //
    // mod grad rho
    //
    derive_lst.add("modgradrho",IndexType::TheCellType(),1,dermodgradrho,grow_box_by_one);
    derive_lst.addComponent("modgradrho",desc_lst,State_Type,Density,1);

#if (BL_SPACEDIM==3)
    //
    // u dot laplacian u
    //
    derive_lst.add("udotlapu",IndexType::TheCellType(),1,derudotlapu,grow_box_by_one);
    derive_lst.addComponent("udotlapu",desc_lst,State_Type,Xvel,BL_SPACEDIM);
#endif
    //
    // kinetic energy
    //
    derive_lst.add("energy",IndexType::TheCellType(),1,derkeng,the_same_box);
    derive_lst.addComponent("energy",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("energy",desc_lst,State_Type,Xvel,BL_SPACEDIM);

    derive_lst.add("mag_vel",IndexType::TheCellType(),1,dermvel,the_same_box);
    derive_lst.addComponent("mag_vel",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // magnitude of vorticity
    //
    derive_lst.add("mag_vort",IndexType::TheCellType(),1,DeriveFunc3D(dermgvort),grow_box_by_one);
    derive_lst.addComponent("mag_vort",desc_lst,State_Type,Xvel,BL_SPACEDIM);
#if (BL_SPACEDIM == 3)
    //
    //  vorticity vector field
    //
    derive_lst.add("vort_x",IndexType::TheCellType(),1,dervortx,grow_box_by_one);
    derive_lst.addComponent("vort_x",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    derive_lst.add("vort_y",IndexType::TheCellType(),1,dervorty,grow_box_by_one);
    derive_lst.addComponent("vort_y",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    derive_lst.add("vort_z",IndexType::TheCellType(),1,dervortz,grow_box_by_one);
    derive_lst.addComponent("vort_z",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    derive_lst.add("DMagVort",IndexType::TheCellType(),1,derdmag,grow_box_by_one);
    derive_lst.addComponent("DMagVort",desc_lst,State_Type,Xvel,BL_SPACEDIM);
#endif
    //
    // divergence of velocity field
    //
    derive_lst.add("diveru",IndexType::TheCellType(),1,DeriveFunc3D(dermgdivu),grow_box_by_one);
    derive_lst.addComponent("diveru",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // average pressure
    //
    derive_lst.add("avg_pressure",IndexType::TheCellType(),1,DeriveFunc3D(deravgpres),
                   the_same_box);
    derive_lst.addComponent("avg_pressure",desc_lst,Press_Type,Pressure,1);
    //
    // pressure gradient in X direction
    //
    derive_lst.add("gradpx",IndexType::TheCellType(),1,DeriveFunc3D(dergrdpx),the_same_box);
    derive_lst.addComponent("gradpx",desc_lst,Press_Type,Pressure,1);
    //
    // pressure gradient in Y direction
    //
    derive_lst.add("gradpy",IndexType::TheCellType(),1,DeriveFunc3D(dergrdpy),the_same_box);
    derive_lst.addComponent("gradpy",desc_lst,Press_Type,Pressure,1);
    //
    // magnitude of pressure gradient 
    //
    derive_lst.add("gradp",IndexType::TheCellType(),1,dergrdp,the_same_box);
    derive_lst.addComponent("gradp",desc_lst,Press_Type,Pressure,1);

#if (BL_SPACEDIM == 3)
    //
    // pressure gradient in Z direction
    //
    derive_lst.add("gradpz",IndexType::TheCellType(),1,DeriveFunc3D(dergrdpz),the_same_box);
    derive_lst.addComponent("gradpz",desc_lst,Press_Type,Pressure,1);
    //
    // radial velocity
    //
    derive_lst.add("radial_velocity",IndexType::TheCellType(),1,derradvel,the_same_box);
    derive_lst.addComponent("radial_velocity",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // azimuthal velocity
    //
    derive_lst.add("azimuthal_velocity",IndexType::TheCellType(),1,derazivel,the_same_box);
    derive_lst.addComponent("azimuthal_velocity",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // x_velocity in laboratory frame for rotating frame of refernce
    //
    derive_lst.add("x_velocity_rot",IndexType::TheCellType(),1,derxvelrot,the_same_box);
    derive_lst.addComponent("x_velocity_rot",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // y_velocity in laboratory frame for rotating frame of refernce
    //
    derive_lst.add("y_velocity_rot",IndexType::TheCellType(),1,deryvelrot,the_same_box);
    derive_lst.addComponent("y_velocity_rot",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // magnitude of velocity in laboratory frame for rotating frame of refernce
    //
    derive_lst.add("mag_velocity_rot",IndexType::TheCellType(),1,dermagvelrot,the_same_box);
    derive_lst.addComponent("mag_velocity_rot",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // magnitude of vorticity in laboratory frame for rotating frame of refernce
    //
    derive_lst.add("mag_vorticity_rot",IndexType::TheCellType(),1,dermagvortrot,grow_box_by_one);
    derive_lst.addComponent("mag_vorticity_rot",desc_lst,State_Type,Xvel,BL_SPACEDIM);
#if defined(DO_IAMR_FORCE)
    //
    // forcing - used to calculate the rate of injection of energy in probtype 14 (HIT)
    //
    derive_lst.add("forcing",IndexType::TheCellType(),1,derforcing,the_same_box);
    derive_lst.addComponent("forcing",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("forcing",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // forcex - used to put the forcing term in the plot file
    //
    derive_lst.add("forcex",IndexType::TheCellType(),1,derforcex,the_same_box);
    derive_lst.addComponent("forcex",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("forcex",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // forcey - used to put the forcing term in the plot file
    //
    derive_lst.add("forcey",IndexType::TheCellType(),1,derforcey,the_same_box);
    derive_lst.addComponent("forcey",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("forcey",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // forcez - used to put the forcing term in the plot file
    //
    derive_lst.add("forcez",IndexType::TheCellType(),1,derforcez,the_same_box);
    derive_lst.addComponent("forcez",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("forcez",desc_lst,State_Type,Xvel,BL_SPACEDIM);
#endif
    //
    // Turbulence Variable - for integrating on the fly
    //
    derive_lst.add("TurbVars",IndexType::TheCellType(),16,derturbvars,grow_box_by_one);
    derive_lst.addComponent("TurbVars",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("TurbVars",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // Pressure stuff for on-the-fly integration
    //
    derive_lst.add("PresVars",IndexType::TheCellType(),4,derpresvars,the_same_box);
    derive_lst.addComponent("PresVars",desc_lst,Press_Type,Pressure,1);
#ifdef SUMJET
    //
    // Jet Variables - for integrating on the fly
    //
    derive_lst.add("JetVars",IndexType::TheCellType(),27,derjetvars,grow_box_by_one);
    derive_lst.addComponent("JetVars",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("JetVars",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    derive_lst.addComponent("JetVars",desc_lst,State_Type,Trac,1);
    if (do_trac2)
	derive_lst.addComponent("JetVars",desc_lst,State_Type,Trac2,1);
    //
    // Pressure stuff for the jet on-the-fly integration (note need to grow for slope reconstruction of *derivatives*)
    //
    derive_lst.add("JetPresVars",IndexType::TheCellType(),16,derjetpresvars,grow_box_by_one);
    derive_lst.addComponent("JetPresVars",desc_lst,Press_Type,Pressure,1); 
#endif
//3D
#endif

#ifdef AMREX_PARTICLES
    //
    // The particle count at this level.
    //
    derive_lst.add("particle_count",IndexType::TheCellType(),1,
                   dernull,the_same_box);
    derive_lst.addComponent("particle_count",desc_lst,State_Type,Density,1);
    //
    // The total # of particles at our level or above.
    //
    derive_lst.add("total_particle_count",IndexType::TheCellType(),1,
                   dernull,the_same_box);
    derive_lst.addComponent("total_particle_count",desc_lst,State_Type,Density,1);
#endif

    //
    // **************  DEFINE ERROR ESTIMATION QUANTITIES  *************
    //
    error_setup();
}
