

#include <NavierStokes.H>
#include <NS_BC.H>
#include <RegType.H>
#include <NS_derive.H>
#include <NS_bcfill.H>
#include <AMReX_FArrayBox.H>

using namespace amrex;

static Box the_same_box (const Box& b)    { return b;                 }
static Box grow_box_by_one (const Box& b) { return amrex::grow(b,1); }
static Box grow_box_by_two (const Box& b) { return amrex::grow(b,2); }

// NOTE: the int arrays that define the mapping from physical BCs to mathematical
// (norm_vel_bc, tang_vel_bc, scalar_bc, temp_bc, press_bc, divu_bc, dsdt_bc)
// are now all defined in IAMR/Source/NS_BC.H

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
set_gradpx_bc (BCRec&       bc,
	       const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,norm_gradp_bc[lo_bc[0]]);
    bc.setHi(0,norm_gradp_bc[hi_bc[0]]);
    bc.setLo(1,tang_gradp_bc[lo_bc[1]]);
    bc.setHi(1,tang_gradp_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_gradp_bc[lo_bc[2]]);
    bc.setHi(2,tang_gradp_bc[hi_bc[2]]);
#endif
}

static
void
set_gradpy_bc (BCRec&       bc,
	       const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_gradp_bc[lo_bc[0]]);
    bc.setHi(0,tang_gradp_bc[hi_bc[0]]);
    bc.setLo(1,norm_gradp_bc[lo_bc[1]]);
    bc.setHi(1,norm_gradp_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_gradp_bc[lo_bc[2]]);
    bc.setHi(2,tang_gradp_bc[hi_bc[2]]);
#endif
}

#if (BL_SPACEDIM == 3)
static
void
set_gradpz_bc (BCRec&       bc,
	       const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_gradp_bc[lo_bc[0]]);
    bc.setHi(0,tang_gradp_bc[hi_bc[0]]);
    bc.setLo(1,tang_gradp_bc[lo_bc[1]]);
    bc.setHi(1,tang_gradp_bc[hi_bc[1]]);
    bc.setLo(2,norm_gradp_bc[lo_bc[2]]);
    bc.setHi(2,norm_gradp_bc[hi_bc[2]]);
}
#endif

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

static
void
set_average_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i,average_bc[lo_bc[i]]);
        bc.setHi(i,average_bc[hi_bc[i]]);
    }
}


typedef StateDescriptor::BndryFunc BndryFunc;

//
// Get EB-aware interpolater when needed
//
#ifdef AMREX_USE_EB
  static auto& cc_interp = eb_cell_cons_interp;
#else
  static auto& cc_interp = cell_cons_interp;
#endif

void
NavierStokes::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    Initialize();

    BCRec bc;
    //
    // Set number of state variables.
    //
    NUM_STATE = Density + 1;
    Tracer = NUM_STATE++;
    if (do_trac2)
        Tracer2 = NUM_STATE++;
    if (do_temp)
        Temp = NUM_STATE++;
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
    bool state_data_extrap = false;
    bool store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
    			   StateDescriptor::Point,NUM_GROW,NUM_STATE,
    			   &cc_interp,state_data_extrap,store_in_checkpoint);
    // TODO: state does not really need to carry ghost cells, since it is the
    // philosophy of IAMR to FillPatch before using.
    // However, changing NUM_GROW here creates a problem for restarting from
    // older checkpoint files with more ghost cells, so a workaround is
    // needed.

    BndryFunc vel_bf(vel_fill);
    vel_bf.setRunOnGPU(true);

    BndryFunc state_bf(state_fill);
    state_bf.setRunOnGPU(true);

    BndryFunc press_bf(press_fill);
    press_bf.setRunOnGPU(true);

    BndryFunc null_bf(dummy_fill);
    null_bf.setRunOnGPU(true);

    set_x_vel_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Xvel,"x_velocity",bc,vel_bf);

    set_y_vel_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Yvel,"y_velocity",bc,vel_bf);

#if (BL_SPACEDIM == 3)
    set_z_vel_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Zvel,"z_velocity",bc,vel_bf);
#endif
    //
    // **************  DEFINE SCALAR VARIABLES  ********************
    //
    set_scalar_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Density,"density",bc,state_bf);

    set_scalar_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Tracer,"tracer",bc,state_bf);

    if (do_trac2)
    {
       set_scalar_bc(bc,phys_bc);
       desc_lst.setComponent(State_Type,Tracer2,"tracer2",bc,state_bf);
    }
    //
    // **************  DEFINE TEMPERATURE  ********************
    //
    if (do_temp)
    {
        set_temp_bc(bc,phys_bc);
        desc_lst.setComponent(State_Type,Temp,"temp",bc,state_bf);
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

    advectionType[Tracer] = NonConservative;
    diffusionType[Tracer] = Laplacian_S; if (do_cons_trac) {
      advectionType[Tracer] = Conservative;
      diffusionType[Tracer] = Laplacian_SoverRho;
      amrex::Print() << "Using conservative advection update for tracer.\n";
    }

    if (do_trac2) {
	advectionType[Tracer2] = NonConservative;
	diffusionType[Tracer2] = Laplacian_S;
	if (do_cons_trac2) {
	  advectionType[Tracer2] = Conservative;
	  diffusionType[Tracer2] = Laplacian_SoverRho;
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
    desc_lst.addDescriptor(Press_Type,IndexType::TheNodeType(),
                           StateDescriptor::Interval,1,1,
                           &node_bilinear_interp);

    set_pressure_bc(bc,phys_bc);
    desc_lst.setComponent(Press_Type,Pressure,"pressure",bc,press_bf);

    //
    // ---- grad P
    //
    //
    // FIXME ----
    // maybe we're better off recomputing rather than reading from chk? But then still have
    // FillPatch issue at coarse fine boundary; if want to interpolate in time, need an old and new
    desc_lst.addDescriptor(Gradp_Type,IndexType::TheCellType(),
    			   StateDescriptor::Interval,gradp_grow,AMREX_SPACEDIM,
    			   &cc_interp,state_data_extrap,store_in_checkpoint);

    Vector<BCRec>       bcs(BL_SPACEDIM);
    Vector<std::string> name(BL_SPACEDIM);

    set_gradpx_bc(bc,phys_bc);
    bcs[0]  = bc;
    name[0] = "gradpx";

    set_gradpy_bc(bc,phys_bc);
    bcs[1]  = bc;
    name[1] = "gradpy";

#if(AMREX_SPACEDIM==3)
    set_gradpz_bc(bc,phys_bc);
    bcs[2]  = bc;
    name[2] = "gradpz";
#endif

    desc_lst.setComponent(Gradp_Type, Gradpx, name, bcs, null_bf);

    //
    // ---- Additions for using Temperature
    //
    if (do_temp)
    {
	// stick Divu_Type on the end of the descriptor list
	Divu_Type = desc_lst.size();
	int nGrowDivu = 1;
	desc_lst.addDescriptor(Divu_Type,IndexType::TheCellType(),
                               StateDescriptor::Point,nGrowDivu,1,
			       &cc_interp);
	set_divu_bc(bc,phys_bc);
	desc_lst.setComponent(Divu_Type,Divu,"divu",bc,null_bf);

	// stick Dsdt_Type on the end of the descriptor list
	Dsdt_Type = desc_lst.size();
	int nGrowDsdt = 0;
	desc_lst.addDescriptor(Dsdt_Type,IndexType::TheCellType(),
                               StateDescriptor::Point,nGrowDsdt,1,
			       &cc_interp);
	set_dsdt_bc(bc,phys_bc);
	desc_lst.setComponent(Dsdt_Type,Dsdt,"dsdt",bc,null_bf);
    }

    //
    // For using on-the-fly averaging
    //
    if (NavierStokesBase::avg_interval > 0)
    {
      Average_Type = desc_lst.size();
      bool state_data_extrap = false;
      bool store_in_checkpoint = true;
      desc_lst.addDescriptor(Average_Type,IndexType::TheCellType(),
                             StateDescriptor::Point,0,BL_SPACEDIM*2,
                             &cc_interp,state_data_extrap,store_in_checkpoint);

      set_average_bc(bc,phys_bc);
      desc_lst.setComponent(Average_Type,Xvel,"xvel_avg_dummy",bc,null_bf);
      desc_lst.setComponent(Average_Type,Xvel+BL_SPACEDIM,"xvel_rms_dummy",bc,null_bf);
      desc_lst.setComponent(Average_Type,Yvel,"yvel_avg_dummy",bc,null_bf);
      desc_lst.setComponent(Average_Type,Yvel+BL_SPACEDIM,"yvel_rms_dummy",bc,null_bf);
#if (BL_SPACEDIM==3)
      desc_lst.setComponent(Average_Type,Zvel,"zvel_avg_dummy",bc,null_bf);
      desc_lst.setComponent(Average_Type,Zvel+BL_SPACEDIM,"zvel_rms_dummy",bc,null_bf);
#endif
    }

    //
    // **************  DEFINE DERIVED QUANTITIES ********************
    //
    using namespace derive_functions;

    if (NavierStokesBase::avg_interval > 0)
    {
      //
      // Average and RMS velocity
      //
      Vector<std::string> var_names_ave(BL_SPACEDIM*2);
      var_names_ave[Xvel] = "x_vel_average";
      var_names_ave[Yvel] = "y_vel_average";
#if (BL_SPACEDIM==3)
      var_names_ave[Zvel] = "z_vel_average";
#endif
      var_names_ave[Xvel+BL_SPACEDIM] = "x_vel_rms";
      var_names_ave[Yvel+BL_SPACEDIM] = "y_vel_rms";
#if (BL_SPACEDIM==3)
      var_names_ave[Zvel+BL_SPACEDIM] = "z_vel_rms";
#endif
      derive_lst.add("velocity_average",IndexType::TheCellType(),BL_SPACEDIM*2,
                     var_names_ave,der_vel_avg,the_same_box);
      derive_lst.addComponent("velocity_average",desc_lst,Average_Type,Xvel,BL_SPACEDIM*2);
    }

    //
    // kinetic energy
    //
    derive_lst.add("energy",IndexType::TheCellType(),1,derkeng,the_same_box);
    derive_lst.addComponent("energy",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("energy",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // magnitude of vorticity
    //
    derive_lst.add("mag_vort",IndexType::TheCellType(),1,dermgvort,grow_box_by_two);
    derive_lst.addComponent("mag_vort",desc_lst,State_Type,Xvel,BL_SPACEDIM);
    //
    // average pressure
    //
    derive_lst.add("avg_pressure",IndexType::TheCellType(),1,deravgpres,
                   the_same_box);
    derive_lst.addComponent("avg_pressure",desc_lst,Press_Type,Pressure,1);

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
