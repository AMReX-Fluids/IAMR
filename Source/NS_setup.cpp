

#include <NavierStokes.H>
#include <NS_BC.H>
#include <RegType.H>
#include <AMReX_ErrorList.H>
#include <PROB_NS_F.H>
#include <DERIVE_F.H>
#include <NS_derive.H>
#include <NS_error_F.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

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
    desc_lst.addDescriptor(Press_Type,IndexType::TheNodeType(),
                           StateDescriptor::Interval,1,1,
                           &node_bilinear_interp);

    set_pressure_bc(bc,phys_bc);
    desc_lst.setComponent(Press_Type,Pressure,"pressure",bc,BndryFunc(FORT_PRESFILL));
 
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
    amrex::StateDescriptor::BndryFunc gradp_bf(dummy_fill);
    gradp_bf.setRunOnGPU(true);

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

    desc_lst.setComponent(Gradp_Type, Gradpx, name, bcs, gradp_bf);

/*
    if (NavierStokesBase::avg_interval > 0)
    {
      Average_Type = desc_lst.size();
      bool state_data_extrap = false;
      bool store_in_checkpoint = true;
      desc_lst.addDescriptor(Average_Type,IndexType::TheCellType(),
                             StateDescriptor::Point,0,BL_SPACEDIM*2,
                             &cc_interp,state_data_extrap,store_in_checkpoint);

      set_average_bc(bc,phys_bc);
      desc_lst.setComponent(Average_Type,Xvel,"xvel_avg_dummy",bc,BndryFunc(FORT_DSDTFILL));
      desc_lst.setComponent(Average_Type,Xvel+BL_SPACEDIM,"xvel_rms_dummy",bc,BndryFunc(FORT_DSDTFILL));
      desc_lst.setComponent(Average_Type,Yvel,"yvel_avg_dummy",bc,BndryFunc(FORT_DSDTFILL));
      desc_lst.setComponent(Average_Type,Yvel+BL_SPACEDIM,"yvel_rms_dummy",bc,BndryFunc(FORT_DSDTFILL));
#if (BL_SPACEDIM==3)
      desc_lst.setComponent(Average_Type,Zvel,"zvel_avg_dummy",bc,BndryFunc(FORT_DSDTFILL));
      desc_lst.setComponent(Average_Type,Zvel+BL_SPACEDIM,"zvel_rms_dummy",bc,BndryFunc(FORT_DSDTFILL));
#endif

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
*/
      
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
	desc_lst.setComponent(Divu_Type,Divu,"divu",bc,BndryFunc(FORT_DIVUFILL));
	
	// stick Dsdt_Type on the end of the descriptor list
	Dsdt_Type = desc_lst.size();
	int nGrowDsdt = 0;
	desc_lst.addDescriptor(Dsdt_Type,IndexType::TheCellType(),
                               StateDescriptor::Point,nGrowDsdt,1,
			       &cc_interp);
	set_dsdt_bc(bc,phys_bc);
	desc_lst.setComponent(Dsdt_Type,Dsdt,"dsdt",bc,BndryFunc(FORT_DSDTFILL));
    }

    if (NavierStokesBase::avg_interval > 0)
    {
      Average_Type = desc_lst.size();
      bool state_data_extrap = false;
      bool store_in_checkpoint = true;
      desc_lst.addDescriptor(Average_Type,IndexType::TheCellType(),
                             StateDescriptor::Point,0,BL_SPACEDIM*2,
                             &cc_interp,state_data_extrap,store_in_checkpoint);

      set_average_bc(bc,phys_bc);
      desc_lst.setComponent(Average_Type,Xvel,"xvel_avg_dummy",bc,BndryFunc(FORT_DSDTFILL));
      desc_lst.setComponent(Average_Type,Xvel+BL_SPACEDIM,"xvel_rms_dummy",bc,BndryFunc(FORT_DSDTFILL));
      desc_lst.setComponent(Average_Type,Yvel,"yvel_avg_dummy",bc,BndryFunc(FORT_DSDTFILL));
      desc_lst.setComponent(Average_Type,Yvel+BL_SPACEDIM,"yvel_rms_dummy",bc,BndryFunc(FORT_DSDTFILL));
#if (BL_SPACEDIM==3)
      desc_lst.setComponent(Average_Type,Zvel,"zvel_avg_dummy",bc,BndryFunc(FORT_DSDTFILL));
      desc_lst.setComponent(Average_Type,Zvel+BL_SPACEDIM,"zvel_rms_dummy",bc,BndryFunc(FORT_DSDTFILL));
#endif

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
    // magnitude of pressure gradient 
    //
    derive_lst.add("gradp",IndexType::TheCellType(),1,dergrdp,the_same_box);
    derive_lst.addComponent("gradp",desc_lst,Press_Type,Pressure,1);

#if (BL_SPACEDIM == 3)
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
    // Pressure stuff for on-the-fly integration
    //
    derive_lst.add("PresVars",IndexType::TheCellType(),4,derpresvars,the_same_box);
    derive_lst.addComponent("PresVars",desc_lst,Press_Type,Pressure,1);
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

    //
    // Dynamically generated error tagging functions
    //
    std::string amr_prefix = "amr";
    ParmParse ppamr(amr_prefix);
    Vector<std::string> refinement_indicators;
    ppamr.queryarr("refinement_indicators",refinement_indicators,0,ppamr.countval("refinement_indicators"));
    for (int i=0; i<refinement_indicators.size(); ++i)
    {
        std::string ref_prefix = amr_prefix + "." + refinement_indicators[i];

        ParmParse ppr(ref_prefix);
        RealBox realbox;
        if (ppr.countval("in_box_lo")) {
            std::vector<Real> box_lo(BL_SPACEDIM), box_hi(BL_SPACEDIM);
            ppr.getarr("in_box_lo",box_lo,0,box_lo.size());
            ppr.getarr("in_box_hi",box_hi,0,box_hi.size());
            realbox = RealBox(&(box_lo[0]),&(box_hi[0]));
        }

        AMRErrorTagInfo info;

        if (realbox.ok()) {
            info.SetRealBox(realbox);
        }
        if (ppr.countval("start_time") > 0) {
            Real min_time; ppr.get("start_time",min_time);
            info.SetMinTime(min_time);
        }
        if (ppr.countval("end_time") > 0) {
            Real max_time; ppr.get("end_time",max_time);
            info.SetMaxTime(max_time);
        }
        if (ppr.countval("max_level") > 0) {
            int max_level; ppr.get("max_level",max_level);
            info.SetMaxLevel(max_level);
        }

        if (ppr.countval("value_greater")) {
            Real value; ppr.get("value_greater",value);
            std::string field; ppr.get("field_name",field);
            errtags.push_back(AMRErrorTag(value,AMRErrorTag::GREATER,field,info));
        }
        else if (ppr.countval("value_less")) {
            Real value; ppr.get("value_less",value);
            std::string field; ppr.get("field_name",field);
            errtags.push_back(AMRErrorTag(value,AMRErrorTag::LESS,field,info));
        }
        else if (ppr.countval("vorticity_greater")) {
            Real value; ppr.get("vorticity_greater",value);
            const std::string field="mag_vort";
            errtags.push_back(AMRErrorTag(value,AMRErrorTag::VORT,field,info));
        }
        else if (ppr.countval("adjacent_difference_greater")) {
            Real value; ppr.get("adjacent_difference_greater",value);
            std::string field; ppr.get("field_name",field);
            errtags.push_back(AMRErrorTag(value,AMRErrorTag::GRAD,field,info));
        }
        else if (realbox.ok())
        {
            errtags.push_back(AMRErrorTag(info));
        }
        else {
            Abort(std::string("Unrecognized refinement indicator for " + refinement_indicators[i]).c_str());
        }
    }

    error_setup();
}
