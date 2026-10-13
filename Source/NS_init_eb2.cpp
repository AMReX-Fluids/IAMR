#include <NavierStokesBase.H>
using namespace amrex;

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <EBUserDefined.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include "AMReX_VisMF.H"
#include "AMReX_PlotFileUtil.H"



#if (AMREX_SPACEDIM == 3)
static
void reentrant_profile(std::vector<amrex::RealVect> &points) {
  amrex::RealVect p;

  p = amrex::RealVect(AMREX_D_DECL(36.193*0.1, 7.8583*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(35.924*0.1, 7.7881*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(35.713*0.1, 7.5773*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(35.643*0.1, 7.3083*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(35.3*0.1, 7.0281*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(35.421*0.1, 6.241*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(34.82*0.1, 5.686*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(30.539*0.1, 3.5043*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(29.677*0.1, 2.6577*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(29.457*0.1, 1.47*0.1, 0.0));
  points.push_back(p);
  // p = amrex::RealVect(AMREX_D_DECL(29.38*0.1, -1.1038*0.1, 0.0));
  // points.push_back(p);
  // p = amrex::RealVect(AMREX_D_DECL(29.3*0.1, -2.7262*0.1, 0.0));
  // points.push_back(p);
  // p = amrex::RealVect(AMREX_D_DECL(29.273*0.1, -4.3428*0.1, 0.0));
  // points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(28.364*0.1, -5.7632*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(27.151*0.1, -6.8407*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(25.694*0.1, -7.5555*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(24.035*0.1, -7.8586*0.1, 0.0));
  points.push_back(p);
  p = amrex::RealVect(AMREX_D_DECL(22.358*0.1, -7.6902*0.1, 0.0));
  points.push_back(p);
}
#endif

// called in main before Amr->init(start,stop)
void
initialize_EB2 (const Geometry& geom, int required_coarsening_level,
                int max_coarsening_level)
{
    // read in EB parameters
    ParmParse ppeb2("eb2");
    std::string geom_type;
    ppeb2.get("geom_type", geom_type);

  if (geom_type == "combustor")
  {
#if (AMREX_SPACEDIM == 2)
    Abort("geom_type 'combustor' only available in 3D");
#elif (AMREX_SPACEDIM == 3)
    ParmParse pp("combustor");

    Real fwl;
    pp.get("far_wall_loc",fwl);

    EB2::PlaneIF farwall({AMREX_D_DECL(fwl,0.,0.)},
                         {AMREX_D_DECL(1. ,0.,0.)});

    Vector<Real> pl1pt, pl2pt, pl2nm, pl3pt;
    pp.getarr("ramp_plane1_point", pl1pt);
    pp.getarr("ramp_plane2_point", pl2pt);
    pp.getarr("ramp_plane2_normal", pl2nm);
    pp.getarr("ramp_plane3_point", pl3pt);

    auto ramp = EB2::makeIntersection(EB2::PlaneIF({pl1pt[0], pl1pt[1], 0.},
                                                   {      0.,      -1., 0.}),
                                      EB2::PlaneIF({pl2pt[0], pl2pt[1], 0.},
                                                   {pl2nm[0], pl2nm[1], 0.}),
                                      EB2::PlaneIF({pl3pt[0], pl3pt[1], 0.},
                                                   {      1.,       0., 0.}));

    Vector<Real> pipelo, pipehi;
    pp.getarr("pipe_lo", pipelo);
    pp.getarr("pipe_hi", pipehi);

    EB2::BoxIF pipe({pipelo[0], pipelo[1], -1.}, {pipehi[0], pipehi[1], 1.}, false);

    // where does plane 1 and plane 2 intersect?
    Real k2 = std::abs(pl2nm[0]/pl2nm[1]);
    Real secty = pl2pt[1] + k2*(pl3pt[0]-pl2pt[0]);
    // How much do we cut?
    Real dx = geom.CellSize(0);
    Real dycut = 4.*(1.+max_coarsening_level)*std::min(dx, k2*dx);
    EB2::BoxIF flat_corner({pl3pt[0], 0., -1.}, {1.e10, secty+dycut, 1.}, false);

    auto polys = EB2::makeUnion(farwall, ramp, pipe, flat_corner);

    // Real lenx = Geometry::ProbLength(0);
    // Real leny = Geometry::ProbLength(1);
    Real lenx = DefaultGeometry().ProbLength(0);
    Real leny = DefaultGeometry().ProbLength(1);
    auto pr = EB2::translate(EB2::lathe(polys), {lenx*0.5, leny*0.5, 0.});

    auto gshop = EB2::makeShop(pr);
    EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
#endif
  }
  else if (geom_type == "Piston-Cylinder")
  {
#if (AMREX_SPACEDIM == 2)
    Abort("geom_type 'Piston-Cylinder' only available in 3D");
#elif (AMREX_SPACEDIM == 3)
    EB2::SplineIF Piston;

    std::vector<amrex::RealVect> splpts;
    reentrant_profile(splpts);
    Piston.addSplineElement(splpts);

    amrex::RealVect p;
    std::vector<amrex::RealVect> lnpts;

    p = amrex::RealVect(AMREX_D_DECL(22.358*0.1, -7.6902*0.1, 0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(AMREX_D_DECL(1.9934*0.1, 3.464*0.1, 0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(AMREX_D_DECL(0.0*0.1, 3.464*0.1, 0.0));
    lnpts.push_back(p);
    Piston.addLineElement(lnpts);
    lnpts.clear();

    p = amrex::RealVect(AMREX_D_DECL(49.0*0.1, 7.8583*0.1,  0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(AMREX_D_DECL(36.193*0.1, 7.8583*0.1, 0.0));
    lnpts.push_back(p);
    Piston.addLineElement(lnpts);
    lnpts.clear();

    EB2::CylinderIF cylinder(48.0*0.1, 70.0*0.1, 2, {0.0, 0.0, -10.0*0.1}, true);

    auto revolvePiston  = EB2::lathe(Piston);
    //auto PistonComplement = EB2::makeComplement(revolvePiston);
    //auto PistonCylinder = EB2::makeIntersection(revolvePiston, cylinder);
    auto PistonCylinder = EB2::makeUnion(revolvePiston, cylinder);
    auto gshop = EB2::makeShop(PistonCylinder);
    EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
#endif
  }
  else if (geom_type == "Line-Piston-Cylinder")
  {
#if (AMREX_SPACEDIM == 2)
    Abort("geom_type 'Line-Piston-Cylinder' only available in 3D");
#elif (AMREX_SPACEDIM == 3)
    EB2::SplineIF Piston;
    std::vector<amrex::RealVect> lnpts;
    amrex::RealVect p;

    Real scaleFact;
    scaleFact = 0.0025; //MKS

    p = amrex::RealVect(AMREX_D_DECL(49.0*0.1*scaleFact, 7.8583*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(AMREX_D_DECL(36.193*0.1*scaleFact, 7.8583*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    Piston.addLineElement(lnpts);
    lnpts.clear();

    p = amrex::RealVect(AMREX_D_DECL(36.193*0.1*scaleFact, 7.8583*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(AMREX_D_DECL(24.035*0.1*scaleFact, -7.8586*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    Piston.addLineElement(lnpts);
    lnpts.clear();

    p = amrex::RealVect(AMREX_D_DECL(24.035*0.1*scaleFact, -7.8586*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(AMREX_D_DECL(20.0*0.1*scaleFact, -7.8586*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    Piston.addLineElement(lnpts);
    lnpts.clear();

    p = amrex::RealVect(AMREX_D_DECL(20.0*0.1*scaleFact, -7.8586*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(AMREX_D_DECL(1.9934*0.1*scaleFact, 3.464*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    Piston.addLineElement(lnpts);
    lnpts.clear();

    p = amrex::RealVect(AMREX_D_DECL(1.9934*0.1*scaleFact, 3.464*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(AMREX_D_DECL(0.09061*0.1*scaleFact, 3.464*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    Piston.addLineElement(lnpts);

    EB2::CylinderIF cylinder(48.0*0.1*scaleFact, 70.0*0.1*scaleFact, 2, {0.0, 0.0, -10.0*0.1*scaleFact}, true);

    auto revolvePiston  = EB2::lathe(Piston);
    auto PistonCylinder = EB2::makeUnion(revolvePiston, cylinder);
    auto gshop = EB2::makeShop(PistonCylinder);
    EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
#endif
  }
  else if (geom_type == "Inflow-Pipe")
  {
#if (AMREX_SPACEDIM == 2)
    Abort("geom_type 'Inflow-Pipe' only available in 3D");
#elif (AMREX_SPACEDIM == 3)
    // Initialise parameters
    int direction1 = 2;
    int direction2 = 2;
    Real radius1 = 0.018;
    Real radius2 = 0.007;
    Real height1 = 0.01;
    Real height2 = 0.01;
    Vector<Real> centervec1(3);
    Vector<Real> centervec2(3);

    // Get information from inputs file.
    ParmParse pp("pipe");

    pp.query("direction1", direction1);
    pp.query("direction2", direction2);
    pp.query("radius1", radius1);
    pp.query("radius2", radius2);
    pp.query("height1", height1);
    pp.query("height2", height2);
    pp.getarr("center1", centervec1, 0, 3);
    pp.getarr("center2", centervec2, 0, 3);
    Array<Real, 3> center1 = {centervec1[0], centervec1[1], centervec1[2]};
    Array<Real, 3> center2 = {centervec2[0], centervec2[1], centervec2[2]};

    // Compute distance between cylinder centres
    Real offset = 0.0;
    for(int i = 0; i < 3; i++)
        offset += pow(center1[i] - center2[i], 2);
    offset = sqrt(offset);

    // Print info about cylinders
    amrex::Print() << " CYLINDER 1" << '\n';
    amrex::Print() << " Direction:       " << direction1 << '\n';
    amrex::Print() << " Radius:    " << radius1 << '\n';
    amrex::Print() << " Center:    "
                   << center1[0] << ", " << center1[1] << ", " << center1[2] << '\n';

    amrex::Print() << " CYLINDER 2" << '\n';
    amrex::Print() << " Direction:       " << direction2 << '\n';
    amrex::Print() << " Radius:    " << radius2 << '\n';
    amrex::Print() << " Center:    "
                   << center2[0] << ", " << center2[1] << ", " << center2[2] << '\n';

    amrex::Print() << "\n Offset:          " << offset << '\n';

        // Build the implicit function as a union of two cylinders
    EB2::CylinderIF cyl1(radius1, height1, direction1, center1, false);
    EB2::CylinderIF cyl2(radius2, height2, direction2, center2, false);

    auto twocylinders = EB2::makeDifference(cyl1, cyl2);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(twocylinders);

    // Build index space
    EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
#endif
  }
  else if (geom_type == "Mixing-Pipe")
  {
#if (AMREX_SPACEDIM == 2)
    Abort("geom_type 'Mixing-Pipe' only available in 3D");
#elif (AMREX_SPACEDIM == 3)
    // Initialise parameters
    int direction = 1;
    Real radius = 0.018;
    Real height = 0.01;
    bool internal_flow = true;
    Vector<Real> centervec(3);

    // Get information from inputs file.
    ParmParse pp("pipe");

    pp.query("direction", direction);
    pp.query("radius", radius);
    pp.query("height", height);
    pp.getarr("center", centervec, 0, 3);
    pp.query("internal_flow", internal_flow);
    Array<Real, 3> center = {centervec[0], centervec[1], centervec[2]};

    // Print info about cylinders
    amrex::Print() << " CYLINDER " << '\n';
    amrex::Print() << " Direction:       " << direction << '\n';
    amrex::Print() << " Radius:    " << radius << '\n';
    amrex::Print() << " Center:    "
                   << center[0] << ", " << center[1] << ", " << center[2] << '\n';



        // Build the implicit function as a union of two cylinders
    EB2::CylinderIF cyl(radius, height, direction, center, internal_flow);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(cyl);
    // Build index space
    EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
#endif
  }
  else if (geom_type == "Square-Grid")
  {
#if (AMREX_SPACEDIM == 2)
    Abort("geom_type 'Square-Grid' only available in 3D");
#elif (AMREX_SPACEDIM == 3)
    // Initialise parameters
    Real dim_L0 = 0.08;
    Real ratio_t0_L0_cross = 0.11;
    Real ratio_t0_stream_thickness = 1.0;

    // Get information from inputs file.
    ParmParse pp("square_grid");

    pp.query("dim_L0", dim_L0);
    pp.query("ratio_t0_L0_cross", ratio_t0_L0_cross);
    pp.query("ratio_t0_stream_thickness", ratio_t0_stream_thickness);


    Real cross_dim_t0 = ratio_t0_L0_cross * dim_L0;
    Real pos_big_square = 0.5 * (dim_L0 + cross_dim_t0);
    Real pos_small_square = 0.5 * (dim_L0 - cross_dim_t0);
    Real stream_length = cross_dim_t0 * ratio_t0_stream_thickness;

    Array<Real, 3> big_square_lo = {0.0, -pos_big_square,-pos_big_square};
    Array<Real, 3> big_square_hi = {stream_length,  pos_big_square, pos_big_square};

    Array<Real, 3> small_square_lo = {0.0, -pos_small_square,-pos_small_square};
    Array<Real, 3> small_square_hi = {stream_length,  pos_small_square, pos_small_square};


    // Print info about the square grid parameters
    amrex::Print() << " SQUARE GRID PARAMETERS " << '\n';
    amrex::Print() << " dim_L0:       " << dim_L0 << '\n';
    amrex::Print() << " computed cross section dim_t0:       " << cross_dim_t0 << '\n';
    amrex::Print() << " computed streamwise section length:       " << stream_length << '\n';
    amrex::Print() << " ratio_t0_L0_cross:    " << ratio_t0_L0_cross << '\n';
    amrex::Print() << " ratio_t0_stream_thickness:    " << ratio_t0_stream_thickness << '\n';
    amrex::Print() << " pos_big_square:    " << pos_big_square << '\n';
    amrex::Print() << " pos_small_square:    " << pos_small_square << '\n';


        // Build the implicit function as a union of two cylinders
    EB2::BoxIF big_square(big_square_lo, big_square_hi,   false);
    EB2::BoxIF small_square(small_square_lo, small_square_hi, false);
    auto square_grid = EB2::makeDifference(big_square, small_square);


    // Generate GeometryShop
    auto gshop = EB2::makeShop(square_grid);
    // Build index space
    EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
#endif
  }
  else if (geom_type == "UserDefined") {
      EBUserDefined(geom, required_coarsening_level, max_coarsening_level);
  }
  else
  {
    EB2::Build(geom, required_coarsening_level, max_coarsening_level);
  }
}

#endif
