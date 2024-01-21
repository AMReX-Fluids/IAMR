#include <RegType.H>
#include <AMReX_BC_TYPES.H>

//
// Components are  Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
//
static int norm_vel_bc[] =
{
    amrex::BCType::int_dir, amrex::BCType::ext_dir, amrex::BCType::foextrap, amrex::BCType::reflect_odd, amrex::BCType::ext_dir, amrex::BCType::ext_dir
};

static int tang_vel_bc[] =
{
    amrex::BCType::int_dir, amrex::BCType::ext_dir, amrex::BCType::foextrap, amrex::BCType::reflect_even, amrex::BCType::hoextrap, amrex::BCType::ext_dir
};

static int scalar_bc[] =
{
    amrex::BCType::int_dir, amrex::BCType::ext_dir, amrex::BCType::foextrap, amrex::BCType::reflect_even, amrex::BCType::foextrap, amrex::BCType::foextrap
};

static int press_bc[] =
{
    amrex::BCType::int_dir, amrex::BCType::foextrap, amrex::BCType::foextrap, amrex::BCType::reflect_even, amrex::BCType::foextrap, amrex::BCType::foextrap
};

static int norm_gradp_bc[] =
{
    amrex::BCType::int_dir, amrex::BCType::foextrap, amrex::BCType::foextrap, amrex::BCType::reflect_odd, amrex::BCType::foextrap, amrex::BCType::foextrap
};

static int tang_gradp_bc[] =
{
    amrex::BCType::int_dir, amrex::BCType::foextrap, amrex::BCType::foextrap, amrex::BCType::reflect_even, amrex::BCType::foextrap, amrex::BCType::foextrap
};

static int temp_bc[] =
{
    amrex::BCType::int_dir, amrex::BCType::ext_dir, amrex::BCType::hoextrap, amrex::BCType::reflect_even, amrex::BCType::reflect_even, amrex::BCType::reflect_even
};

static int divu_bc[] =
{
    amrex::BCType::int_dir, amrex::BCType::reflect_even, amrex::BCType::reflect_even, amrex::BCType::reflect_even, amrex::BCType::reflect_even, amrex::BCType::reflect_even
};

static int dsdt_bc[] =
{
    amrex::BCType::int_dir, amrex::BCType::ext_dir, amrex::BCType::ext_dir, amrex::BCType::reflect_even, amrex::BCType::reflect_even, amrex::BCType::reflect_even
};

static int average_bc[] =
{
    amrex::BCType::int_dir, amrex::BCType::int_dir, amrex::BCType::int_dir, amrex::BCType::int_dir,      amrex::BCType::int_dir,      amrex::BCType::int_dir
};
