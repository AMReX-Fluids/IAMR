#include <RegType.H>

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
    INT_DIR, EXT_DIR, HOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, FOEXTRAP
};

static int divu_bc[] =
{
    INT_DIR, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static int dsdt_bc[] =
{
    INT_DIR, EXT_DIR, EXT_DIR, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static RegType project_bc [] =
{
    interior, inflow, outflow, refWall, refWall, refWall
};
