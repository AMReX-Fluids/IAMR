//BL_COPYRIGHT_NOTICE


#ifndef _MACOPDRIVERS_H_

//====================================================================
// This header file contains declarations for the MacOperator Projection
// driver routines
//====================================================================

// scale the variables for a projection solve
void proj_scale_var( MultiFab *rho,         MultiFab *vel,
                     const BoxArray &grids, int level );

// unscale the variables for a projection solve
void proj_unscale_var( MultiFab *rho,         MultiFab *vel,
                       const BoxArray &grids, int level );

void ProjFirstGuess( MultiFab &U_new, MultiFab &P_new,
                     int level, const BoxArray &grids );
#endif
