#ifndef _MACOPERATOR_H_
#define _MACOPERATOR_H_

//
// $Id: MacOperator.H,v 1.3 1997-09-26 16:57:01 lijewski Exp $
//

//==========================================================================
// This class defines how mass continuity is defined for the projection operators
// in class NavierStokes.  For the Mac project, this is done through setting a
// linop and a few driver functions, while for the mac project this is done
// through a few driver functions which scale variables.
//==========================================================================

#include <ABecLaplacian.H>
#include <BndryData.H>

class MacOperator : public ABecLaplacian
{
public:
    MacOperator(const BoxArray &ba, const BndryData& mgb, const REAL * _h);

    virtual ~MacOperator() {}

    // define the gradient coeficients used by the multigrid object.
    void setCoefficients(MultiFab* area, MultiFab& rho, int rho_comp,
			 const REAL* dx);

    // initialize the rhs for the level solves
    void defRHS(MultiFab* area, MultiFab& vol, MultiFab& rhs,
                MultiFab *vel, REAL scale);

    // apply the pressure gradient to the mac velocities
    void velUpdate(MultiFab* vel, MultiFab& phi, const MultiFab& Rho,
                   int rho_comp, const REAL* dx, REAL scale);

    // scale the rhs for the sync solves
    void syncRhs(const MultiFab& Volume, MultiFab& Rhs,
                 REAL rhs_scale, const REAL* dx);
};


#endif
