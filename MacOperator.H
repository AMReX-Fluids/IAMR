#ifndef _MACOPERATOR_H_
#define _MACOPERATOR_H_

//
// $Id: MacOperator.H,v 1.10 2006-09-08 21:21:02 almgren Exp $
//

//==========================================================================
// This class defines how mass continuity is defined for the projection operators
// in class PorousMedia.  For the Mac project, this is done through setting a
// linop and a few driver functions, while for the mac project this is done
// through a few driver functions which scale variables.
//==========================================================================

#include <ABecLaplacian.H>
#include <BndryData.H>

class MacOperator
    :
    public ABecLaplacian
{
public:

    MacOperator (const BndryData& mgb,
                 const Real*      _h);

    virtual ~MacOperator ();
    //
    // Define the gradient coeficients used by the multigrid object.
    //
    void setCoefficients (MultiFab*   area,
                          MultiFab&   coef,
                          const Real* dx);
    //
    // Initialize the rhs for the level solves.
    //
    void defRHS (MultiFab* area,
                 MultiFab& vol,
                 MultiFab& rhs,
                 MultiFab* vel,
                 Real      scale);
    //
    // Apply the pressure gradient to the mac velocities.
    //
    void velUpdate (MultiFab*       vel,
                    MultiFab&       phi,
                    const MultiFab& coef,
                    const Real*     dx,
                    Real            scale);
    //
    // Scale the rhs for the sync solves.
    //
    void syncRhs (const MultiFab& Volume,
                  MultiFab&       Rhs,
                  Real            rhs_scale,
                  const Real*     dx);
};

#endif
