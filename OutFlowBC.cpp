//
// $Id: OutFlowBC.cpp,v 1.2 2000-11-02 23:11:17 lijewski Exp $
//

#include "OutFlowBC.H"
#include "ParmParse.H"

OutFlowBC::OutFlowBC () {}

OutFlowBC::~OutFlowBC () {}

Box
OutFlowBC::SemiGrow (const Box& baseBox,
                         int        nGrow,
                         int        direction)
{
    IntVect grow_factor(D_DECL(nGrow,nGrow,nGrow));
    grow_factor[direction] = 0;
    return ::grow(baseBox,grow_factor);
}

Box
OutFlowBC::SemiCoarsen (const Box& baseBox,
                        int        ref_factor,
                        int        direction)
{
    IntVect ref_ratio(D_DECL(ref_factor,ref_factor,ref_factor));
    ref_ratio[direction] = 1;
    return ::coarsen(baseBox,ref_ratio);
}

OutFlowBC_MG::OutFlowBC_MG (const Box& Domain,
                            FArrayBox* Phi,
                            FArrayBox* Rhs,
                            FArrayBox* Resid,
                            FArrayBox* Beta,
                            Real*      H,
                            int*       IsPeriodic,
                            bool       is_scalar)
    :
    domain(Domain),
    phi(Phi),
    rhs(Rhs),
    resid(Resid),
    cgwork(0),
    beta(Beta),
    next(0),
    beta_is_scalar(is_scalar)
{
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        h[dir]          = H[dir];
        isPeriodic[dir] = IsPeriodic[dir];
    }
}

OutFlowBC_MG::~OutFlowBC_MG ()
{
    if (next != 0)
    {
        delete next->phi;
        delete next->rhs;
        delete next->resid;
        if (beta_is_scalar)
            delete next->beta;
        else
            delete [] next->beta;
        delete next;
    }
    delete cgwork;
}

void 
OutFlowBC_MG::solve (Real tolerance,
                     Real abs_tolerance,
                     int  i1,
                     int  i2,
                     int  maxIters,
                     int  verbose)
{
    int  iter  = 1;
    Real rlast = residual();
    Real res   = rlast;
    Real goal  = Max(rlast*tolerance,abs_tolerance);

    if (verbose)
    {
        cout << "OutFlowBC:Initial Residual: " << rlast << endl;
    }
    if (rlast > goal)
    {
        while (((res = vcycle(i1,i2)) > goal) && (iter < maxIters))
        {
            iter++;
            if (verbose)
                cout << "OutFlowBC: Residual: " << res << " at iteration " << iter << endl;
        }
    }
  
    if (iter >= maxIters)
    {
        cout << "OutFlowBC: solver reached maxIter" << endl;
        cout << "goal was: " << goal << " && res = " << res << endl;
    }

    if (verbose)
    {
        cout << "OutFlowBC: Final Residual: " << res << " after " 
             << iter << " cycles" << endl;
        cout << " " << endl;
    }
}

Real 
OutFlowBC_MG::vcycle (int downiter,
                      int upiter)
{
    Real rnorm = residual();
    step(downiter);
    rnorm = residual();

    if (next != 0)
    {
        restrict();
        next->phi->setVal(0);
        next->vcycle(downiter,upiter);
        interpolate();
        step(upiter);
    } 
    
    return rnorm;
}

void
OutFlowBC::GetOutFlowFace (bool&        haveOutFlow,
                           Orientation& outFace,
                           BCRec*       _phys_bc)
{
    haveOutFlow = false;

    int numOutFlowBC = 0;

    for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
        if (_phys_bc->lo(idir) == Outflow)
        {
            haveOutFlow = true;
            outFace = Orientation(idir,Orientation::low);
            numOutFlowBC++;
        }

        if (_phys_bc->hi(idir) == Outflow)
        {
            haveOutFlow = true;
            outFace = Orientation(idir,Orientation::high);
            numOutFlowBC++;
        }

    }

    if (numOutFlowBC > 1)
        //
        // True signals low-D solve for outflow.
        // False will enforce Div(U) == 0.
        //
        haveOutFlow = false;
}

bool
OutFlowBC::HasOutFlowBC (BCRec* _phys_bc)
{
    bool has_out_flow = false;
    int  numOutFlowBC = 0;

    for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
        if (_phys_bc->lo(idir) == Outflow)
        {
            has_out_flow = true;
            numOutFlowBC++;
        }

        if (_phys_bc->hi(idir) == Outflow)
        {
            has_out_flow = true;
            numOutFlowBC++;
        }
    }

    if (numOutFlowBC > 1)
        //
        // True signals low-D solve for outflow.
        // False will enforce Div(U) == 0.
        //
        has_out_flow = false;

    return has_out_flow;
}
