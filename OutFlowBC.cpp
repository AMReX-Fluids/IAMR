//
// $Id: OutFlowBC.cpp,v 1.6 2003-02-19 19:13:09 almgren Exp $
//
#include <winstd.H>

#include <algorithm>

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
    return BoxLib::grow(baseBox,grow_factor);
}

Box
OutFlowBC::SemiCoarsen (const Box& baseBox,
                        int        ref_factor,
                        int        direction)
{
    IntVect ref_ratio(D_DECL(ref_factor,ref_factor,ref_factor));
    ref_ratio[direction] = 1;
    return BoxLib::coarsen(baseBox,ref_ratio);
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
    Real goal  = std::max(rlast*tolerance,abs_tolerance);

    if (verbose)
    {
        std::cout << "OutFlowBC:Initial Residual: " << rlast << std::endl;
    }
    if (rlast > goal)
    {
        while (((res = vcycle(i1,i2)) > goal) && (iter < maxIters))
        {
            iter++;
            if (verbose)
                std::cout << "OutFlowBC: Residual: "
                          << res
                          << " at iteration "
                          << iter << std::endl;
        }
    }
  
    if (iter >= maxIters)
    {
        std::cout << "OutFlowBC: solver reached maxIter\n"
                  << "goal was: " << goal << " && res = " << res << std::endl;
    }

    if (verbose)
    {
        std::cout << "OutFlowBC: Final Residual: " << res << " after " 
                  << iter << " cycles\n\n";
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
        Restrict();
        next->phi->setVal(0);
        next->vcycle(downiter,upiter);
        interpolate();
        step(upiter);
    } 
    
    return rnorm;
}

void
OutFlowBC::GetOutFlowFaces (bool&        haveOutFlow,
                            Orientation* outFaces,
                            BCRec*       _phys_bc,
                            int&        numOutFlowBC)
{
    haveOutFlow = false;

    numOutFlowBC = 0;

    for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
        if (_phys_bc->lo(idir) == Outflow)
        {
            haveOutFlow = true;
            outFaces[numOutFlowBC] = Orientation(idir,Orientation::low);
            numOutFlowBC++;
        }

        if (_phys_bc->hi(idir) == Outflow)
        {
            haveOutFlow = true;
            outFaces[numOutFlowBC] = Orientation(idir,Orientation::high);
            numOutFlowBC++;
        }

    }

#if 0
//  We used to constrain our outflow bc to only one wall.  We are now
//  in the process of fixing that (10/2/02)
    if (numOutFlowBC > 1)
        //
        // True signals low-D solve for outflow.
        // False will enforce Div(U) == 0.
        //
        haveOutFlow = false;
#endif
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

#if 0
//  We used to constrain our outflow bc to only one wall.  We are now
//  in the process of fixing that (10/2/02)
    if (numOutFlowBC > 1)
        //
        // True signals low-D solve for outflow.
        // False will enforce Div(U) == 0.
        //
        has_out_flow = false;
#endif

    return has_out_flow;
}
