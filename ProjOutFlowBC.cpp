//
// $Id: ProjOutFlowBC.cpp,v 1.16 2001-08-09 22:42:00 marc Exp $
//
#include <winstd.H>

#include "ProjOutFlowBC.H"
#include "PROJOUTFLOWBC_F.H"
#include "ParmParse.H"

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

#define DEF_BOX_LIMITS(box,boxlo,boxhi)      \
const int* boxlo = (box).loVect();           \
const int* boxhi = (box).hiVect();

Real    ProjOutFlowBC::tol     = 1.0e-10; 
Real    ProjOutFlowBC::abs_tol = 5.0e-10;

int  ProjOutFlowBC_MG::verbose           = 0;
bool ProjOutFlowBC_MG::useCGbottomSolver = true;
Real ProjOutFlowBC_MG::cg_tol            = 1.0e-2;
Real ProjOutFlowBC_MG::cg_abs_tol        = 5.0e-12;
Real ProjOutFlowBC_MG::cg_max_jump       = 10.0;
int  ProjOutFlowBC_MG::cg_maxiter        = 40;
int  ProjOutFlowBC_MG::maxIters          = 40;

static
Box
semiSurroundingNodes (const Box& baseBox,
                      int        direction)
{
    Box sBox = BoxLib::surroundingNodes(baseBox);
    sBox.growHi(direction,-1);
    return sBox;
}

ProjOutFlowBC::ProjOutFlowBC ()
{
    ParmParse pp("projoutflow");
    pp.query("tol",tol);
    pp.query("abs_tol",abs_tol);
    solver = (BL_SPACEDIM == 2) ? BC_BACK : BC_MG;
    int solver_type = -1;
    pp.query("solver_type",solver_type);
    if (solver_type != -1)
        solver = OutFlow_Solver_Type(solver_type);
}

void 
ProjOutFlowBC::computeBC (FArrayBox*         velFab,
                          FArrayBox&         divuFab,
                          FArrayBox&         rhoFab,
                          FArrayBox&         phiFab,
                          const Geometry&    geom, 
                          const Orientation& outFace)
{
    const Real* dx    = geom.CellSize();
    const Box& domain = geom.Domain();
    int face          = int(outFace);
    int outDir        = outFace.coordDir();

    int isPeriodic[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
        isPeriodic[dir] = geom.isPeriodic(dir);

    IntVect loFiltered, hiFiltered;
    int isPeriodicFiltered[BL_SPACEDIM];
    Real dxFiltered[BL_SPACEDIM];
    //
    // Filter out the direction we don't care about.
    //
    int ncStripWidth = 1;
    Box origBox = BoxLib::adjCell(domain,outFace,ncStripWidth);
    IntVect lo = origBox.smallEnd();
    IntVect hi = origBox.bigEnd();
    //
    // Rearrange the box, dx, and isPeriodic so that the dimension that is 1
    // is the last dimension.
    //
    int cnt = 0;
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        if (dir != outDir)
	{
            loFiltered[cnt] = lo[dir];
            hiFiltered[cnt] = hi[dir];
            dxFiltered[cnt] = dx[dir];
            isPeriodicFiltered[cnt] = isPeriodic[dir];
            cnt++;
	}
        else
        {
            loFiltered[BL_SPACEDIM-1] = lo[dir];
            hiFiltered[BL_SPACEDIM-1] = hi[dir];
            dxFiltered[BL_SPACEDIM-1] = dx[dir];
            isPeriodicFiltered[BL_SPACEDIM-1] = isPeriodic[dir];
	}
    }

    Box       faceBox(loFiltered,hiFiltered);
    FArrayBox divuExt(faceBox,1);
    FArrayBox rhoExt(faceBox,1);
    FArrayBox uExt(faceBox,BL_SPACEDIM-1);

#if (BL_SPACEDIM == 2)
    //
    // Make cc r (set = 1 if cartesian)
    //
    int R_DIR = 0;
    int Z_DIR = 1;
    int perpDir = (outDir == Z_DIR) ? R_DIR : Z_DIR;
    Box region = BoxLib::adjCell(domain,outFace,1);
    region.shift(outDir, (outFace.faceDir() == Orientation::high) ? -1 : 1);
    int r_lo = region.smallEnd(perpDir);
    int r_hi = region.bigEnd(perpDir);

    Array<Real> rcen(region.length(perpDir), 1.0);
    if (CoordSys::IsRZ() && perpDir == R_DIR) 
        geom.GetCellLoc(rcen, region, perpDir);
#else
    Array<Real> rcen;
    int r_lo = 0;
    int r_hi = 0;
#endif

    DEF_BOX_LIMITS(origBox,origLo,origHi);
    DEF_LIMITS(uExt, uEPtr, uElo, uEhi);
    DEF_LIMITS(divuExt,divuEPtr,divuElo,divuEhi);
    DEF_LIMITS(rhoExt ,rhoEPtr ,rhoElo , rhoEhi );
    DEF_LIMITS(divuFab,divuPtr,divulo,divuhi);
    DEF_LIMITS(rhoFab ,rhoPtr ,rholo , rhohi );
    DEF_LIMITS(*velFab , velPtr, vello, velhi);
    //
    // Extrapolate the velocities, divu, and rho to the outflow edge in
    // the shifted coordinate system (where the last dimension is 1).
    //
    int zeroIt;
    FORT_EXTRAP_PROJ(ARLIM(vello),  ARLIM(velhi), velPtr,
                     ARLIM(divulo), ARLIM(divuhi), divuPtr,
                     ARLIM(rholo),    ARLIM(rhohi),rhoPtr,
#if (BL_SPACEDIM == 2)
                     &r_lo,&r_hi,rcen.dataPtr(),
#endif
                     ARLIM(uElo),ARLIM(uEhi),uEPtr,
                     ARLIM(divuElo),ARLIM(divuEhi),divuEPtr,
                     ARLIM(rhoElo),ARLIM(rhoEhi),rhoEPtr,
                     origLo,origHi,&face,&zeroIt);

    if (zeroIt)
    {
        //
        // No perturbations, set homogeneous bc.
        //
        phiFab.setVal(0);
    }
    else if (solver == BC_MG)
    {
        Box nodeBox = semiSurroundingNodes(faceBox,BL_SPACEDIM-1);
        FArrayBox phiFiltered(nodeBox,1);

        DEF_LIMITS(phiFab,phiFabPtr,phiFab_lo,phiFab_hi);
        DEF_LIMITS(phiFiltered,phiFilteredPtr,phiFiltered_lo,phiFiltered_hi);
        int face = int(outFace);
      
        FORT_HG_SHIFT_PHI(ARLIM(phiFiltered_lo),ARLIM(phiFiltered_hi),
                          phiFilteredPtr,
                          ARLIM(phiFab_lo),ARLIM(phiFab_hi),phiFabPtr,
                          &face);
      
        FArrayBox rhs_temp, beta;
      
        computeCoefficients(rhs_temp,beta,uExt,divuExt,rhoExt,rcen,
                            r_lo,r_hi,faceBox,dxFiltered,isPeriodicFiltered);
        //
        // Need phi to have ghost cells.
        //
        Box phiGhostBox = OutFlowBC::SemiGrow(phiFiltered.box(),1,BL_SPACEDIM-1);
        FArrayBox phi(phiGhostBox,1);
        phi.setVal(0);
        phi.copy(phiFiltered);
      
        Box grownRhs = OutFlowBC::SemiGrow(rhs_temp.box(),1,BL_SPACEDIM-1);
        FArrayBox rhs(grownRhs,1);
        rhs.setVal(0);
        rhs.copy(rhs_temp);
        FArrayBox resid(rhs.box(),1);
        ProjOutFlowBC_MG proj_mg(faceBox,&phi,&rhs,&resid,&beta,
                                 dxFiltered,isPeriodicFiltered);

        proj_mg.solve(tol,abs_tol,2,2,proj_mg.MaxIters(),proj_mg.Verbose());
      
        DEF_LIMITS(phi,phiPtr,phi_lo,phi_hi);
        DEF_BOX_LIMITS(faceBox,lo,hi);
        //
        // Subtract the average phi.
        //
        FORT_HGSUBTRACTAVGPHI(ARLIM(phi_lo),ARLIM(phi_hi),phiPtr,
#if (BL_SPACEDIM == 2)
                              &r_lo,&r_hi,rcen.dataPtr(),
#endif
                              lo,hi,isPeriodicFiltered);
        //
        // Translate the solution back to the original coordinate system.
        //
        FORT_HG_RESHIFT_PHI(ARLIM(phiFab_lo),ARLIM(phiFab_hi),phiFabPtr,
                            ARLIM(phi_lo),ARLIM(phi_hi),phiPtr,&face);
      
#if (BL_SPACEDIM == 2)  
    }
    else if (solver == BC_BACK)
    {
        solveBackSubstitution(phiFab,divuExt,uExt,rhoExt,rcen,r_lo,r_hi,
                              isPeriodic,dxFiltered,faceBox,outFace);
#endif
    }
    else
    {
        BoxLib::Error("ProjOutFlowBC::unknown solver_type");
    }
}

#if (BL_SPACEDIM == 2)
void
ProjOutFlowBC::solveBackSubstitution (FArrayBox&         phi,
                                      FArrayBox&         divuExt,
                                      FArrayBox&         uExt,
                                      FArrayBox&         rhoExt,
                                      Array<Real>&       rcen,
                                      int                r_lo,
                                      int                r_hi,
                                      int*               isPeriodicFiltered,
                                      Real*              dxFiltered,
                                      Box&               faceBox,
                                      const Orientation& outFace)
{
    int face   = int(outFace);
    int outDir = outFace.coordDir();
    int length = divuExt.length()[0];

    BL_ASSERT(length == uExt.length()[0]);
    BL_ASSERT(length == rhoExt.length()[0]);
    BL_ASSERT(length == rcen.size());

    DEF_BOX_LIMITS(faceBox,faceLo,faceHi);
    DEF_LIMITS(uExt, uEPtr, uElo, uEhi);
    DEF_LIMITS(divuExt,divuEPtr,divuElo,divuEhi);
    DEF_LIMITS(rhoExt ,rhoEPtr ,rhoElo , rhoEhi );
    DEF_LIMITS(phi, phiPtr,philo,phihi);

    FORT_HGPHIBC(dxFiltered,rcen.dataPtr(),
                 uEPtr,divuEPtr,rhoEPtr,&length,
                 ARLIM(philo), ARLIM(phihi),phiPtr, 
                 faceLo,faceHi,&face,isPeriodicFiltered);
}
#endif

void 
ProjOutFlowBC::computeCoefficients (FArrayBox&   rhs,
                                    FArrayBox&   beta,
                                    FArrayBox&   uExt,
                                    FArrayBox&   divuExt,
                                    FArrayBox&   rhoExt,
                                    Array<Real>& rcen,
                                    int          r_lo,
                                    int          r_hi,
                                    Box&         faceBox,
                                    Real*        dxFiltered,
                                    int*         isPeriodicFiltered)
{
    Box rhsBox  = semiSurroundingNodes(faceBox,BL_SPACEDIM-1);
    Box betaBox = OutFlowBC::SemiGrow(faceBox,1,BL_SPACEDIM-1);
  
    beta.resize(betaBox,1);
    rhs.resize(rhsBox,1);

    DEF_BOX_LIMITS(faceBox,faceLo,faceHi);
    DEF_LIMITS(beta,  betaPtr, betalo,betahi);
    DEF_LIMITS(rhs, rhsPtr, rhslo,rhshi);
    DEF_LIMITS(uExt, uEPtr, uElo, uEhi);
    DEF_LIMITS(divuExt,divuEPtr,divuElo,divuEhi);
    DEF_LIMITS(rhoExt ,rhoEPtr ,rhoElo , rhoEhi );

    FORT_COMPUTE_COEFF(ARLIM(rhslo),ARLIM(rhshi),rhsPtr,
                       ARLIM(betalo),ARLIM(betahi),betaPtr,
                       ARLIM(uElo), ARLIM(uEhi), uEPtr,
                       ARLIM(divuElo),ARLIM(divuEhi), divuEPtr,
                       ARLIM(rhoElo),ARLIM(rhoEhi),rhoEPtr,
#if (BL_SPACEDIM == 2)
                       &r_lo,&r_hi,rcen.dataPtr(),
#endif
                       faceLo,faceHi,
                       dxFiltered,isPeriodicFiltered);
}

ProjOutFlowBC_MG::ProjOutFlowBC_MG (const Box& Domain,
                                    FArrayBox* Phi,
                                    FArrayBox* Rhs,
                                    FArrayBox* Resid,
                                    FArrayBox* Beta,
                                    Real*      H,
                                    int*       IsPeriodic)
    :
    OutFlowBC_MG(Domain,Phi,Rhs,Resid,Beta,H,IsPeriodic,true)
{
    static int first = true;

    if (first)
    {
        first = false;

        ParmParse pp("proj_mg");

        pp.query("v",verbose);
        int use_cg;
        pp.query("useCGbottomSolver",use_cg);

        useCGbottomSolver = (use_cg > 0) ? true : false;

        pp.query("cg_tol",cg_tol);
        pp.query("cg_abs_tol",cg_abs_tol);
        pp.query("cg_max_jump",cg_max_jump);
        pp.query("cg_maxiter",cg_maxiter);
        pp.query("maxIters",maxIters);
    }

    const IntVect& len = domain.length();

    int min_length = 4;
    bool test_side[BL_SPACEDIM-1];
    for (int dir = 0; dir < BL_SPACEDIM-1; dir++)
        test_side[dir] = (len[dir]&1) != 0 || len[dir] < min_length;

    if (D_TERM(1 && ,test_side[0], || test_side[1]))
    {
        if (useCGbottomSolver)
        {
            Box temp1Box = OutFlowBC::SemiGrow(domain,1,BL_SPACEDIM-1);
            Box temp2Box = semiSurroundingNodes(temp1Box,BL_SPACEDIM-1);
            cgwork = new FArrayBox(temp2Box,4);
        }
    }
    else
    {
        Real newh[BL_SPACEDIM];
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            newh[dir] = 2*h[dir];

        Box newdomain = OutFlowBC::SemiCoarsen(domain,2,BL_SPACEDIM-1);
        Box grownBox  = OutFlowBC::SemiGrow(newdomain,1,BL_SPACEDIM-1);
        Box nodes     = semiSurroundingNodes(newdomain,BL_SPACEDIM-1);
        Box newp_size = OutFlowBC::SemiGrow(nodes,1,BL_SPACEDIM-1);

        FArrayBox* newphi    = new FArrayBox(newp_size,1);
        FArrayBox* newresid  = new FArrayBox(newp_size,1);
        FArrayBox* newrhs    = new FArrayBox(newp_size,1);
        FArrayBox* newbeta   = new FArrayBox(grownBox,1);
        newphi->setVal(0);
        newresid->setVal(0);
        newbeta->setVal(0);
   
        DEF_BOX_LIMITS(domain,dom_lo,dom_hi);
        DEF_BOX_LIMITS(newdomain,new_lo,new_hi);
        DEF_LIMITS(*beta,betaPtr,beta_lo,beta_hi);
        DEF_LIMITS(*newbeta,newbetaPtr,newbeta_lo,newbeta_hi);

        FORT_COARSIG(betaPtr,ARLIM(beta_lo),ARLIM(beta_hi),
                     newbetaPtr,ARLIM(newbeta_lo),ARLIM(newbeta_hi),
                     dom_lo,dom_hi,new_lo,new_hi,isPeriodic);

        next = new ProjOutFlowBC_MG(newdomain,newphi,newrhs,newresid, 
                                    newbeta,newh,isPeriodic);
    }
}

ProjOutFlowBC_MG::~ProjOutFlowBC_MG () {}


Real
ProjOutFlowBC_MG::residual ()
{
    Real rnorm;

    FArrayBox dgphi(semiSurroundingNodes(domain,BL_SPACEDIM-1),1);

    DEF_BOX_LIMITS(domain,lo,hi);
    DEF_LIMITS(*rhs,rhsPtr,rhslo,rhshi);
    DEF_LIMITS(*beta,betaPtr,betalo,betahi);
    DEF_LIMITS(*phi,phiPtr,philo,phihi);
    DEF_LIMITS(*resid,residPtr,residlo,residhi);
    DEF_LIMITS(dgphi,dgphiPtr,dglo,dghi);

    resid->setVal(0);

    FORT_HGRESID (ARLIM(rhslo), ARLIM(rhshi),  rhsPtr,
                  ARLIM(betalo), ARLIM(betahi),  betaPtr,
                  ARLIM(philo), ARLIM(phihi),  phiPtr,
                  ARLIM(residlo), ARLIM(residhi),  residPtr,
                  ARLIM(dglo), ARLIM(dghi), dgphiPtr,
                  lo, hi, h, isPeriodic, &rnorm);

    return rnorm;
}

void 
ProjOutFlowBC_MG::step (int nGSRB)
{
    if (cgwork != 0)
    {
        Real resnorm  = 0.0;

        FArrayBox dest0(phi->box(),1);

        DEF_BOX_LIMITS(domain,lo,hi);
        DEF_LIMITS(*phi,phiPtr,phi_lo,phi_hi);
        DEF_LIMITS(*resid,residPtr,resid_lo,resid_hi);
        DEF_LIMITS(dest0,dest0Ptr,dest0_lo,dest0_hi);
        DEF_LIMITS(*rhs,rhsPtr,rhs_lo,rhs_hi);
        DEF_LIMITS(*beta, betaPtr, beta_lo,beta_hi); 
        DEF_LIMITS(*cgwork,dummPtr,cg_lo,cg_hi);

        FORT_SOLVEHG(phiPtr,ARLIM(phi_lo),ARLIM(phi_hi),
                     dest0Ptr, ARLIM(dest0_lo),ARLIM(dest0_hi),
                     rhsPtr,ARLIM(rhs_lo),ARLIM(rhs_hi),
                     betaPtr, ARLIM(beta_lo),ARLIM(beta_hi),
                     cgwork->dataPtr(0),ARLIM(cg_lo),ARLIM(cg_hi),
                     cgwork->dataPtr(1), ARLIM(cg_lo),ARLIM(cg_hi),
                     cgwork->dataPtr(2), ARLIM(cg_lo),ARLIM(cg_hi),
                     cgwork->dataPtr(3),ARLIM(cg_lo),ARLIM(cg_hi),
                     residPtr, ARLIM(resid_lo),ARLIM(resid_hi),
                     lo,hi,h,isPeriodic,&cg_maxiter,&cg_tol,
                     &cg_abs_tol,&cg_max_jump,&resnorm);
    }
    else
    {
        gsrb(nGSRB);
    }
}

void 
ProjOutFlowBC_MG::restrict ()
{
    DEF_BOX_LIMITS(domain,lo,hi);
    DEF_BOX_LIMITS(next->theDomain(),loc,hic);
    DEF_LIMITS(*resid,residPtr,resid_lo,resid_hi);
    DEF_LIMITS(*(next->theRhs()),rescPtr,resc_lo,resc_hi);

    next->theRhs()->setVal(0);
  
    FORT_RESTRICT(residPtr,ARLIM(resid_lo),ARLIM(resid_hi), 
                  rescPtr, ARLIM(resc_lo),ARLIM(resc_hi),
                  lo,hi,loc,hic,isPeriodic);
}

void
ProjOutFlowBC_MG::interpolate ()
{
    FArrayBox temp(phi->box(),1);

    temp.setVal(0);
    
    DEF_BOX_LIMITS(domain,lo,hi);
    DEF_BOX_LIMITS(next->theDomain(),loc,hic);
    DEF_LIMITS(*phi,phiPtr,phi_lo,phi_hi);
    DEF_LIMITS(temp,tempPtr,temp_lo,temp_hi);
    DEF_LIMITS(*(next->thePhi()),deltacPtr,deltac_lo,deltac_hi);
    DEF_LIMITS(*beta,betaPtr,beta_lo,beta_hi);

    FORT_INTERP(phiPtr, ARLIM(phi_lo),ARLIM(phi_hi),
                tempPtr, ARLIM(temp_lo),ARLIM(temp_hi), 
                deltacPtr, ARLIM(deltac_lo),ARLIM(deltac_hi), 
                betaPtr, ARLIM(beta_lo),ARLIM(beta_hi), 
                lo,hi,loc,hic,isPeriodic);
}

void
ProjOutFlowBC_MG::gsrb (int nstep)
{
    FArrayBox dgphi(semiSurroundingNodes(domain,BL_SPACEDIM-1),1);

    DEF_BOX_LIMITS(domain,lo,hi);
    DEF_LIMITS(*phi,phiPtr,philo,phihi);
    DEF_LIMITS(*beta,  betaPtr, betalo,betahi);
    DEF_LIMITS(*rhs, rhsPtr, rhslo,rhshi);
    DEF_LIMITS(dgphi,dgphiPtr,dglo,dghi);

    FORT_HGRELAX(ARLIM(rhslo),ARLIM(rhshi),rhsPtr,
                 ARLIM(betalo),ARLIM(betahi),betaPtr,
                 ARLIM(philo),ARLIM(phihi),phiPtr,
                 ARLIM(dglo),ARLIM(dghi),dgphiPtr,
                 lo,hi,h,isPeriodic,&nstep);
}
