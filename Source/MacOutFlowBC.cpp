

#include <MacOutFlowBC.H>
#include <MACOUTFLOWBC_F.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

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

#if (BL_SPACEDIM == 3)
//
// Set defaults for these in the appropriate Initialize() routine!!!
//
Real MacOutFlowBC::tol;
Real MacOutFlowBC::abs_tol;

Real MacOutFlowBC_MG::cg_tol;
int  MacOutFlowBC_MG::verbose;
int  MacOutFlowBC_MG::maxIters;
int  MacOutFlowBC_MG::cg_maxiter;
Real MacOutFlowBC_MG::cg_abs_tol;
Real MacOutFlowBC_MG::cg_max_jump;
bool MacOutFlowBC_MG::useCGbottomSolver;
#endif

namespace
{
    bool outflow_initialized    = false;

#if (BL_SPACEDIM == 3)
    bool outflow_mg_initialized = false;
#endif
}

MacOutFlowBC::MacOutFlowBC ()
{
    Initialize();
}

void
MacOutFlowBC::Initialize ()
{
#if (BL_SPACEDIM == 3)
    if (outflow_initialized) return;
    //
    // Set defaults here!!!
    //
    MacOutFlowBC::tol     = 1.0e-10;
    MacOutFlowBC::abs_tol = 5.0e-10;

    ParmParse pp("macoutflow");

    pp.query("tol",     tol);
    pp.query("abs_tol", abs_tol);

    amrex::ExecOnFinalize(MacOutFlowBC::Finalize);

    outflow_initialized = true;
#endif
}

void
MacOutFlowBC::Finalize ()
{
    outflow_initialized = false;
}

void 
MacOutFlowBC::computeBC (FArrayBox         velMF[][2*BL_SPACEDIM],
                         FArrayBox         divuMF[2*BL_SPACEDIM],
                         FArrayBox         rhoMF[2*BL_SPACEDIM],
                         FArrayBox         phiMF[2*BL_SPACEDIM],
                         const Geometry&   geom, 
                         Orientation*      outFaces,
                         int               numOutFlowFaces,
                         const int*        lo_bc,
                         const int*        hi_bc,
                         Real              gravity)
{
    const Real small_udiff = 1.e-10;
    computeBC(velMF,divuMF,rhoMF,phiMF,geom,outFaces,numOutFlowFaces,lo_bc,hi_bc,small_udiff,gravity);
}

void 
MacOutFlowBC::computeBC (FArrayBox         velMF[][2*BL_SPACEDIM],
                         FArrayBox         divuMF[2*BL_SPACEDIM],
                         FArrayBox         rhoMF[2*BL_SPACEDIM],
                         FArrayBox         phiMF[2*BL_SPACEDIM],
                         const Geometry&   geom, 
                         Orientation*      outFaces,
                         int               numOutFlowFaces,
                         const int*        lo_bc,
                         const int*        hi_bc,
                         Real              small_udiff,
                         Real              gravity)
{
    BL_ASSERT(numOutFlowFaces <= 2*BL_SPACEDIM);

    int i,iface;

    int faces[2*BL_SPACEDIM];
    for (i = 0; i < numOutFlowFaces; i++)
        faces[i] = int(outFaces[i]);

    const Real* dx     = geom.CellSize();
    const Box&  domain = geom.Domain();

    int zeroIt[2*BL_SPACEDIM];
    for (i = 0; i < numOutFlowFaces; i++)
        zeroIt[i] = 0;

#if (BL_SPACEDIM == 2)
    Vector<Vector<Real> > redge(2*BL_SPACEDIM);
#endif

    FArrayBox ccExt[2*BL_SPACEDIM];
  
    int isPeriodic[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
        isPeriodic[dir] = geom.isPeriodic(dir);
  
    IntVect loFiltered, hiFiltered;
    int isPeriodicFiltered[2*BL_SPACEDIM][BL_SPACEDIM];
    Real dxFiltered[2*BL_SPACEDIM][BL_SPACEDIM];

    for (iface = 0; iface < numOutFlowFaces; iface++)
    {
        const int   outDir = outFaces[iface].coordDir();
        //
        // Filter out direction we don't care about.
        //
        int ncStripWidth = 1;
        Box origBox = amrex::adjCell(domain,outFaces[iface],ncStripWidth);
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
                dxFiltered[iface][cnt] = dx[dir];
                isPeriodicFiltered[iface][cnt] = isPeriodic[dir];
                cnt++;
            }
            else
            {
                loFiltered[BL_SPACEDIM-1] = lo[dir];
                hiFiltered[BL_SPACEDIM-1] = hi[dir];
                dxFiltered[iface][BL_SPACEDIM-1] = dx[dir];
                isPeriodicFiltered[iface][BL_SPACEDIM-1] = isPeriodic[dir];
            }
        }

        Box       faceBox(loFiltered,hiFiltered);
        //
        //  One for rho, one for divu.
        //
        ccExt[iface].resize(faceBox,2);
  
#if (BL_SPACEDIM == 2)
        //
        // Make edge-centered r (set = 1 if cartesian).
        //
        int perpDir = 1 - outDir;
        int r_len = domain.length(perpDir)+1;
        redge[iface].resize(r_len);
        //
        // Here we know the ordering of faces is XLO,YLO,XHI,YHI.
        //
        if (geom.IsRZ())
        {
            if (faces[iface] == 0)
            {
                for (i=0;i<r_len;i++)
                    redge[iface][i] = geom.ProbLo()[0];
            }
            else if (faces[iface] == 2)
            {
                for (i=0;i<r_len;i++)
                    redge[iface][i] = geom.ProbHi()[0];
            }
            else if (faces[iface] == 1 || faces[iface]== 3)
            {
                for (i=0;i<r_len;i++)
                    redge[iface][i] = geom.ProbLo()[0] +i*dx[0];
            }
        }
        else
        {
            for (i = 0; i < r_len; i++)
                redge[iface][i] = 1.;
        }
#endif
   
        DEF_BOX_LIMITS(origBox,origLo,origHi);

        const int* ccElo = ccExt[iface].loVect();
        const int* ccEhi = ccExt[iface].hiVect();
        const Real*  rhoEPtr = ccExt[iface].dataPtr(0);
        const Real* divuEPtr = ccExt[iface].dataPtr(1);

        DEF_LIMITS(divuMF[iface],divuPtr,divulo, divuhi);
        DEF_LIMITS( rhoMF[iface], rhoPtr, rholo, rhohi);
        DEF_LIMITS(velMF[0][iface],velXPtr,velXlo,velXhi);
        DEF_LIMITS(velMF[1][iface],velYPtr,velYlo,velYhi);
#if (BL_SPACEDIM == 3)
        DEF_LIMITS(velMF[2][iface],velZPtr,velZlo,velZhi);
#endif
        //
        // Extrapolate divu, and rho to the outflow edge in
        // the shifted coordinate system (where the last dimension is 1),
        // and replace (divu) by (divu - d/dperpdir (vel)).
        //
        extrap_mac(
            ARLIM(velXlo), ARLIM(velXhi), velXPtr,
            ARLIM(velYlo), ARLIM(velYhi), velYPtr,
#if (BL_SPACEDIM == 3)
            ARLIM(velZlo), ARLIM(velZhi), velZPtr,
#endif
            ARLIM(divulo),ARLIM(divuhi),divuPtr,
            ARLIM(rholo), ARLIM(rhohi), rhoPtr,
#if (BL_SPACEDIM == 2)
            &r_len, redge[iface].dataPtr(),
#endif
            ARLIM(ccElo),ARLIM(ccEhi),divuEPtr,
            ARLIM(ccElo),ARLIM(ccEhi),rhoEPtr,
            dx,
            origLo,origHi,&faces[iface],isPeriodicFiltered[iface],&zeroIt[iface],&small_udiff);
    }

    int connected = 0;
    //
    //  Test for whether multiple faces are touching.
    //
    if (numOutFlowFaces == 2)
    {
        if (outFaces[0].coordDir() != outFaces[1].coordDir())
            connected = 1;
    }
    else if (numOutFlowFaces > 2)
    {
        connected = 1;
    }
  
    if (connected == 0)
    {
	Vector<Real> x;

        for (iface = 0; iface < numOutFlowFaces; iface++)
        {
            zeroIt[iface] = 1; // HACK HACK

            if (zeroIt[iface])
            {
                phiMF[iface].setVal<RunOn::Host>(0);
            }
            else
            {
#if (BL_SPACEDIM == 2)
                int length = ccExt[iface].length()[0];

                const Real*  rhoEPtr = ccExt[iface].dataPtr(0);
                const Real* divuEPtr = ccExt[iface].dataPtr(1);
  
                Box faceBox(ccExt[iface].box());
                DEF_LIMITS(phiMF[iface], phiPtr,philo,phihi);

                x.resize(length);
                macphibc(x.dataPtr(), &length, divuEPtr, rhoEPtr,
                              redge[iface].dataPtr(), dxFiltered[iface],
                              isPeriodicFiltered[iface]);

                macphi_from_x(ARLIM(philo),ARLIM(phihi),phiPtr,
                                   &length,x.dataPtr());
#elif (BL_SPACEDIM == 3)
                Box faceBox(ccExt[iface].box());

                DEF_LIMITS(phiMF[iface],phiFabPtr,phiFab_lo,phiFab_hi);

                FArrayBox rhs, beta[BL_SPACEDIM-1];
      
                computeCoefficients(rhs,beta,
                                    ccExt[iface],
                                    faceBox,dxFiltered[iface],isPeriodicFiltered[iface]);
                //
                // Need phi to have ghost cells.
                //
                Box phiGhostBox = OutFlowBC::SemiGrow(faceBox,1,BL_SPACEDIM-1);
                FArrayBox phi(phiGhostBox,1), resid(rhs.box(),1);
                phi.setVal<RunOn::Host>(0);
            
                MacOutFlowBC_MG mac_mg(faceBox,&phi,&rhs,&resid,beta,
                                       dxFiltered[iface],isPeriodicFiltered[iface]);
      
                mac_mg.solve(tol,abs_tol,2,2,mac_mg.MaxIters(),mac_mg.Verbose());
      
                DEF_LIMITS(phi,phiPtr,phi_lo,phi_hi);
                //DEF_BOX_LIMITS(faceBox,lo,hi);
                //
                // Subtract the average phi.
                //
                //macsubtractavgphi(ARLIM(phi_lo),ARLIM(phi_hi),phiPtr,
                //             &r_len,redge,
                //             lo,hi,isPeriodicFiltered[iface]);
                //
                // Translate the solution back to the original coordinate system.
                //
                int face   = int(outFaces[iface]);
                mac_reshift_phi(ARLIM(phiFab_lo),ARLIM(phiFab_hi),phiFabPtr,
                                     ARLIM(phi_lo),ARLIM(phi_hi),phiPtr,&face);
#endif
            }
        }
    }

    if (connected == 1)
    {
        //
        // Define connected region.  In both 2-d and 3-d, if there are
        // multiple outflow faces and it's not just two across from
        // each other, then the multiple faces form a *single*
        // connected region.
        //
        int zeroAll = 1;
        for (i = 0; i < numOutFlowFaces; i++)
            if (zeroIt[i] == 0) zeroAll = 0;

        zeroAll = 1;  // HACK HACK

        if (zeroAll)
        {
            for ( i = 0; i < numOutFlowFaces; i++ )
                phiMF[i].setVal<RunOn::Host>(0);
        }
        else
        {
            //
            // Since we only use a constant dx in the Fortran,
            //  we'll assume for now we can choose either one.
            //
            BL_ASSERT(dx[0] == dx[1]);

            int lenx = domain.length(0);
            int leny = domain.length(1);

            int length = 0;

#if (BL_SPACEDIM == 2)
            //
            // Here we know the ordering of faces is XLO,YLO,XHI,YHI.
            //
            Real *ccEptr0,*ccEptr1,*ccEptr2,*ccEptr3;
            for (i=0; i < numOutFlowFaces; i++)
            {
                if (faces[i] == 0)
                {
                    ccEptr0 = ccExt[i].dataPtr();
                    length = length + leny;
                }
                else if (faces[i] == 1)
                {
                    ccEptr1 = ccExt[i].dataPtr();
                    length = length + lenx;
                }
                else if (faces[i] == 2)
                {
                    ccEptr2 = ccExt[i].dataPtr();
                    length = length + leny;
                }
                else if (faces[i] == 3)
                {
                    ccEptr3 = ccExt[i].dataPtr();
                    length = length + lenx;
                }
            }

            Real *r0,*r1,*r2,*r3;
            for (i=0; i < numOutFlowFaces; i++)
            {
                if (faces[i] == 0)
                {
                    r0 = redge[i].dataPtr();
                }
                else if (faces[i] == 1)
                {
                    r1 = redge[i].dataPtr();
                }
                else if (faces[i] == 2)
                {
                    r2 = redge[i].dataPtr();
                }
                else if (faces[i] == 3)
                {
                    r3 = redge[i].dataPtr();
                }
            }

#elif (BL_SPACEDIM == 3)
            //
            // Here we know the ordering of faces is XLO,YLO,ZLO,XHI,YHI,ZHI.
            //
            int lenz = domain.length(2);

            Real *ccEptr0,*ccEptr1,*ccEptr2,*ccEptr3,*ccEptr4,*ccEptr5;
            for (i=0; i < numOutFlowFaces; i++)
            {
                if (faces[i] == 0)
                {
                    ccEptr0 = ccExt[i].dataPtr();
                    length = length + leny*lenz;
                }
                else if (faces[i] == 1)
                {
                    ccEptr1 = ccExt[i].dataPtr();
                    length = length + lenx*lenz;
                }
                else if (faces[i] == 2)
                {
                    ccEptr2 = ccExt[i].dataPtr();
                    length = length + lenx*leny;
                }
                else if (faces[i] == 3)
                {
                    ccEptr3 = ccExt[i].dataPtr();
                    length = length + leny*lenz;
                }
                else if (faces[i] == 4)
                {
                    ccEptr4 = ccExt[i].dataPtr();
                    length = length + lenx*lenz;
                }
                else if (faces[i] == 5)
                {
                    ccEptr5 = ccExt[i].dataPtr();
                    length = length + lenx*leny;
                }
            }
#endif
            IntVect loconn, hiconn;

            loconn[0] = 0;
            hiconn[0] = length-1;
#if (BL_SPACEDIM == 3)
            loconn[1] = 0;
            hiconn[1] = lenz-1;
#endif
            loconn[BL_SPACEDIM-1] = 0;
            hiconn[BL_SPACEDIM-1] = 0;
            Box connected_region(loconn,hiconn);
            FArrayBox ccE_conn(connected_region,BL_SPACEDIM+1);
            ccE_conn.setVal<RunOn::Host>(1.e200);

#if (BL_SPACEDIM == 2)
            FArrayBox x(connected_region,1);

            int per = (numOutFlowFaces == 2*BL_SPACEDIM) ? 1 : 0;

	    Vector<Real> redge_conn(length+1);

            macfill_oned(&lenx,&leny,&length,faces,&numOutFlowFaces,
                              ccEptr0,ccEptr1,ccEptr2,ccEptr3,
                              r0,r1,r2,r3,
                              ccE_conn.dataPtr(),redge_conn.dataPtr());

            macphibc(x.dataPtr(), &length, 
                          ccE_conn.dataPtr(1), ccE_conn.dataPtr(0),
                          redge_conn.dataPtr(), dx, &per);

            Real *phiptr0, *phiptr1, *phiptr2, *phiptr3;
            for (int i=0; i < numOutFlowFaces; i++)
            {
                if (faces[i] == 0)
                {
                    phiptr0 = phiMF[i].dataPtr();
                }
                if (faces[i] == 1)
                {
                    phiptr1 = phiMF[i].dataPtr();
                }
                if (faces[i] == 2)
                {
                    phiptr2 = phiMF[i].dataPtr();
                }
                if (faces[i] == 3)
                {
                    phiptr3 = phiMF[i].dataPtr();
                }
            }

            macallphi_from_x(&lenx,&leny,&length,faces,&numOutFlowFaces,
                                  phiptr0, phiptr1, phiptr2, phiptr3,
                                  x.dataPtr());
#elif (BL_SPACEDIM == 3)

            int width = lenz;
            macfill_twod(&lenx,&leny,&length,&width,faces,&numOutFlowFaces,
                              ccEptr0,ccEptr1,ccEptr2,ccEptr3,ccEptr4,ccEptr5,
                              ccE_conn.dataPtr());
            ccE_conn.clear();

            Box faceBox(ccExt[iface].box());

            DEF_LIMITS(phiMF[iface],phiFabPtr,phiFab_lo,phiFab_hi);
      
            FArrayBox rhs, beta[BL_SPACEDIM-1];
      
            computeCoefficients(rhs,beta,
                                ccExt[iface],
                                faceBox,dxFiltered[iface],isPeriodicFiltered[iface]);
            //
            // Need phi to have ghost cells.
            //
            Box phiGhostBox = OutFlowBC::SemiGrow(faceBox,1,BL_SPACEDIM-1);
            FArrayBox phi(phiGhostBox,1);
            phi.setVal<RunOn::Host>(0);
      
            FArrayBox resid(rhs.box(),1);
      
            MacOutFlowBC_MG mac_mg(faceBox,&phi,&rhs,&resid,beta,
                                   dxFiltered[iface],isPeriodicFiltered[iface]);
      
            mac_mg.solve(tol,abs_tol,2,2,mac_mg.MaxIters(),mac_mg.Verbose());
      
            DEF_LIMITS(phi,phiPtr,phi_lo,phi_hi);
            //DEF_BOX_LIMITS(faceBox,lo,hi);
            //
            // Subtract the average phi.
            //
            //macsubtractavgphi(ARLIM(phi_lo),ARLIM(phi_hi),phiPtr,
            //                   &r_len,redge,
            //                   lo,hi,isPeriodicFiltered[iface]);
            //
            // Translate the solution back to the original coordinate system.
            //
            int face   = int(outFaces[iface]);
            mac_reshift_phi(ARLIM(phiFab_lo),ARLIM(phiFab_hi),phiFabPtr,
                                 ARLIM(phi_lo),ARLIM(phi_hi),phiPtr,&face);
#endif
        }
    }
}

#if (BL_SPACEDIM == 3)
void 
MacOutFlowBC::computeCoefficients (FArrayBox&   rhs,
                                   FArrayBox*   beta,
                                   FArrayBox&   ccExt,
                                   Box&         faceBox,
                                   Real*        dxFiltered,
				   int*         isPeriodicFiltered)
{
    rhs.resize(faceBox,1);
    beta[0].resize(amrex::surroundingNodes(faceBox,0),1);
    beta[1].resize(amrex::surroundingNodes(faceBox,1),1);
   
    DEF_BOX_LIMITS(faceBox,faceLo,faceHi);
    const int* ccElo = ccExt.loVect();
    const int* ccEhi = ccExt.hiVect();
    const Real*  rhoEPtr = ccExt.dataPtr(0);
    const Real* divuEPtr = ccExt.dataPtr(1);
    DEF_LIMITS(rhs, rhsPtr, rhslo,rhshi);
    DEF_LIMITS(beta[0],beta0Ptr, beta0lo, beta0hi);
    DEF_LIMITS(beta[1],beta1Ptr, beta1lo, beta1hi);


    compute_maccoeff(ARLIM(rhslo),ARLIM(rhshi),rhsPtr,
                          ARLIM(beta0lo),ARLIM(beta0hi),beta0Ptr,
                          ARLIM(beta1lo),ARLIM(beta1hi),beta1Ptr,
                          ARLIM(ccElo),ARLIM(ccEhi),divuEPtr,
                          ARLIM(ccElo),ARLIM(ccEhi),rhoEPtr,
                          faceLo,faceHi,
                          dxFiltered,isPeriodicFiltered);
}

void
MacOutFlowBC_MG::Initialize ()
{
    if (outflow_mg_initialized) return;
    //
    // Set defaults here!!!
    //
    MacOutFlowBC_MG::cg_tol            = 1.0e-2;
    MacOutFlowBC_MG::verbose           = 0;
    MacOutFlowBC_MG::maxIters          = 40;
    MacOutFlowBC_MG::cg_maxiter        = 40;
    MacOutFlowBC_MG::cg_abs_tol        = 5.0e-12;
    MacOutFlowBC_MG::cg_max_jump       = 10.0;
    MacOutFlowBC_MG::useCGbottomSolver = true;

    ParmParse pp("mac_mg");

    pp.query("v",                 verbose);
    pp.query("cg_tol",            cg_tol);
    pp.query("maxIters",          maxIters);
    pp.query("cg_abs_tol",        cg_abs_tol);
    pp.query("cg_maxiter",        cg_maxiter);
    pp.query("cg_max_jump",       cg_max_jump);
    pp.query("useCGbottomSolver", useCGbottomSolver);

    amrex::ExecOnFinalize(MacOutFlowBC_MG::Finalize);

    outflow_mg_initialized = true;
}

void
MacOutFlowBC_MG::Finalize ()
{
    outflow_mg_initialized = false;
}

MacOutFlowBC_MG::MacOutFlowBC_MG (Box&       Domain,
                                  FArrayBox* Phi,
                                  FArrayBox* Rhs,
                                  FArrayBox* Resid,
                                  FArrayBox* Beta,
                                  Real*      H,
                                  int*       IsPeriodic)
    :
    OutFlowBC_MG(Domain,Phi,Rhs,Resid,Beta,H,IsPeriodic,false)
{
    Initialize();

    IntVect len = domain.size();

    int min_length = 4;
    bool test_side[BL_SPACEDIM-1];
    for (int dir = 0; dir < BL_SPACEDIM-1; dir++)
        test_side[dir] = (len[dir]&1) != 0 || len[dir] < min_length;

#if BL_SPACEDIM==3
    bool doit = test_side[0] || test_side[1];
#elif BL_SPACEDIM==2
    bool doit = test_side[0];
#else
    bool doit = true;
#endif

    if (doit)
    {
        if (useCGbottomSolver)
            cgwork.resize(domain,4);
    }
    else
    {
        Real newh[BL_SPACEDIM];
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            newh[dir] = 2*h[dir];

        Box newdomain = OutFlowBC::SemiCoarsen(domain,2,BL_SPACEDIM-1);
        Box grownBox  = OutFlowBC::SemiGrow(newdomain,1,BL_SPACEDIM-1);
      
        FArrayBox* newphi    = new FArrayBox(grownBox,1);
        FArrayBox* newresid  = new FArrayBox(newdomain,1);
        FArrayBox* newrhs    = new FArrayBox(newdomain,1);
        FArrayBox* newbeta = new FArrayBox[BL_SPACEDIM-1];
        newbeta[0].resize(amrex::surroundingNodes(newdomain,0),1);
        newbeta[1].resize(amrex::surroundingNodes(newdomain,1),1);

        // FIXME MSD: Combine these
        newphi->setVal<RunOn::Host>(0);
        newresid->setVal<RunOn::Host>(0);
        newbeta[0].setVal<RunOn::Host>(0);
        newbeta[1].setVal<RunOn::Host>(0);

        DEF_BOX_LIMITS(domain,dom_lo,dom_hi);
        DEF_BOX_LIMITS(newdomain,new_lo,new_hi);
        DEF_LIMITS(beta[0],beta0Ptr,beta0_lo,beta0_hi);
        DEF_LIMITS(beta[1],beta1Ptr,beta1_lo,beta1_hi);
        DEF_LIMITS(newbeta[0],newbeta0Ptr,newbeta0_lo,newbeta0_hi);
        DEF_LIMITS(newbeta[1],newbeta1Ptr,newbeta1_lo,newbeta1_hi);
      
        coarsigma(beta0Ptr,ARLIM(beta0_lo),ARLIM(beta0_hi),
#if (BL_SPACEDIM == 3)
                       beta1Ptr,ARLIM(beta1_lo),ARLIM(beta1_hi),
#endif
                       newbeta0Ptr,ARLIM(newbeta0_lo),ARLIM(newbeta0_hi),
#if (BL_SPACEDIM == 3)
                       newbeta1Ptr,ARLIM(newbeta1_lo),ARLIM(newbeta1_hi),
#endif
                       dom_lo,dom_hi,new_lo,new_hi);
      
        next = new MacOutFlowBC_MG(newdomain,newphi,newrhs,
                                   newresid,newbeta,newh,isPeriodic);
    }
}

int
MacOutFlowBC_MG::Verbose ()
{
    return verbose;
}

int
MacOutFlowBC_MG::MaxIters ()
{
    return maxIters;
}

MacOutFlowBC_MG::~MacOutFlowBC_MG () {}

Real
MacOutFlowBC_MG::residual ()
{
    Real rnorm;

    DEF_BOX_LIMITS(domain,lo,hi);
    DEF_LIMITS(*rhs,rhsPtr,rhslo,rhshi);
    DEF_LIMITS(*resid,residPtr,residlo,residhi);
    DEF_LIMITS(*phi,phiPtr,philo,phihi);
    DEF_LIMITS(beta[0],beta0Ptr,beta0lo,beta0hi);
    DEF_LIMITS(beta[1],beta1Ptr,beta1lo,beta1hi);

    macresid(ARLIM(rhslo),ARLIM(rhshi),rhsPtr,
                  ARLIM(beta0lo),ARLIM(beta0hi),beta0Ptr,
                  ARLIM(beta1lo),ARLIM(beta1hi),beta1Ptr,
                  ARLIM(philo), ARLIM(phihi),phiPtr,
                  ARLIM(residlo),ARLIM(residhi), residPtr,
                  lo,hi,h,isPeriodic,&rnorm);
  
    return rnorm;
}

void 
MacOutFlowBC_MG::step (int nGSRB)
{
    if (cgwork.isAllocated())
    {
        Real resnorm = 0.0;

        FArrayBox dest0(phi->box(),1);

        DEF_BOX_LIMITS(domain,lo,hi);
        DEF_LIMITS(*phi,phiPtr,phi_lo,phi_hi);
        DEF_LIMITS(*resid,residPtr,resid_lo,resid_hi);
        DEF_LIMITS(dest0,dest0Ptr,dest0_lo,dest0_hi);
        DEF_LIMITS(*rhs,rhsPtr,rhs_lo,rhs_hi);
        DEF_LIMITS(beta[0], beta0Ptr, beta0_lo,beta0_hi); 
        DEF_LIMITS(beta[1], beta1Ptr, beta1_lo,beta1_hi);
        DEF_BOX_LIMITS(cgwork,cg_lo,cg_hi);
    
        solvemac(phiPtr, ARLIM(phi_lo),ARLIM(phi_hi),
                      dest0Ptr,ARLIM(dest0_lo),ARLIM(dest0_hi),
                      rhsPtr, ARLIM(rhs_lo),ARLIM(rhs_hi),
                      beta0Ptr, ARLIM(beta0_lo),ARLIM(beta0_hi),
                      beta1Ptr,ARLIM(beta1_lo),ARLIM(beta1_hi),
                      cgwork.dataPtr(0), ARLIM(cg_lo),ARLIM(cg_hi),
                      cgwork.dataPtr(1), ARLIM(cg_lo),ARLIM(cg_hi),
                      cgwork.dataPtr(2), ARLIM(cg_lo),ARLIM(cg_hi),
                      cgwork.dataPtr(3), ARLIM(cg_lo),ARLIM(cg_hi),
                      residPtr, ARLIM(resid_lo), ARLIM(resid_hi),
                      lo,hi,h,isPeriodic,&cg_maxiter,&cg_tol,
                      &cg_abs_tol,&cg_max_jump,&resnorm);
    }
    else
    {
        gsrb(nGSRB);
    }
}

void 
MacOutFlowBC_MG::Restrict ()
{
    DEF_BOX_LIMITS(domain,lo,hi);
    DEF_BOX_LIMITS(next->theDomain(),loc,hic);
    DEF_LIMITS(*resid,residPtr,resid_lo,resid_hi);
    DEF_LIMITS(*(next->theRhs()),rescPtr,resc_lo,resc_hi);

    outflowbc_restrict(residPtr, ARLIM(resid_lo),ARLIM(resid_hi), 
                  rescPtr, ARLIM(resc_lo),ARLIM(resc_hi), 
                  lo,hi,loc,hic);
}

void 
MacOutFlowBC_MG::interpolate ()
{
    DEF_BOX_LIMITS(domain,lo,hi);
    DEF_BOX_LIMITS(next->theDomain(),loc,hic);
    DEF_LIMITS(*phi,phiPtr,phi_lo,phi_hi);
    DEF_LIMITS(*(next->thePhi()),deltacPtr,deltac_lo,deltac_hi);

    fort_interpolate(phiPtr, ARLIM(phi_lo),ARLIM(phi_hi), 
                     deltacPtr,ARLIM(deltac_lo),ARLIM(deltac_hi), 
                     lo,hi,loc,hic);
}

void 
MacOutFlowBC_MG::gsrb (int nstep)
{
    DEF_BOX_LIMITS(domain,lo,hi);
    DEF_LIMITS(*rhs, rhsPtr, rhslo,rhshi);
    DEF_LIMITS(beta[0], beta0Ptr, beta0lo, beta0hi);
    DEF_LIMITS(beta[1], beta1Ptr, beta1lo, beta1hi);
    DEF_LIMITS(*phi,phiPtr,philo,phihi);

    macrelax(ARLIM(rhslo),ARLIM(rhshi),rhsPtr,
                  ARLIM(beta0lo),ARLIM(beta0hi),beta0Ptr,
                  ARLIM(beta1lo),ARLIM(beta1hi),beta1Ptr,
                  ARLIM(philo),ARLIM(phihi),phiPtr,
                  lo,hi,h,isPeriodic,&nstep);
}
#endif
