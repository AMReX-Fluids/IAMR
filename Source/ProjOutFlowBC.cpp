

#include <iostream>

#include <ProjOutFlowBC.H>
#include <PROJOUTFLOWBC_F.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

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

namespace
{
    bool outflowbc_initialized    = false;

#if (BL_SPACEDIM == 3)
    bool outflowbc_mg_initialized = false;
#endif
}

#if (BL_SPACEDIM == 3)
//
// Set defaults in appropriate Initialize() routine!!!
//
Real ProjOutFlowBC::tol;
Real ProjOutFlowBC::abs_tol;

Real ProjOutFlowBC_MG::cg_tol;
int  ProjOutFlowBC_MG::verbose;
int  ProjOutFlowBC_MG::maxIters;
Real ProjOutFlowBC_MG::cg_abs_tol;
int  ProjOutFlowBC_MG::cg_maxiter;
Real ProjOutFlowBC_MG::cg_max_jump;
bool ProjOutFlowBC_MG::useCGbottomSolver;

static
Box
semiSurroundingNodes (const Box& baseBox,
                      int        direction)
{
    Box sBox = amrex::surroundingNodes(baseBox);
    sBox.growHi(direction,-1);
    return sBox;
}
#endif

void
ProjOutFlowBC::Initialize ()
{
#if (BL_SPACEDIM == 3)
    if (outflowbc_initialized) return;
    //
    // Set defaults here !!!
    //
    ProjOutFlowBC::tol     = 1.0e-10; 
    ProjOutFlowBC::abs_tol = 5.0e-10;

    ParmParse pp("projoutflow");

    pp.query("tol",tol);
    pp.query("abs_tol",abs_tol);

    amrex::ExecOnFinalize(ProjOutFlowBC::Finalize);

    outflowbc_initialized = true;
#endif
}

ProjOutFlowBC::ProjOutFlowBC ()
{
    Initialize();
}

void
ProjOutFlowBC::Finalize ()
{
    outflowbc_initialized = false;
}

void 
ProjOutFlowBC::computeBC (FArrayBox       velMF[][2*BL_SPACEDIM],
                          FArrayBox        divuMF[2*BL_SPACEDIM],
                          FArrayBox        rhoMF[2*BL_SPACEDIM],
                          FArrayBox        phiMF[2*BL_SPACEDIM],
                          const Geometry&   geom, 
                          Orientation*      outFaces,
                          int               numOutFlowFaces,
                          const int*        lo_bc,
                          const int*        hi_bc,
                          Real              small_udiff,
                          Real              gravity)
{
  computeBC (velMF,divuMF,rhoMF,phiMF,geom,outFaces,numOutFlowFaces,lo_bc,hi_bc,gravity);
}

void 
ProjOutFlowBC::computeBC (FArrayBox       velMF[][2*BL_SPACEDIM],
                          FArrayBox        divuMF[2*BL_SPACEDIM],
                          FArrayBox        rhoMF[2*BL_SPACEDIM],
                          FArrayBox        phiMF[2*BL_SPACEDIM],
                          const Geometry&   geom, 
                          Orientation*      outFaces,
                          int               numOutFlowFaces,
                          const int*        lo_bc,
                          const int*        hi_bc,
                          Real              gravity)
{
    BL_ASSERT(numOutFlowFaces <= 2*BL_SPACEDIM);

    int faces[2*BL_SPACEDIM];
    for (int i = 0; i < numOutFlowFaces; i++)
        faces[i] = int(outFaces[i]);

    const Real* dx    = geom.CellSize();
    const Box& domain = geom.Domain();

    int lenx = domain.length(0);
    int leny = domain.length(1);

    int zeroIt[2*BL_SPACEDIM];
    for (int i = 0; i < numOutFlowFaces; i++)
        zeroIt[i] = 0;

#if (BL_SPACEDIM == 2)
    Vector<Vector<Real> > rcen(2*BL_SPACEDIM);
    Vector<Vector<Real> > redge(2*BL_SPACEDIM);
#endif

    FArrayBox ccExt[2*BL_SPACEDIM];

    int isPeriodic[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
        isPeriodic[dir] = geom.isPeriodic(dir);

    IntVect loFiltered, hiFiltered;
#if (BL_SPACEDIM == 2)
    int isPeriodicFiltered[2*BL_SPACEDIM][BL_SPACEDIM];
#endif
#if (BL_SPACEDIM == 3)
    Real dxFiltered[2*BL_SPACEDIM][BL_SPACEDIM];
#endif

    for (int iface = 0; iface < numOutFlowFaces; iface++)
    {
        int outDir        = outFaces[iface].coordDir();
        //
        // Filter out the direction we don't care about.
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
#if (BL_SPACEDIM == 3)
                dxFiltered[iface][cnt] = dx[dir];
#endif
#if (BL_SPACEDIM == 2)
                isPeriodicFiltered[iface][cnt] = isPeriodic[dir];
#endif
                cnt++;
            }
            else
            {
                loFiltered[BL_SPACEDIM-1] = lo[dir];
                hiFiltered[BL_SPACEDIM-1] = hi[dir];
#if (BL_SPACEDIM == 3)
                dxFiltered[iface][BL_SPACEDIM-1] = dx[dir];
#endif
#if (BL_SPACEDIM == 2)
                isPeriodicFiltered[iface][BL_SPACEDIM-1] = isPeriodic[dir];
#endif
            }
        }

        Box       faceBox(loFiltered,hiFiltered);
        for (int dir = 0; dir < BL_SPACEDIM-1; dir++)
            faceBox.grow(dir,1);
        //
        //  One for rho, one for divu, (BL_SPACEDIM-1) for velocity.
        //
        ccExt[iface].resize(faceBox,BL_SPACEDIM+1);

#if (BL_SPACEDIM == 2)
        //
        // Make edge-centered and cc r (set = 1 if cartesian)
        //
        int perpDir = 1 - outDir;
        int r_len = domain.length(perpDir)+1;
        rcen[iface].resize(r_len-1);
        redge[iface].resize(r_len);
        //
        // Here we know the ordering of faces is XLO,YLO,XHI,YHI.
        //
        if (geom.IsRZ())
        {
            if (faces[iface] == 0)
            {
                for (int i=0;i<r_len  ;i++)
                    redge[iface][i] = geom.ProbLo()[0];
                for (int i=0;i<r_len-1;i++)
                    rcen[iface][i] = geom.ProbLo()[0];
            }
            else if (faces[iface] == 2)
            {
                for (int i=0;i<r_len  ;i++)
                    redge[iface][i] = geom.ProbHi()[0];
                for (int i=0;i<r_len-1;i++)
                    rcen[iface][i] = geom.ProbHi()[0];
            }
            else if (faces[iface] == 1 || faces[iface]== 3)
            {
                for (int i=0;i<r_len  ;i++)
                    redge[iface][i] = geom.ProbLo()[0] + i     *dx[0];
                for (int i=0;i<r_len-1;i++)
                    rcen[iface][i] = geom.ProbLo()[0] +(i+0.5)*dx[0];
            }
        }
        else
        {
            for (int i = 0; i < r_len  ; i++)
                redge[iface][i] = 1.;
            for (int i = 0; i < r_len-1; i++)
                rcen[iface][i] = 1.;
        }
#endif

        DEF_BOX_LIMITS(origBox,origLo,origHi);

        const int* ccElo = ccExt[iface].loVect();
        const int* ccEhi = ccExt[iface].hiVect();
        const Real*  rhoEPtr = ccExt[iface].dataPtr(0);
        const Real* divuEPtr = ccExt[iface].dataPtr(1);
        const Real*    uEPtr = ccExt[iface].dataPtr(2);

        DEF_LIMITS(divuMF[iface], divuPtr, divulo, divuhi);
        DEF_LIMITS( rhoMF[iface],  rhoPtr,  rholo,  rhohi);
        DEF_LIMITS( velMF[0][iface],  velPtr,  vello,  velhi);
        //
        // Extrapolate the velocities, divu, and rho to the outflow edge in
        // the shifted coordinate system (where the last dimension is 1).
        //
        extrap_proj(ARLIM(vello),  ARLIM(velhi), velPtr,
                         ARLIM(divulo), ARLIM(divuhi), divuPtr,
                         ARLIM(rholo),  ARLIM(rhohi),rhoPtr,
#if (BL_SPACEDIM == 2)
                         &r_len,redge[iface].dataPtr(),
#endif
                         ARLIM(ccElo),ARLIM(ccEhi),uEPtr,
                         ARLIM(ccElo),ARLIM(ccEhi),divuEPtr,
                         ARLIM(ccElo),ARLIM(ccEhi),rhoEPtr,
                         origLo,origHi,&faces[iface],&zeroIt[iface]);
    }
    //
    //  Test for whether multiple faces are touching.
    //  therefore not touching.
    //
    int numRegions = 1;
    if ( (numOutFlowFaces == 2) &&
         (outFaces[0].coordDir() == outFaces[1].coordDir()) )
        numRegions = 2;
    //
    // Since we only use a constant dx in the Fortran,
    //  we'll assume for now we can choose either one.
    //
    if (numRegions == 1 && numOutFlowFaces > 1)
        BL_ASSERT(dx[0] == dx[1]);
    //
    //   Note numRegions = 1 or 2, those are the only possibilities.
    //
    for (int ireg = 0; ireg < numRegions; ireg++) 
    {
        //
        // Define connected region.  In both 2-d and 3-d, if there are
        // multiple outflow faces and it's not just two across from
        // each other, then the multiple faces form a *single* 
        // connected region.
        //
        int zeroAll = zeroIt[ireg];
        if (numRegions == 1)
            for (int i = 0; i < numOutFlowFaces; i++)
                if (zeroIt[i] == 0) zeroAll = 0;

        zeroAll = 1; // HACK HACK

        if (zeroAll)
        {
            for (int i=0; i < numOutFlowFaces; i++)
                phiMF[i].setVal<RunOn::Host>(0);
        }
        else
        {
            int faces2[2*BL_SPACEDIM];
#if (BL_SPACEDIM == 2)
            int numOutFlowFacesInRegion;
#endif
            if (numRegions == 1)
            {
                for (int i=0; i < numOutFlowFaces; i++) 
                    faces2[i] = int(outFaces[i]);
#if (BL_SPACEDIM == 2)
                numOutFlowFacesInRegion = numOutFlowFaces;
#endif
            }
            else if (numRegions == 2)
            {
                faces2[0] = int(outFaces[ireg]);
#if (BL_SPACEDIM == 2)
                numOutFlowFacesInRegion = 1;
#endif
            }

#if (BL_SPACEDIM == 2)
            //
            // Here we know the ordering of faces2 is XLO,XHI,YLO,YHI.
            //
            int length = 0;
            Real *ccEptr0,*ccEptr1,*ccEptr2,*ccEptr3;
            Real *r0,*r1,*r2,*r3;
            for (int i=0; i < numOutFlowFacesInRegion; i++) 
            {
                if (faces2[i] == 0)
                {
                    ccEptr0 = ccExt[i].dataPtr();
                    r0 = rcen[i].dataPtr();
                    length = length + leny;
                }
                else if (faces2[i] == 1)
                {
                    ccEptr1 = ccExt[i].dataPtr();
                    r1 = rcen[i].dataPtr();
                    length = length + lenx;
                }
                else if (faces2[i] == 2)
                {
                    ccEptr2 = ccExt[i].dataPtr();
                    r2 = rcen[i].dataPtr();
                    length = length + leny;
                }
                else if (faces2[i] == 3) {
                    ccEptr3 = ccExt[i].dataPtr();
                    r3 = rcen[i].dataPtr();
                    length = length + lenx;
                }
            }

            IntVect loconn;
            IntVect hiconn;

            loconn[0] = 0;
            hiconn[0] = length-1;
            loconn[BL_SPACEDIM-1] = 0;
            hiconn[BL_SPACEDIM-1] = 0;
            Box connected_region(loconn,hiconn);
            FArrayBox ccE_conn(connected_region,1);
  
            hiconn[0] = length;
            Box nodal_connected_region(loconn,hiconn);
            FArrayBox x(nodal_connected_region,1);
            FArrayBox s(nodal_connected_region,1);
            s.setVal<RunOn::Host>(0.);

            ccE_conn.setVal<RunOn::Host>(1.e200);

            int per = 0;
            if ( (numOutFlowFaces == 1) || 
                 (numRegions == 2) ) per = isPeriodicFiltered[ireg][0];

            fill_oned(&lenx,&leny,&length,faces2,&numOutFlowFacesInRegion,
                           ccEptr0, ccEptr1, ccEptr2, ccEptr3,
                           r0,r1,r2,r3,
                           ccE_conn.dataPtr(),s.dataPtr(),&per,
                           &(dx[0]),&(dx[1]));

            if (numOutFlowFaces == 2*BL_SPACEDIM) per = 1;

            hgphibc(dx,
                         ccE_conn.dataPtr(0),
                         s.dataPtr(),
                         x.dataPtr(),
                         &length,&per);

            Real *phiptr0,*phiptr1,*phiptr2,*phiptr3;

            for (int i=0; i < numOutFlowFacesInRegion; i++) 
            {
                if (faces2[i] == 0) {
                    phiptr0 = phiMF[i].dataPtr();
                }
                if (faces2[i] == 1) {
                    phiptr1 = phiMF[i].dataPtr();
                }
                if (faces2[i] == 2) {
                    phiptr2 = phiMF[i].dataPtr();
                }
                if (faces2[i] == 3) {
                    phiptr3 = phiMF[i].dataPtr();
                }
            }

            allphi_from_x(&lenx,&leny,&length,faces2,&numOutFlowFaces,
                               phiptr0, phiptr1, phiptr2, phiptr3,
                               x.dataPtr());
#else

#ifdef AMREX_DEBUG
            //
            // Assert that, if faces are connected, one of the coordinate
            // directions has no outflow faces.
            //
            int outx = 0, outy = 0, outz = 0;
            for (int iface = 0; iface < numOutFlowFaces; iface++)
            {
                int outDir = outFaces[iface].coordDir();
                if (outDir == 0) outx = 1;
                if (outDir == 1) outy = 1;
                if (outDir == 2) outz = 1;
            }
            BL_ASSERT ((outx + outy + outz) > 0);
            BL_ASSERT ((outx + outy + outz) < 3);
            // FOR NOW: ASSERT THAT NO OUTFLOW FACES IN Z-DIR!
            BL_ASSERT (outz == 0);
#endif
            BL_ASSERT(dx[1] == dx[2]);

            // Here we know the ordering of faces2 is XLO,YLO,ZLO,XHI,YHI,ZHI.

            int lenz = domain.length(2);

            int length = 0;
            int  width = lenz;
            Real *ccEptr0,*ccEptr1,*ccEptr2,*ccEptr3,*ccEptr4,*ccEptr5;
            for (int i=0; i < numOutFlowFaces; i++) 
            {
                if (faces2[i] == 0) {
                    ccEptr0 = ccExt[i].dataPtr();
                    length = length + leny*lenz;
                } else if (faces2[i] == 1) {
                    ccEptr1 = ccExt[i].dataPtr();
                    length = length + lenx*lenz;
                } else if (faces2[i] == 2) {
                    ccEptr2 = ccExt[i].dataPtr();
                    length = length + lenx*leny;
                } else if (faces2[i] == 3) {
                    ccEptr3 = ccExt[i].dataPtr();
                    length = length + leny*lenz;
                } else if (faces2[i] == 4) {
                    ccEptr4 = ccExt[i].dataPtr();
                    length = length + lenx*lenz;
                } else if (faces2[i] == 5) {
                    ccEptr5 = ccExt[i].dataPtr();
                    length = length + lenx*leny;
                } else {
   		    amrex::Print() << "OOPS - DIDNT PROGRAM FOR Z-OUTFLOW FACES! "
				   << i << " " << faces2[i] << std::endl;
		    exit(0);
                }
            }

            IntVect loconn;
            IntVect hiconn;

            loconn[0] = 0;
            hiconn[0] = length-1;
            loconn[1] = 0;
            hiconn[1] = width-1;
            loconn[BL_SPACEDIM-1] = 0;
            hiconn[BL_SPACEDIM-1] = 0;
            Box connected_region(loconn,hiconn);
            FArrayBox ccE_conn(connected_region,BL_SPACEDIM+1);
 
            hiconn[0] = length;
            hiconn[1] = width;
            Box nodal_connected_region(loconn,hiconn);
            FArrayBox phiFiltered(nodal_connected_region,1);
            phiFiltered.setVal<RunOn::Host>(0.);

            fill_twod(&lenx,&leny,&lenz,&length,&width,
                           faces2,&numOutFlowFaces,
                           ccEptr0, ccEptr1, ccEptr2, ccEptr3, ccEptr4, ccEptr5,
                           ccE_conn.dataPtr());
      
            FArrayBox rhs_temp, beta;
  
            int* per = new int[2];
            per[0] = (numOutFlowFaces == 2*BL_SPACEDIM) ? 1 : 0;
            per[1] = isPeriodic[BL_SPACEDIM-1];
        
            computeCoefficients(rhs_temp,beta,ccE_conn,connected_region,dxFiltered[0],per);
            //
            // Need phi to have ghost cells.
            //
            Box phiGhostBox = OutFlowBC::SemiGrow(phiFiltered.box(),1,BL_SPACEDIM-1);
            FArrayBox phi(phiGhostBox,1);
            phi.setVal<RunOn::Host>(0);
            phi.copy<RunOn::Host>(phiFiltered);
      
            Box grownRhs = OutFlowBC::SemiGrow(rhs_temp.box(),1,BL_SPACEDIM-1);
            FArrayBox rhs(grownRhs,1);
            rhs.setVal<RunOn::Host>(0);
            rhs.copy<RunOn::Host>(rhs_temp);
            FArrayBox resid(rhs.box(),1);
            ProjOutFlowBC_MG proj_mg(connected_region,&phi,&rhs,&resid,&beta,
                                     dxFiltered[0],per);

            proj_mg.solve(tol,abs_tol,2,2,proj_mg.MaxIters(),proj_mg.Verbose());
      
            DEF_BOX_LIMITS(phi,phi_lo,phi_hi);

            Real *phiptr0,*phiptr1,*phiptr2,*phiptr3,*phiptr4,*phiptr5;

            for (int i=0; i < numOutFlowFaces; i++) 
            {
                if (faces2[i] == 0) {
                    phiptr0 = phiMF[i].dataPtr();
                } else 
                    if (faces2[i] == 1) {
                        phiptr1 = phiMF[i].dataPtr();
                    } else 
                        if (faces2[i] == 2) {
                            phiptr2 = phiMF[i].dataPtr();
                        } else 
                            if (faces2[i] == 3) {
                                phiptr3 = phiMF[i].dataPtr();
                            } else 
                                if (faces2[i] == 4) {
                                    phiptr4 = phiMF[i].dataPtr();
                                } else 
                                    if (faces2[i] == 5) {
                                        phiptr5 = phiMF[i].dataPtr();
                                    }
            }
            allphi_from_x(&lenx,&leny,&lenz,&length,&width,faces2,&numOutFlowFaces,
                               phiptr0, phiptr1, phiptr2, phiptr3, phiptr4, phiptr5,
                               phi.dataPtr(),ARLIM(phi_lo),ARLIM(phi_hi));
#endif
        }
    }

    if (std::fabs(gravity) > 0.)
    {
        const int* domlo  = domain.loVect();
        const int* domhi  = domain.hiVect();
        for (int iface = 0; iface < numOutFlowFaces; iface++) 
        {
            int face          = int(outFaces[iface]);
            int outDir        = outFaces[iface].coordDir();

            DEF_LIMITS(phiMF[iface], phiPtr,philo,phihi);
            DEF_LIMITS(rhoMF[iface], rhoPtr,rholo,rhohi);
            if (outDir != (BL_SPACEDIM-1))
                rhogbc(rhoPtr,ARLIM(rholo),ARLIM(rhohi),
                            phiPtr,ARLIM(philo),ARLIM(phihi),
                            &face,&gravity,dx,domlo,domhi,
                            lo_bc,hi_bc);
  
        }
    }
}

void 
ProjOutFlowBC::computeRhoG (FArrayBox*         rhoMF,
                            FArrayBox*         phiMF,
                            const Geometry&    geom, 
                            Orientation*       outFaces,
                            int                numOutFlowFaces,
                            Real               gravity,
                            const int*         lo_bc,
                            const int*         hi_bc)

{
    const Real* dx    = geom.CellSize();
    const Box& domain = geom.Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    if (std::fabs(gravity) > 0.)
    {
      for (int iface = 0; iface < numOutFlowFaces; iface++) {

	int face          = int(outFaces[iface]);
	int outDir        = outFaces[iface].coordDir();
	
	DEF_LIMITS(phiMF[iface], phiPtr,philo,phihi);
	DEF_LIMITS(rhoMF[iface], rhoPtr,rholo,rhohi);
	
	if (outDir != (BL_SPACEDIM-1))
	  rhogbc(rhoPtr,ARLIM(rholo),ARLIM(rhohi),
		 phiPtr,ARLIM(philo),ARLIM(phihi),
		 &face,&gravity,dx,domlo,domhi,
		 lo_bc,hi_bc);
      }
    }
}

#if (BL_SPACEDIM == 3)
void 
ProjOutFlowBC::computeCoefficients (FArrayBox&   rhs,
                                    FArrayBox&   beta,
                                    FArrayBox&   ccExt,
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
    const int* ccElo = ccExt.loVect();
    const int* ccEhi = ccExt.hiVect();
    Real*  rhoEPtr = ccExt.dataPtr(0);
    Real* divuEPtr = ccExt.dataPtr(1);
    Real*    uEPtr = ccExt.dataPtr(2);

    compute_coeff(ARLIM(rhslo),ARLIM(rhshi),rhsPtr,
                       ARLIM(betalo),ARLIM(betahi),betaPtr,
                       ARLIM(ccElo),ARLIM(ccEhi),uEPtr,
                       ARLIM(ccElo),ARLIM(ccEhi),divuEPtr,
                       ARLIM(ccElo),ARLIM(ccEhi),rhoEPtr,
                       faceLo,faceHi,
                       dxFiltered,isPeriodicFiltered);
}

void
ProjOutFlowBC_MG::Initialize ()
{
    if (outflowbc_mg_initialized) return;
    //
    // Set defaults here!!!
    //
    ProjOutFlowBC_MG::cg_tol            = 1.0e-2;
    ProjOutFlowBC_MG::verbose           = 0;
    ProjOutFlowBC_MG::maxIters          = 40;
    ProjOutFlowBC_MG::cg_abs_tol        = 5.0e-12;
    ProjOutFlowBC_MG::cg_maxiter        = 40;
    ProjOutFlowBC_MG::cg_max_jump       = 10.0;
    ProjOutFlowBC_MG::useCGbottomSolver = true;

    ParmParse pp("proj_mg");

    pp.query("v",                 verbose);
    pp.query("cg_tol",            cg_tol);
    pp.query("maxIters",          maxIters);
    pp.query("cg_maxiter",        cg_maxiter);
    pp.query("cg_abs_tol",        cg_abs_tol);
    pp.query("cg_max_jump",       cg_max_jump);
    pp.query("useCGbottomSolver", useCGbottomSolver);

    amrex::ExecOnFinalize(ProjOutFlowBC_MG::Finalize);

    outflowbc_mg_initialized = true;
}

void
ProjOutFlowBC_MG::Finalize ()
{
    outflowbc_mg_initialized = false;
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
        {
            Box temp1Box = OutFlowBC::SemiGrow(domain,1,BL_SPACEDIM-1);
            Box temp2Box = semiSurroundingNodes(temp1Box,BL_SPACEDIM-1);
            cgwork.resize(temp2Box,4);
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
        // FIXME MSD: Combine these
        newphi->setVal<RunOn::Host>(0);
        newresid->setVal<RunOn::Host>(0);
        newbeta->setVal<RunOn::Host>(0);
   
        DEF_BOX_LIMITS(domain,dom_lo,dom_hi);
        DEF_BOX_LIMITS(newdomain,new_lo,new_hi);
        DEF_LIMITS(*beta,betaPtr,beta_lo,beta_hi);
        DEF_LIMITS(*newbeta,newbetaPtr,newbeta_lo,newbeta_hi);

        coarsig(betaPtr,ARLIM(beta_lo),ARLIM(beta_hi),
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

    resid->setVal<RunOn::Host>(0);

    hgresid (ARLIM(rhslo), ARLIM(rhshi),  rhsPtr,
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
    if (cgwork.isAllocated())
    {
        Real resnorm  = 0.0;

        FArrayBox dest0(phi->box(),1);

        DEF_BOX_LIMITS(domain,lo,hi);
        DEF_LIMITS(*phi,phiPtr,phi_lo,phi_hi);
        DEF_LIMITS(*resid,residPtr,resid_lo,resid_hi);
        DEF_LIMITS(dest0,dest0Ptr,dest0_lo,dest0_hi);
        DEF_LIMITS(*rhs,rhsPtr,rhs_lo,rhs_hi);
        DEF_LIMITS(*beta, betaPtr, beta_lo,beta_hi); 
        DEF_BOX_LIMITS(cgwork,cg_lo,cg_hi);

        solvehg(phiPtr,ARLIM(phi_lo),ARLIM(phi_hi),
                     dest0Ptr, ARLIM(dest0_lo),ARLIM(dest0_hi),
                     rhsPtr,ARLIM(rhs_lo),ARLIM(rhs_hi),
                     betaPtr, ARLIM(beta_lo),ARLIM(beta_hi),
                     cgwork.dataPtr(0),ARLIM(cg_lo),ARLIM(cg_hi),
                     cgwork.dataPtr(1), ARLIM(cg_lo),ARLIM(cg_hi),
                     cgwork.dataPtr(2), ARLIM(cg_lo),ARLIM(cg_hi),
                     cgwork.dataPtr(3),ARLIM(cg_lo),ARLIM(cg_hi),
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
ProjOutFlowBC_MG::Restrict ()
{
    DEF_BOX_LIMITS(domain,lo,hi);
    DEF_BOX_LIMITS(next->theDomain(),loc,hic);
    DEF_LIMITS(*resid,residPtr,resid_lo,resid_hi);
    DEF_LIMITS(*(next->theRhs()),rescPtr,resc_lo,resc_hi);

    next->theRhs()->setVal<RunOn::Host>(0);
  
    fort_restrict(residPtr,ARLIM(resid_lo),ARLIM(resid_hi), 
                  rescPtr, ARLIM(resc_lo),ARLIM(resc_hi),
                  lo,hi,loc,hic,isPeriodic);
}

void
ProjOutFlowBC_MG::interpolate ()
{
    FArrayBox temp(phi->box(),1);

    temp.setVal<RunOn::Host>(0);
    
    DEF_BOX_LIMITS(domain,lo,hi);
    DEF_BOX_LIMITS(next->theDomain(),loc,hic);
    DEF_LIMITS(*phi,phiPtr,phi_lo,phi_hi);
    DEF_LIMITS(temp,tempPtr,temp_lo,temp_hi);
    DEF_LIMITS(*(next->thePhi()),deltacPtr,deltac_lo,deltac_hi);
    DEF_LIMITS(*beta,betaPtr,beta_lo,beta_hi);

    interp(phiPtr, ARLIM(phi_lo),ARLIM(phi_hi),
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

    hgrelax(ARLIM(rhslo),ARLIM(rhshi),rhsPtr,
                 ARLIM(betalo),ARLIM(betahi),betaPtr,
                 ARLIM(philo),ARLIM(phihi),phiPtr,
                 ARLIM(dglo),ARLIM(dghi),dgphiPtr,
                 lo,hi,h,isPeriodic,&nstep);
}
#endif
