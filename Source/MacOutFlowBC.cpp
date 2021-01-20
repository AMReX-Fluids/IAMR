

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

namespace
{
    bool outflow_initialized    = false;

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

    for (iface = 0; iface < numOutFlowFaces; iface++)
    {
      phiMF[iface].setVal<RunOn::Gpu>(0);
    }
}
