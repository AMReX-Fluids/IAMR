

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

    for (int i=0; i < numOutFlowFaces; i++)
        phiMF[i].setVal<RunOn::Gpu>(0.);

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

//void
//ProjOutFlowBC::rhogbc(){
// TODO: port code from PROJOUTFLOW_3D.F90
//                const Box& bx = ;
//                const auto& du       = divu->array();
//                const auto& dsdt_arr = dsdt->array();
//                const auto& rhcc_arr = rhcclev->array();
//                amrex::ParallelFor(bx, [du,dsdt_arr,rhcc_arr,dt_inv]
//                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//                {
//                   // do stuff
//                });
//             }
// }


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
