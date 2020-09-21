#include <iamr_ppm_godunov.H>

using namespace amrex;

void
PPM::PredictVelOnFaces (Box const& bx, int ncomp,
                        AMREX_D_DECL( Array4<Real> const& Imx,
                                      Array4<Real> const& Imy,
                                      Array4<Real> const& Imz),
                        AMREX_D_DECL( Array4<Real> const& Ipx,
                                      Array4<Real> const& Ipy,
                                      Array4<Real> const& Ipz),
                        Array4<Real const> const& q,
                        Array4<Real const> const& vel,
                        Geometry geom,
                        Real dt,
                        BCRec const* pbc)
{
    const Box& domain = geom.Domain();
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);

    const auto dx = geom.CellSizeArray();
    AMREX_D_TERM( Real l_dtdx = dt / dx[0];,
                  Real l_dtdy = dt / dx[1];,
                  Real l_dtdz = dt / dx[2];);

    amrex::ParallelFor(bx, AMREX_SPACEDIM,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        PredictVelOnXFace(i,j,k,n,l_dtdx,vel(i,j,k,0),q,Imx,Ipx,pbc[n],dlo.x,dhi.x);
        PredictVelOnYFace(i,j,k,n,l_dtdy,vel(i,j,k,1),q,Imy,Ipy,pbc[n],dlo.y,dhi.y);
#if (AMREX_SPACEDIM==3)
        PredictVelOnZFace(i,j,k,n,l_dtdz,vel(i,j,k,2),q,Imz,Ipz,pbc[n],dlo.z,dhi.z);
#endif
    });
}
