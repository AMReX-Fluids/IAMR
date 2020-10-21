#include <iamr_constants.H>
#include <iamr_mol.H>
#include <AMReX_Slopes_K.H>

using namespace amrex;

void
MOL::PredictVelOnFaces (  D_DECL( Box const& ubx,
                                  Box const& vbx,
                                  Box const& wbx),
                          D_DECL( Array4<Real> const& u,
                                  Array4<Real> const& v,
                                  Array4<Real> const& w),
                          Array4<Real const> const& vcc,
                          const Geometry&  geom,
                          const BCRec* bc )
{

    const Box& domain_box = geom.Domain();
    const int  domain_ilo = domain_box.smallEnd(0);
    const int  domain_ihi = domain_box.bigEnd(0);
    const int  domain_jlo = domain_box.smallEnd(1);
    const int  domain_jhi = domain_box.bigEnd(1);
#if (AMREX_SPACEDIM==3)
    const int  domain_klo = domain_box.smallEnd(2);
    const int  domain_khi = domain_box.bigEnd(2);
#endif

    bool extdir_or_ho_ilo = (bc[0].lo(0) == BCType::ext_dir) || (bc[0].lo(0) == BCType::hoextrap);
    bool extdir_or_ho_ihi = (bc[0].hi(0) == BCType::ext_dir) || (bc[0].hi(0) == BCType::hoextrap);
    bool extdir_or_ho_jlo = (bc[1].lo(1) == BCType::ext_dir) || (bc[1].lo(1) == BCType::hoextrap);
    bool extdir_or_ho_jhi = (bc[1].lo(1) == BCType::ext_dir) || (bc[1].lo(1) == BCType::hoextrap);
#if (AMREX_SPACEDIM==3)
    bool extdir_or_ho_klo = (bc[2].lo(2) == BCType::ext_dir) || (bc[2].lo(2) == BCType::hoextrap);
    bool extdir_or_ho_khi = (bc[2].lo(2) == BCType::ext_dir) || (bc[2].lo(2) == BCType::hoextrap);
#endif

    bool extdir_ilo = (bc[0].lo(0) == BCType::ext_dir);
    bool extdir_ihi = (bc[0].hi(0) == BCType::ext_dir);
    bool extdir_jlo = (bc[1].lo(1) == BCType::ext_dir);
    bool extdir_jhi = (bc[1].lo(1) == BCType::ext_dir);
#if (AMREX_SPACEDIM==3)
    bool extdir_klo = (bc[2].lo(2) == BCType::ext_dir);
    bool extdir_khi = (bc[2].lo(2) == BCType::ext_dir);
#endif

    //
    // X direction
    //

    //slope order
    int order = 2;

    // At an ext_dir boundary, the boundary value is on the face, not cell center
    if ((extdir_or_ho_ilo and domain_ilo >= ubx.smallEnd(0)-1) or
        (extdir_or_ho_ihi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(ubx, [vcc,order,extdir_ilo,extdir_ihi,extdir_or_ho_ilo,extdir_or_ho_ihi,domain_ilo,domain_ihi,u]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
	    if (extdir_ilo and i == domain_ilo) {
	        u(i,j,k) = vcc(i-1,j,k,0);
            } else if (extdir_ihi and i == domain_ihi+1) {
                u(i,j,k) = vcc(i,j,k,0);
            } else {

	      Real upls = vcc(i,j,k,0) - 0.5 * amrex_calc_xslope_extdir
		(i,j,k,0,order,vcc, extdir_or_ho_ilo, extdir_or_ho_ihi, domain_ilo, domain_ihi);
	      Real umns = vcc(i-1,j,k,0) + 0.5 * amrex_calc_xslope_extdir
		(i-1,j,k,0,order,vcc, extdir_or_ho_ilo, extdir_or_ho_ihi, domain_ilo, domain_ihi);
	      if (umns < 0.0 and upls > 0.0) {
                u(i,j,k) = 0.0;
	      } else {
                Real avg = 0.5 * (upls + umns);
                if (std::abs(avg) < small_vel) {
                    u(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    u(i,j,k) = umns;
                } else {
                    u(i,j,k) = upls;
                }
	      }
            }
        });
    }
    else
    {
        amrex::ParallelFor(ubx, [vcc,order,u]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real upls = vcc(i  ,j,k,0) - 0.5 * amrex_calc_xslope(i  ,j,k,0,order,vcc);
	    Real umns = vcc(i-1,j,k,0) + 0.5 * amrex_calc_xslope(i-1,j,k,0,order,vcc);
            if (umns < 0.0 and upls > 0.0) {
                u(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (upls + umns);
                if (std::abs(avg) < small_vel) {
                    u(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    u(i,j,k) = umns;
                } else {
                    u(i,j,k) = upls;
                }
            }
        });
    }


    //
    // Y direction
    //

    if ((extdir_or_ho_jlo and domain_jlo >= vbx.smallEnd(1)-1) or
        (extdir_or_ho_jhi and domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(vbx, [vcc,order,extdir_jlo,extdir_jhi,extdir_or_ho_jlo,extdir_or_ho_jhi,domain_jlo,domain_jhi,v]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
	    if (extdir_jlo and j == domain_jlo) {
	        v(i,j,k) = vcc(i,j-1,k,1);
            } else if (extdir_jhi and j == domain_jhi+1) {
                v(i,j,k) = vcc(i,j,k,1);
            } else {

	      Real vpls = vcc(i,j,k,1) - 0.5 * amrex_calc_yslope_extdir
		(i,j,k,1,order,vcc, extdir_or_ho_jlo, extdir_or_ho_jhi, domain_jlo, domain_jhi);
	      Real vmns = vcc(i,j-1,k,1) + 0.5 * amrex_calc_yslope_extdir
                (i,j-1,k,1,order,vcc, extdir_or_ho_jlo, extdir_or_ho_jhi, domain_jlo, domain_jhi);
            if (vmns < 0.0 and vpls > 0.0) {
                v(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (vpls + vmns);
                if (std::abs(avg) < small_vel) {
                    v(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    v(i,j,k) = vmns;
                } else {
                    v(i,j,k) = vpls;
                }
            }
            }
        });
    }
    else
    {
        amrex::ParallelFor(vbx, [vcc,order,v]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
	    Real vpls = vcc(i,j  ,k,1) - 0.5 * amrex_calc_yslope(i,j  ,k,1,order,vcc);
            Real vmns = vcc(i,j-1,k,1) + 0.5 * amrex_calc_yslope(i,j-1,k,1,order,vcc);
            if (vmns < 0.0 and vpls > 0.0) {
                v(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (vpls + vmns);
                if (std::abs(avg) < small_vel) {
                    v(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    v(i,j,k) = vmns;
                } else {
                    v(i,j,k) = vpls;
                }
            }
        });
    }

#if (AMREX_SPACEDIM==3)
    //
    // Y direction
    //

    if ((extdir_or_ho_klo and domain_klo >= wbx.smallEnd(2)-1) or
        (extdir_or_ho_khi and domain_khi <= wbx.bigEnd(2)))
    {
      amrex::ParallelFor(wbx, [vcc,order,extdir_klo,extdir_khi,extdir_or_ho_klo,extdir_or_ho_khi,domain_klo,domain_khi,w]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
	    if (extdir_klo and k == domain_klo) {
                w(i,j,k) = vcc(i,j,k-1,2);
            } else if (extdir_khi and k == domain_khi+1) {
                w(i,j,k) = vcc(i,j,k,2);
            } else {

	      Real wpls = vcc(i,j,k,2) - 0.5 * amrex_calc_zslope_extdir
                (i,j,k,2,order,vcc, extdir_or_ho_klo, extdir_or_ho_khi, domain_klo, domain_khi);
	      Real wmns = vcc(i,j,k-1,2) + 0.5 * amrex_calc_zslope_extdir
		(i,j,k-1,2,order,vcc, extdir_or_ho_klo, extdir_or_ho_khi, domain_klo, domain_khi);
            if (wmns < 0.0 and wpls > 0.0) {
                w(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (wpls + wmns);
                if (std::abs(avg) < small_vel) {
                    w(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    w(i,j,k) = wmns;
                } else {
                    w(i,j,k) = wpls;
                }
            }
            }
        });
    }
    else
    {
        amrex::ParallelFor(wbx, [vcc,order,w]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real wpls = vcc(i,j,k  ,2) - 0.5 * amrex_calc_zslope(i,j,k  ,2,order,vcc);
            Real wmns = vcc(i,j,k-1,2) + 0.5 * amrex_calc_zslope(i,j,k-1,2,order,vcc);
            if (wmns < 0.0 and wpls > 0.0) {
                w(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (wpls + wmns);
                if (std::abs(avg) < small_vel) {
                    w(i,j,k) = 0.0;
            } else if (avg > 0.0) {
                    w(i,j,k) = wmns;
                } else {
                    w(i,j,k) = wpls;
                }
            }
        });
    }
#endif
}
