#ifndef _TRANSVERSE_2D_K_H_
#define _TRANSVERSE_2D_K_H_

#include <EBGodunov.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
linear_interp (Real xout,
               Real x1in, Real val1,
               Real x2in, Real val2)
{
   Real val_out =  val1 + (xout-x1in)/(x2in-x1in) * (val2-val1);
   return val_out;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void create_transverse_terms_for_xface(int ic, int j, int k, 
                                       Array4<Real const> const& v_ad,
                                       Array4<Real const> const& yhat,
                                       Array4<Real const> const& apy,
                                       Array4<Real const> const& fcy,
                                       Real& trans_y, 
                                       const Real l_dt, const Real dy)
{
    Real v_tmp_j, v_tmp_jp1;
    Real y_hat_j, y_hat_jp1;

    // Tangential extrapolation in the x-direction at j-1/2
    if (apy(ic,j,k) == 1.0)
    {
        v_tmp_j = v_ad(ic,j,k);
        y_hat_j = yhat(ic,j,k);
    } 
    else 
    {
        // Tangential extrapolation in x-direction at j-1/2
        int jc    = j;
        int ii    = (fcy(ic,jc,k,0) <= 0.0) ? ic-1 : ic+1;
        Real xoff = (fcy(ic,jc,k,0) <= 0.0) ?  -1. : 1.;

        v_tmp_j = linear_interp(0.0, fcy(ic,jc,k), v_ad(ic,jc,k), (xoff+fcy(ii,jc,k)), v_ad(ii,jc,k));
        y_hat_j = linear_interp(0.0, fcy(ic,jc,k), yhat(ic,jc,k), (xoff+fcy(ii,jc,k)), yhat(ii,jc,k));
    }

    // Tangential extrapolation in the x-direction at j+1/2
    if (apy(ic,j+1,k) == 1.0)
    {
        v_tmp_jp1 = v_ad(ic,j+1,k);
        y_hat_jp1 = yhat(ic,j+1,k);
    }
    else 
    {
        // Tangential extrapolation in x-direction at j+1/2
        int jc    = j+1;
        int ii    = (fcy(ic,jc,k,0) <= 0.0) ? ic-1 : ic+1;
        Real xoff = (fcy(ic,jc,k,0) <= 0.0) ?  -1. : 1.;

        v_tmp_jp1 = linear_interp(0.0, fcy(ic,jc,k), v_ad(ic,jc,k), (xoff+fcy(ii,jc,k)), v_ad(ii,jc,k));
        y_hat_jp1 = linear_interp(0.0, fcy(ic,jc,k), yhat(ic,jc,k), (xoff+fcy(ii,jc,k)), yhat(ii,jc,k));
    }

    trans_y = 0.5 * (v_tmp_jp1 + v_tmp_j) * (y_hat_jp1 - y_hat_j);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void create_transverse_terms_for_yface(int i, int jc, int k, 
                                       Array4<Real const> const& u_ad,
                                       Array4<Real const> const& xhat,
                                       Array4<Real const> const& apx,
                                       Array4<Real const> const& fcx,
                                       Real& trans_x, 
                                       const Real l_dt, const Real dx)
{
    Real u_tmp_i, u_tmp_ip1;
    Real x_hat_i, x_hat_ip1;

    // Tangential extrapolation in the y-direction at i-1/2
    if (apx(i,jc,k) == 1.0)
    {
        u_tmp_i = u_ad(i,jc,k);
        x_hat_i = xhat(i,jc,k);
    } 
    else 
    {
        // Tangential extrapolation in y-direction at i-1/2
        int ic    = i;
        int jj    = (fcx(ic,jc,k,0) <= 0.0) ? jc-1 : jc+1;
        Real yoff = (fcx(ic,jc,k,0) <= 0.0) ?  -1. : 1.;

        u_tmp_i = linear_interp(0.0, fcx(ic,jc,k), u_ad(ic,jc,k), (yoff+fcx(ic,jj,k)), u_ad(ic,jj,k));
        x_hat_i = linear_interp(0.0, fcx(ic,jc,k), xhat(ic,jc,k), (yoff+fcx(ic,jj,k)), xhat(ic,jj,k));
    }

    // Tangential extrapolation in the x-direction at j+1/2
    if (apx(i+1,jc,k) == 1.0)
    {
        u_tmp_ip1 = u_ad(i+1,jc,k);
        x_hat_ip1 = xhat(i+1,jc,k);
    }
    else 
    {
        // Tangential extrapolation in y-direction at i+1/2
        int ic = i+1;
        int jj = (fcx(ic,jc,k,0) <= 0.0) ? jc-1 : jc+1;
        Real yoff = (fcx(ic,jc,k,0) <= 0.0) ?  -1. : 1.;

        u_tmp_ip1 = linear_interp(0.0, fcx(ic,jc,k), u_ad(ic,jc,k), (yoff+fcx(ic,jj,k)), u_ad(ic,jj,k));
        x_hat_ip1 = linear_interp(0.0, fcx(ic,jc,k), xhat(ic,jc,k), (yoff+fcx(ic,jj,k)), xhat(ic,jj,k));
    }

    trans_x = 0.5 * (u_tmp_ip1 + u_tmp_i) * (x_hat_ip1 - x_hat_i);
}

#endif
