#ifndef _TRANSVERSE_3D_K_H_
#define _TRANSVERSE_3D_K_H_

#include <EBGodunov.H>
#include <AMReX_EBMultiFabUtil_3D_C.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void create_transverse_terms_for_xface(int ic, int j, int k, 
                                       Array4<Real const> const& v_ad,
                                       Array4<Real const> const& w_ad,
                                       Array4<Real const> const& yzlo,
                                       Array4<Real const> const& zylo,
                                       Array4<Real const> const& apy,
                                       Array4<Real const> const& apz,
                                       Array4<Real const> const& fcy,
                                       Array4<Real const> const& fcz,
                                       Real& trans_y, Real& trans_z,
                                       const Real l_dt, const Real dy, const Real dz)
{
    Real v_tmp_j, v_tmp_jp1;
    Real w_tmp_k, w_tmp_kp1;
    Real y_hat_j, y_hat_jp1;
    Real z_hat_k, z_hat_kp1;

    int ii, jj, kk;

    // Tangential extrapolation in the x- and z-directions at j-1/2
    if (apy(ic,j,k) == 1.0)
    {
        v_tmp_j = v_ad(ic,j,k);
        y_hat_j = yzlo(ic,j,k);

    } else {

        int jc = j;
        int kc = k;

        int ii = (fcy(ic,jc,kc,0) <= 0.0) ? ic-1 : ic+1;
        int kk = (fcy(ic,jc,kc,1) <= 0.0) ? kc-1 : kc+1;

        Real xoff = (fcy(ic,jc,kc,0) <= 0.0) ? -1. : 1.;
        Real zoff = (fcy(ic,jc,kc,1) <= 0.0) ? -1. : 1.;

        Real x0 =  fcy(ic,jc,kc,0);
        Real z0 =  fcy(ic,jc,kc,1);
        Real v0 = v_ad(ic,jc,kc  ); Real h0 = yzlo(ic,jc,kc);

        Real x1 =  fcy(ii,jc,kc,0)+xoff;
        Real z1 =  fcy(ii,jc,kc,1);
        Real v1 = v_ad(ii,jc,kc  ); Real h1 = yzlo(ii,jc,kc);

        Real x2 =  fcy(ic,jc,kk,0);
        Real z2 =  fcy(ic,jc,kk,1)+zoff;
        Real v2 = v_ad(ic,jc,kk  ); Real h2 = yzlo(ic,jc,kk);

        Real x3 =  fcy(ii,jc,kk,0)+xoff;
        Real z3 =  fcy(ii,jc,kk,1)+zoff;
        Real v3 = v_ad(ii,jc,kk  ); Real h3 = yzlo(ii,jc,kk);

        v_tmp_j = EB_interp_in_quad(0.,0.,v0,v1,v2,v3,x0,z0,x1,z1,x2,z2,x3,z3);
        y_hat_j = EB_interp_in_quad(0.,0.,h0,h1,h2,h3,x0,z0,x1,z1,x2,z2,x3,z3);
    }

    // Tangential extrapolation in the x- and z-directions at j+1/2
    if (apy(ic,j+1,k) == 1.0)
    {
        v_tmp_jp1 = v_ad(ic,j+1,k);
        y_hat_jp1 = yzlo(ic,j+1,k);

    } else {

        int jc = j+1;
        int kc = k;

        int ii = (fcy(ic,jc,kc,0) <= 0.0) ? ic-1 : ic+1;
        int kk = (fcy(ic,jc,kc,1) <= 0.0) ? kc-1 : kc+1;

        Real xoff = (fcy(ic,jc,kc,0) <= 0.0) ? -1. : 1.;
        Real zoff = (fcy(ic,jc,kc,1) <= 0.0) ? -1. : 1.;

        Real x0 =  fcy(ic,jc,kc,0);
        Real z0 =  fcy(ic,jc,kc,1);
        Real v0 = v_ad(ic,jc,kc  ); Real h0 = yzlo(ic,jc,kc);

        Real x1 =  fcy(ii,jc,kc,0)+xoff;
        Real z1 =  fcy(ii,jc,kc,1);
        Real v1 = v_ad(ii,jc,kc  ); Real h1 = yzlo(ii,jc,kc);

        Real x2 =  fcy(ic,jc,kk,0);
        Real z2 =  fcy(ic,jc,kk,1)+zoff;
        Real v2 = v_ad(ic,jc,kk  ); Real h2 = yzlo(ic,jc,kk);

        Real x3 =  fcy(ii,jc,kk,0)+xoff;
        Real z3 =  fcy(ii,jc,kk,1)+zoff;
        Real v3 = v_ad(ii,jc,kk  ); Real h3 = yzlo(ii,jc,kk);

        v_tmp_jp1 = EB_interp_in_quad(0.,0.,v0,v1,v2,v3,x0,z0,x1,z1,x2,z2,x3,z3);
        y_hat_jp1 = EB_interp_in_quad(0.,0.,h0,h1,h2,h3,x0,z0,x1,z1,x2,z2,x3,z3);
    }

    // Tangential extrapolation in the x- and y-directions at k-1/2
    if (apz(ic,j,k) == 1.0)
    {
        w_tmp_k = w_ad(ic,j,k);
        z_hat_k = zylo(ic,j,k);

    } else {

        int jc = j;
        int kc = k;

        int ii = (fcz(ic,jc,kc,0) <= 0.0) ? ic-1 : ic+1;
        int jj = (fcz(ic,jc,kc,1) <= 0.0) ? jc-1 : jc+1;

        Real xoff = (fcz(ic,jc,kc,0) <= 0.0) ? -1. : 1.;
        Real yoff = (fcz(ic,jc,kc,1) <= 0.0) ? -1. : 1.;

        Real x0 =  fcz(ic,jc,kc,0);
        Real y0 =  fcz(ic,jc,kc,1);
        Real v0 = w_ad(ic,jc,kc  ); Real h0 = zylo(ic,jc,kc);

        Real x1 =  fcz(ii,jc,kc,0)+xoff;
        Real y1 =  fcz(ii,jc,kc,1);
        Real v1 = w_ad(ii,jc,kc  ); Real h1 = zylo(ii,jc,kc);

        Real x2 =  fcz(ic,jj,kc,0);
        Real y2 =  fcz(ic,jj,kc,1)+yoff;
        Real v2 = w_ad(ic,jj,kc  ); Real h2 = zylo(ic,jj,kc);

        Real x3 =  fcz(ii,jj,kc,0)+xoff;
        Real y3 =  fcz(ii,jj,kc,1)+yoff;
        Real v3 = w_ad(ii,jj,kc  ); Real h3 = zylo(ii,jj,kc);

        w_tmp_k = EB_interp_in_quad(0.,0.,v0,v1,v2,v3,x0,y0,x1,y1,x2,y2,x3,y3);
        z_hat_k = EB_interp_in_quad(0.,0.,h0,h1,h2,h3,x0,y0,x1,y1,x2,y2,x3,y3);
    }

    // Tangential extrapolation in the x- and y-directions at k+1/2
    if (apz(ic,j,k+1) == 1.0)
    {
        w_tmp_kp1 = w_ad(ic,j,k+1);
        z_hat_kp1 = zylo(ic,j,k+1);

    } else {

        int jc = j;
        int kc = k+1;

        int ii = (fcz(ic,jc,kc,0) <= 0.0) ? ic-1 : ic+1;
        int jj = (fcz(ic,jc,kc,1) <= 0.0) ? jc-1 : jc+1;

        Real xoff = (fcz(ic,jc,kc,0) <= 0.0) ? -1. : 1.;
        Real yoff = (fcz(ic,jc,kc,1) <= 0.0) ? -1. : 1.;

        Real x0 =  fcz(ic,jc,kc,0);
        Real y0 =  fcz(ic,jc,kc,1);
        Real v0 = w_ad(ic,jc,kc  ); Real h0 = zylo(ic,jc,kc);

        Real x1 =  fcz(ii,jc,kc,0)+xoff;
        Real y1 =  fcz(ii,jc,kc,1);
        Real v1 = w_ad(ii,jc,kc  ); Real h1 = zylo(ii,jc,kc);

        Real x2 =  fcz(ic,jj,kc,0);
        Real y2 =  fcz(ic,jj,kc,1)+yoff;
        Real v2 = w_ad(ic,jj,kc  ); Real h2 = zylo(ic,jj,kc);

        Real x3 =  fcz(ii,jj,kc,0)+xoff;
        Real y3 =  fcz(ii,jj,kc,1)+yoff;
        Real v3 = w_ad(ii,jj,kc  ); Real h3 = zylo(ii,jj,kc);

        w_tmp_kp1 = EB_interp_in_quad(0.,0.,v0,v1,v2,v3,x0,y0,x1,y1,x2,y2,x3,y3);
        z_hat_kp1 = EB_interp_in_quad(0.,0.,h0,h1,h2,h3,x0,y0,x1,y1,x2,y2,x3,y3);
    }

    Real v_tmp_l = 0.5 * (v_tmp_jp1 + v_tmp_j);
    Real  dudy_l =       (y_hat_jp1 - y_hat_j);

    Real w_tmp_l = 0.5 * (w_tmp_kp1 + w_tmp_k);
    Real  dudz_l =       (z_hat_kp1 - z_hat_k);

    trans_y = v_tmp_l * dudy_l / dy;
    trans_z = w_tmp_l * dudz_l / dz;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void create_transverse_terms_for_yface(int i, int jc, int k, 
                                       Array4<Real const> const& u_ad,
                                       Array4<Real const> const& w_ad,
                                       Array4<Real const> const& xzlo,
                                       Array4<Real const> const& zxlo,
                                       Array4<Real const> const& apx,
                                       Array4<Real const> const& apz,
                                       Array4<Real const> const& fcx,
                                       Array4<Real const> const& fcz,
                                       Real& trans_x, Real& trans_z,
                                       const Real l_dt, const Real dx, const Real dz)
{
    Real u_tmp_i, u_tmp_ip1;
    Real w_tmp_k, w_tmp_kp1;
    Real x_hat_i, x_hat_ip1;
    Real z_hat_k, z_hat_kp1;

    int ii, jj, kk;

    // Tangential extrapolation in the y- and z-directions at i-1/2
    if (apx(i,jc,k) == 1.0)
    {
        u_tmp_i = u_ad(i,jc,k);
        x_hat_i = xzlo(i,jc,k);

    } else {

        int ic = i;
        int kc = k;

        int jj = (fcx(ic,jc,kc,0) <= 0.0) ? jc-1 : jc+1;
        int kk = (fcx(ic,jc,kc,1) <= 0.0) ? kc-1 : kc+1;

        Real yoff = (fcx(ic,jc,kc,0) <= 0.0) ? -1. : 1.;
        Real zoff = (fcx(ic,jc,kc,1) <= 0.0) ? -1. : 1.;

        Real y0 =  fcx(ic,jc,kc,0);
        Real z0 =  fcx(ic,jc,kc,1);
        Real v0 = u_ad(ic,jc,kc  ); Real h0 = xzlo(ic,jc,kc);

        Real y1 =  fcx(ic,jj,kc,0)+yoff;
        Real z1 =  fcx(ic,jj,kc,1);
        Real v1 = u_ad(ic,jj,kc  ); Real h1 = xzlo(ic,jj,kc);

        Real y2 =  fcx(ic,jc,kk,0);
        Real z2 =  fcx(ic,jc,kk,1)+zoff;
        Real v2 = u_ad(ic,jc,kk  ); Real h2 = xzlo(ic,jc,kk);

        Real y3 =  fcx(ic,jj,kk,0)+yoff;
        Real z3 =  fcx(ic,jj,kk,1)+zoff;
        Real v3 = u_ad(ic,jj,kk  ); Real h3 = xzlo(ic,jj,kk);

        u_tmp_i = EB_interp_in_quad(0.,0.,v0,v1,v2,v3,y0,z0,y1,z1,y2,z2,y3,z3);
        x_hat_i = EB_interp_in_quad(0.,0.,h0,h1,h2,h3,y0,z0,y1,z1,y2,z2,y3,z3);
    }

    // Tangential extrapolation in the y- and z-directions at i+1/2
    if (apx(i+1,jc,k) == 1.0)
    {
        u_tmp_ip1 = u_ad(i+1,jc,k);
        x_hat_ip1 = xzlo(i+1,jc,k);

    } else {

        int ic = i+1;
        int kc = k;

        int jj = (fcx(ic,jc,kc,0) <= 0.0) ? jc-1 : jc+1;
        int kk = (fcx(ic,jc,kc,1) <= 0.0) ? kc-1 : kc+1;

        Real yoff = (fcx(ic,jc,kc,0) <= 0.0) ? -1. : 1.;
        Real zoff = (fcx(ic,jc,kc,1) <= 0.0) ? -1. : 1.;

        Real y0 =  fcx(ic,jc,kc,0);
        Real z0 =  fcx(ic,jc,kc,1);
        Real v0 = u_ad(ic,jc,kc  ); Real h0 = xzlo(ic,jc,kc);

        Real y1 =  fcx(ic,jj,kc,0)+yoff;
        Real z1 =  fcx(ic,jj,kc,1);
        Real v1 = u_ad(ic,jj,kc  ); Real h1 = xzlo(ic,jj,kc);

        Real y2 =  fcx(ic,jc,kk,0);
        Real z2 =  fcx(ic,jc,kk,1)+zoff;
        Real v2 = u_ad(ic,jc,kk  ); Real h2 = xzlo(ic,jc,kk);

        Real y3 =  fcx(ic,jj,kk,0)+yoff;
        Real z3 =  fcx(ic,jj,kk,1)+zoff;
        Real v3 = u_ad(ic,jj,kk  ); Real h3 = xzlo(ic,jj,kk);

        u_tmp_ip1 = EB_interp_in_quad(0.,0.,v0,v1,v2,v3,y0,z0,y1,z1,y2,z2,y3,z3);
        x_hat_ip1 = EB_interp_in_quad(0.,0.,h0,h1,h2,h3,y0,z0,y1,z1,y2,z2,y3,z3);
    }

    // Tangential extrapolation in the x- and y-directions at k-1/2
    if (apz(i,jc,k) == 1.0)
    {
        w_tmp_k = w_ad(i,jc,k);
        z_hat_k = zxlo(i,jc,k);

    } else {

        int ic = i;
        int kc = k;

        int ii = (fcz(ic,jc,kc,0) <= 0.0) ? ic-1 : ic+1;
        int jj = (fcz(ic,jc,kc,1) <= 0.0) ? jc-1 : jc+1;

        Real xoff = (fcz(ic,jc,kc,0) <= 0.0) ? -1. : 1.;
        Real yoff = (fcz(ic,jc,kc,1) <= 0.0) ? -1. : 1.;

        Real x0 =  fcz(ic,jc,kc,0);
        Real y0 =  fcz(ic,jc,kc,1);
        Real v0 = w_ad(ic,jc,kc  ); Real h0 = zxlo(ic,jc,kc);

        Real x1 =  fcz(ii,jc,kc,0)+xoff;
        Real y1 =  fcz(ii,jc,kc,1);
        Real v1 = w_ad(ii,jc,kc  ); Real h1 = zxlo(ii,jc,kc);

        Real x2 =  fcz(ic,jj,kc,0);
        Real y2 =  fcz(ic,jj,kc,1)+yoff;
        Real v2 = w_ad(ic,jj,kc  ); Real h2 = zxlo(ic,jj,kc);

        Real x3 =  fcz(ii,jj,kc,0)+xoff;
        Real y3 =  fcz(ii,jj,kc,1)+yoff;
        Real v3 = w_ad(ii,jj,kc  ); Real h3 = zxlo(ii,jj,kc);

        w_tmp_k = EB_interp_in_quad(0.,0.,v0,v1,v2,v3,x0,y0,x1,y1,x2,y2,x3,y3);
        z_hat_k = EB_interp_in_quad(0.,0.,h0,h1,h2,h3,x0,y0,x1,y1,x2,y2,x3,y3);
    }

    // Tangential extrapolation in the x- and y-directions at k+1/2
    if (apz(i,jc,k+1) == 1.0)
    {
        w_tmp_kp1 = w_ad(i,jc,k+1);
        z_hat_kp1 = zxlo(i,jc,k+1);

    } else {

        int ic = i;
        int kc = k+1;

        int ii = (fcz(ic,jc,kc,0) <= 0.0) ? ic-1 : ic+1;
        int jj = (fcz(ic,jc,kc,1) <= 0.0) ? jc-1 : jc+1;

        Real xoff = (fcz(ic,jc,kc,0) <= 0.0) ? -1. : 1.;
        Real yoff = (fcz(ic,jc,kc,1) <= 0.0) ? -1. : 1.;

        Real x0 =  fcz(ic,jc,kc,0);
        Real y0 =  fcz(ic,jc,kc,1);
        Real v0 = w_ad(ic,jc,kc  ); Real h0 = zxlo(ic,jc,kc);

        Real x1 =  fcz(ii,jc,kc,0)+xoff;
        Real y1 =  fcz(ii,jc,kc,1);
        Real v1 = w_ad(ii,jc,kc  ); Real h1 = zxlo(ii,jc,kc);

        Real x2 =  fcz(ic,jj,kc,0);
        Real y2 =  fcz(ic,jj,kc,1)+yoff;
        Real v2 = w_ad(ic,jj,kc  ); Real h2 = zxlo(ic,jj,kc);

        Real x3 =  fcz(ii,jj,kc,0)+xoff;
        Real y3 =  fcz(ii,jj,kc,1)+yoff;
        Real v3 = w_ad(ii,jj,kc  ); Real h3 = zxlo(ii,jj,kc);

        w_tmp_kp1 = EB_interp_in_quad(0.,0.,v0,v1,v2,v3,x0,y0,x1,y1,x2,y2,x3,y3);
        z_hat_kp1 = EB_interp_in_quad(0.,0.,h0,h1,h2,h3,x0,y0,x1,y1,x2,y2,x3,y3);
    }

    Real u_tmp_l = 0.5 * (u_tmp_ip1 + u_tmp_i);
    Real  dudx_l =       (x_hat_ip1 - x_hat_i);

    Real w_tmp_l = 0.5 * (w_tmp_kp1 + w_tmp_k);
    Real  dudz_l =       (z_hat_kp1 - z_hat_k);

    trans_x = u_tmp_l * dudx_l / dx;
    trans_z = w_tmp_l * dudz_l / dz;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void create_transverse_terms_for_zface(int i, int j, int kc, 
                                       Array4<Real const> const& u_ad,
                                       Array4<Real const> const& v_ad,
                                       Array4<Real const> const& xylo,
                                       Array4<Real const> const& yxlo,
                                       Array4<Real const> const& apx,
                                       Array4<Real const> const& apy,
                                       Array4<Real const> const& fcx,
                                       Array4<Real const> const& fcy,
                                       Real& trans_x, Real& trans_y,
                                       const Real l_dt, const Real dx, const Real dy)
{
    Real u_tmp_i, u_tmp_ip1;
    Real v_tmp_j, v_tmp_jp1;
    Real x_hat_i, x_hat_ip1;
    Real y_hat_j, y_hat_jp1;

    int ii, jj, kk;

    // Tangential extrapolation in the y- and z-directions at i-1/2
    if (apx(i,j,kc) == 1.0)
    {
        u_tmp_i = u_ad(i,j,kc);
        x_hat_i = xylo(i,j,kc);

    } else {

        int ic = i;
        int jc = j;

        int jj = (fcx(ic,jc,kc,0) <= 0.0) ? jc-1 : jc+1;
        int kk = (fcx(ic,jc,kc,1) <= 0.0) ? kc-1 : kc+1;

        Real yoff = (fcx(ic,jc,kc,0) <= 0.0) ? -1. : 1.;
        Real zoff = (fcx(ic,jc,kc,1) <= 0.0) ? -1. : 1.;

        Real y0 =  fcx(ic,jc,kc,0); 
        Real z0 =  fcx(ic,jc,kc,1); 
        Real v0 = u_ad(ic,jc,kc  ); Real h0 = xylo(ic,jc,kc  );

        Real y1 =  fcx(ic,jj,kc,0)+yoff;
        Real z1 =  fcx(ic,jj,kc,1);
        Real v1 = u_ad(ic,jj,kc  ); Real h1 = xylo(ic,jj,kc  );

        Real y2 =  fcx(ic,jc,kk,0);
        Real z2 =  fcx(ic,jc,kk,1)+zoff;
        Real v2 = u_ad(ic,jc,kk  ); Real h2 = xylo(ic,jc,kk  );

        Real y3 =  fcx(ic,jj,kk,0)+yoff;
        Real z3 =  fcx(ic,jj,kk,1)+zoff;
        Real v3 = u_ad(ic,jj,kk  ); Real h3 = xylo(ic,jj,kk  );

        u_tmp_i = EB_interp_in_quad(0.,0.,v0,v1,v2,v3,y0,z0,y1,z1,y2,z2,y3,z3);
        x_hat_i = EB_interp_in_quad(0.,0.,h0,h1,h2,h3,y0,z0,y1,z1,y2,z2,y3,z3);
    }

    // Tangential extrapolation in the y- and z-directions at i+1/2
    if (apx(i+1,j,kc) == 1.0)
    {
        u_tmp_ip1 = u_ad(i+1,j,kc);
        x_hat_ip1 = xylo(i+1,j,kc);

    } else {

        int ic = i+1;
        int jc = j;

        int jj = (fcx(ic,jc,kc,0) <= 0.0) ? jc-1 : jc+1;
        int kk = (fcx(ic,jc,kc,1) <= 0.0) ? kc-1 : kc+1;

        Real yoff = (fcx(ic,jc,kc,0) <= 0.0) ? -1. : 1.;
        Real zoff = (fcx(ic,jc,kc,1) <= 0.0) ? -1. : 1.;

        Real y0 =  fcx(ic,jc,kc,0); 
        Real z0 =  fcx(ic,jc,kc,1); 
        Real v0 = u_ad(ic,jc,kc  ); Real h0 = xylo(ic,jc,kc  );

        Real y1 =  fcx(ic,jj,kc,0)+yoff;
        Real z1 =  fcx(ic,jj,kc,1);
        Real v1 = u_ad(ic,jj,kc  ); Real h1 = xylo(ic,jj,kc  );

        Real y2 =  fcx(ic,jc,kk,0);
        Real z2 =  fcx(ic,jc,kk,1)+zoff;
        Real v2 = u_ad(ic,jc,kk  ); Real h2 = xylo(ic,jc,kk  );

        Real y3 =  fcx(ic,jj,kk,0)+yoff;
        Real z3 =  fcx(ic,jj,kk,1)+zoff;
        Real v3 = u_ad(ic,jj,kk  ); Real h3 = xylo(ic,jj,kk  );

        u_tmp_ip1 = EB_interp_in_quad(0.,0.,v0,v1,v2,v3,y0,z0,y1,z1,y2,z2,y3,z3);
        x_hat_ip1 = EB_interp_in_quad(0.,0.,h0,h1,h2,h3,y0,z0,y1,z1,y2,z2,y3,z3);
    }

    // Tangential extrapolation in the x- and z-directions at j-1/2
    if (apy(i,j,kc) == 1.0)
    {
        v_tmp_j = v_ad(i,j,kc);
        y_hat_j = yxlo(i,j,kc);

    } else {

        int ic = i;
        int jc = j;

        int ii = (fcy(ic,jc,kc,0) <= 0.0) ? ic-1 : ic+1;
        int kk = (fcy(ic,jc,kc,1) <= 0.0) ? kc-1 : kc+1;

        Real xoff = (fcy(ic,jc,kc,0) <= 0.0) ? -1. : 1.;
        Real zoff = (fcy(ic,jc,kc,1) <= 0.0) ? -1. : 1.;

        Real x0 =  fcy(ic,jc,kc,0);
        Real z0 =  fcy(ic,jc,kc,1);
        Real v0 = v_ad(ic,jc,kc  ); Real h0 = yxlo(ic,jc,kc);

        Real x1 =  fcy(ii,jc,kc,0)+xoff;
        Real z1 =  fcy(ii,jc,kc,1);
        Real v1 = v_ad(ii,jc,kc  ); Real h1 = yxlo(ii,jc,kc);

        Real x2 =  fcy(ic,jc,kk,0);
        Real z2 =  fcy(ic,jc,kk,1)+zoff;
        Real v2 = v_ad(ic,jc,kk  ); Real h2 = yxlo(ic,jc,kk);

        Real x3 =  fcy(ii,jc,kk,0)+xoff;
        Real z3 =  fcy(ii,jc,kk,1)+zoff;
        Real v3 = v_ad(ii,jc,kk  ); Real h3 = yxlo(ii,jc,kk);

        v_tmp_j = EB_interp_in_quad(0.,0.,v0,v1,v2,v3,x0,z0,x1,z1,x2,z2,x3,z3);
        y_hat_j = EB_interp_in_quad(0.,0.,h0,h1,h2,h3,x0,z0,x1,z1,x2,z2,x3,z3);
    }

    // Tangential extrapolation in the x- and z-directions at j+1/2
    if (apy(i,j+1,kc) == 1.0)
    {
        v_tmp_jp1 = v_ad(i,j+1,kc);
        y_hat_jp1 = yxlo(i,j+1,kc);

    } else {

        int ic = i;
        int jc = j+1;

        int ii = (fcy(ic,jc,kc,0) <= 0.0) ? ic-1 : ic+1;
        int kk = (fcy(ic,jc,kc,1) <= 0.0) ? kc-1 : kc+1;

        Real xoff = (fcy(ic,jc,kc,0) <= 0.0) ? -1. : 1.;
        Real zoff = (fcy(ic,jc,kc,1) <= 0.0) ? -1. : 1.;

        Real x0 =  fcy(ic,jc,kc,0);
        Real z0 =  fcy(ic,jc,kc,1);
        Real v0 = v_ad(ic,jc,kc  ); Real h0 = yxlo(ic,jc,kc);

        Real x1 =  fcy(ii,jc,kc,0)+xoff;
        Real z1 =  fcy(ii,jc,kc,1);
        Real v1 = v_ad(ii,jc,kc  ); Real h1 = yxlo(ii,jc,kc);

        Real x2 =  fcy(ic,jc,kk,0);
        Real z2 =  fcy(ic,jc,kk,1)+zoff;
        Real v2 = v_ad(ic,jc,kk  ); Real h2 = yxlo(ic,jc,kk);

        Real x3 =  fcy(ii,jc,kk,0)+xoff;
        Real z3 =  fcy(ii,jc,kk,1)+zoff;
        Real v3 = v_ad(ii,jc,kk  ); Real h3 = yxlo(ii,jc,kk);

        v_tmp_jp1 = EB_interp_in_quad(0.,0.,v0,v1,v2,v3,x0,z0,x1,z1,x2,z2,x3,z3);
        y_hat_jp1 = EB_interp_in_quad(0.,0.,h0,h1,h2,h3,x0,z0,x1,z1,x2,z2,x3,z3);
    }

    Real u_tmp_l = 0.5 * (u_tmp_ip1 + u_tmp_i);
    Real  dudx_l =       (x_hat_ip1 - x_hat_i);

    Real v_tmp_l = 0.5 * (v_tmp_jp1 + v_tmp_j);
    Real  dudy_l =       (y_hat_jp1 - y_hat_j);

    trans_x = u_tmp_l * dudx_l / dx;
    trans_y = v_tmp_l * dudy_l / dy;
}
#endif
