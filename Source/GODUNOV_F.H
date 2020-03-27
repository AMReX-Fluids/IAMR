
#ifndef _GODUNOV_F_H_
#define _GODUNOV_F_H_

#ifndef BL_LANG_FORT

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
   void fort_test_umac_rho(const amrex::Real* um, ARLIM_P(ulo), ARLIM_P(uhi),
			   const amrex::Real* vm, ARLIM_P(vlo), ARLIM_P(vhi),
#if (BL_SPACEDIM == 3)
			   const amrex::Real* wm, ARLIM_P(wlo), ARLIM_P(whi),
#endif
			   const amrex::Real* rh, ARLIM_P(rlo), ARLIM_P(rhi),
			   const int* lo, const int* hi, const amrex::Real* dt,
			   const amrex::Real* dx, amrex::Real* cfl, const amrex::Real* u_max);


   void fort_maxchng_velmag(const amrex::Real* Uodat,  ARLIM_P(uo_lo), ARLIM_P(uo_hi),
		   const amrex::Real* Undat,  ARLIM_P(un_lo), ARLIM_P(un_hi),
		   const int* lo, const int* hi, amrex::Real* max_change);


   void extrap_vel_to_faces(const int* lo, const int* hi,
                            const amrex::Real* u_dat, const int* u_lo, const int* u_hi,
                            const int* ubc,
                            const amrex::Real* tfx_dat, const int* tfx_lo, const int* tfx_hi,
                            amrex::Real*       umac_dat, const int* umac_lo, const int* umac_hi,
                            const int* vbc,
                            const amrex::Real* tfy_dat, const int* tfy_lo, const int* tfy_hi,
                            amrex::Real*       vmac_dat, const int* vmac_lo, const int* vmac_hi,
#if (AMREX_SPACEDIM == 3)
                            const int* wbc,
                            const amrex::Real* tfz_dat, const int* tfz_lo, const int* tfz_hi,
                            amrex::Real*       wmac_dat, const int* wmac_lo, const int* wmac_hi,
														int* corner_couple,
#endif
                            const amrex::Real* dt, const amrex::Real* dx,
                            int* use_forces_in_trans, int* ppm_type);

    void extrap_state_to_faces(const int* lo, const int* hi,
                               const amrex::Real* s_dat, const int* s_lo, const int* s_hi,  const int* nc,
                               const amrex::Real* tf_dat, const int* tf_lo, const int* tf_hi,
                               const amrex::Real* du_dat, const int* du_lo, const int* du_hi,
                               const amrex::Real* u_dat, const int* u_lo, const int* u_hi,
                               amrex::Real* xstate_dat, const int* xstate_lo, const int* xstate_hi,
                               const amrex::Real* v_dat, const int* v_lo, const int* v_hi,
                               amrex::Real* ystate_dat, const int* ystate_lo, const int* ystate_hi,
#if (AMREX_SPACEDIM == 3)
                               const amrex::Real* w_dat, const int* w_lo, const int* w_hi,
                               amrex::Real* zstate_dat, const int* zstate_lo, const int* zstate_hi,
															 int* corner_couple,
#endif
                               const amrex::Real* dt, const amrex::Real* dx, const int* bc, const int* state_ind,
                               const int* use_forces_in_trans, const int* ppm_type, const int* iconserv);


   void adv_forcing(const amrex::Real* aofs_dat, ARLIM_P(a_lo), ARLIM_P(a_hi),

			 const amrex::Real* xflux_dat, ARLIM_P(xflux_lo),ARLIM_P(xflux_hi),
                         const amrex::Real* uedge_dat, ARLIM_P(uedge_lo),ARLIM_P(uedge_hi),
			 const amrex::Real* areax_dat, ARLIM_P(ax_lo), ARLIM_P(ax_hi),

			 const amrex::Real* yflux_dat, ARLIM_P(yflux_lo),ARLIM_P(yflux_hi),
                         const amrex::Real* vedge_dat, ARLIM_P(vedge_lo),ARLIM_P(vedge_hi),
			 const amrex::Real* areay_dat, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
#if (BL_SPACEDIM == 3)
			 const amrex::Real* zflux_dat, ARLIM_P(zflux_lo),ARLIM_P(zflux_hi),
                         const amrex::Real* wedge_dat, ARLIM_P(wedge_lo),ARLIM_P(wedge_hi),
			 const amrex::Real* areaz_dat, ARLIM_P(az_lo), ARLIM_P(az_hi),
#endif
			 const amrex::Real* vol_dat, ARLIM_P(v_lo),ARLIM_P(v_hi),
			 const int* lo, const int* hi, int* iconserv);

   void sync_adv_forcing(const amrex::Real* sync_dat,ARLIM_P(s_lo),ARLIM_P(s_hi),

			      const amrex::Real* xflux_dat, ARLIM_P(xflux_lo),ARLIM_P(xflux_hi),
                              const amrex::Real* ucorr_dat, ARLIM_P(ucorr_lo),ARLIM_P(ucorr_hi),
			      const amrex::Real* areax_dat, ARLIM_P(ax_lo), ARLIM_P(ax_hi),

			      const amrex::Real* yflux_dat, ARLIM_P(yflux_lo),ARLIM_P(yflux_hi),
                              const amrex::Real* vcorr_dat, ARLIM_P(vcorr_lo),ARLIM_P(vcorr_hi),
			      const amrex::Real* areay_dat, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
#if (BL_SPACEDIM == 3)
			      const amrex::Real* zflux_dat, ARLIM_P(zflux_lo),ARLIM_P(zflux_hi),
                              const amrex::Real* wcorr_dat, ARLIM_P(wcorr_lo),ARLIM_P(wcorr_hi),
			      const amrex::Real* areaz_dat, ARLIM_P(az_lo), ARLIM_P(az_hi),
#endif
			      const amrex::Real* vol_dat, ARLIM_P(v_lo),ARLIM_P(v_hi),
			      const int* lo, const int* hi);

   void convscalminmax (const amrex::Real* Sold_dat,
                             ARLIM_P(slo), ARLIM_P(shi),
                             const amrex::Real* Snew_dat,
                             ARLIM_P(snlo), ARLIM_P(snhi),
#if (BL_SPACEDIM == 3)
			     const amrex::Real* smin_dat, const amrex::Real* smax_dat,
			     ARLIM_P(lo),ARLIM_P(hi),
#endif
			     const int* lo, const int* hi, const int* bc);

   void consscalminmax (const amrex::Real* Sold_dat,
			     const amrex::Real*  Rho_dat,
                             ARLIM_P(solo), ARLIM_P(sohi),
			     const amrex::Real* Snew_dat,
			     const amrex::Real* Rhon_dat,
			     ARLIM_P(snlo), ARLIM_P(snhi),
#if (BL_SPACEDIM == 3)
			     const amrex::Real* smin_dat, const amrex::Real* smax_dat,
			     ARLIM_P(lo),ARLIM_P(hi),
#endif
			     const int* lo, const int* hi, const int* bc);

   void  update_tf(const amrex::Real* SOdat, ARLIM_P(slo), ARLIM_P(shi),
			const amrex::Real* SNdat, ARLIM_P(slo), ARLIM_P(shi),
			const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			const int* lo, const int* hi, amrex::Real* dt, int* num_comp);

   void fort_correct_tf(const amrex::Real* SSdat, ARLIM_P(slo), ARLIM_P(shi),
			const amrex::Real* SPdat, ARLIM_P(splo), ARLIM_P(sphi),
			const amrex::Real* TSdat, ARLIM_P(tlo), ARLIM_P(thi),
			const amrex::Real* TNdat, ARLIM_P(tnlo), ARLIM_P(tnhi),
			const int* lo, const int* hi, amrex::Real* dt, int* num_comp);

   void update_aofs_tf(const amrex::Real* SOdat, ARLIM_P(slo), ARLIM_P(shi),
			    const amrex::Real* SNdat, ARLIM_P(slo), ARLIM_P(shi),
			    const amrex::Real* AOdat, ARLIM_P(alo), ARLIM_P(ahi),
			    const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			    const int* lo, const int* hi,
			    amrex::Real* dt, int* num_comp);

   void update_aofs_tf_gp(const amrex::Real* UOdat, ARLIM_P(ulo), ARLIM_P(uhi),
			       const amrex::Real* UNdat, ARLIM_P(ulo), ARLIM_P(uhi),
			       const amrex::Real* AOdat, ARLIM_P(alo), ARLIM_P(ahi),
			       const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			       const amrex::Real* GPdat, ARLIM_P(glo), ARLIM_P(ghi),
			       const amrex::Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
			       const int* lo, const int* hi, amrex::Real* dt);

   void fort_sum_tf_gp_visc(const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			    const amrex::Real* VIdat, ARLIM_P(vlo), ARLIM_P(vhi),
			    const amrex::Real* GPdat, ARLIM_P(glo), ARLIM_P(ghi),
			    const amrex::Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
			    const int* tlo, const int* thi);

   void fort_sum_tf_divu_visc(const amrex::Real* Sdat,  ARLIM_P(slo), ARLIM_P(shi),
			      const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			      const amrex::Real* DUdat, ARLIM_P(dlo), ARLIM_P(dhi),
			      const amrex::Real* VIdat, ARLIM_P(vlo), ARLIM_P(vhi),
			      const amrex::Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
			      const int* tlo, const int* thi,
			      int* num_comp, int* iconserv);

   void fort_sum_tf_gp(const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
		       const amrex::Real* GPdat, ARLIM_P(glo), ARLIM_P(ghi),
		       const amrex::Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
		       const int* tlo, const int* thi);

   void fort_sum_tf_divu(const amrex::Real* Sdat,  ARLIM_P(slo), ARLIM_P(shi),
			 const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			 const amrex::Real* DUdat, ARLIM_P(dlo), ARLIM_P(dhi),
			 const amrex::Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
			 const int* tlo, const int* thi,
			 int* num_comp, int* iconserv);

   void set_params(const int& slope_order,
			const int& use_unlimited_slopes);
#ifdef __cplusplus
}
#endif

#endif
#endif /*_GODUNOV_F_H_*/
