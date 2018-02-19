
#ifndef _GODUNOV_F_H_
#define _GODUNOV_F_H_

#ifdef BL_LANG_FORT
#    define FORT_TEST_U_RHO         testurho
#    define FORT_TEST_UMAC_RHO      testumacrho
#    define FORT_ESTDT              estmdt
#    define FORT_MAXCHNG_VELMAG     maxchngvelmag
#    define FORT_TRANSVEL           transvel
#    define FORT_ESTATE             estate
#    define FORT_ESTATE_FPU         estatefpu
#    define FORT_ESTATE_BDS         estatebds
#    define FORT_ADV_FORCING        advforcing
#    define FORT_SYNC_ADV_FORCING   syncadvforcing
#    define FORT_CONSSCALMINMAX     consscalminmax
#    define FORT_CONVSCALMINMAX     convscalminmax
#    define FORT_UPDATE_TF          updatetf
#    define FORT_CORRECT_TF         correcttf
#    define FORT_UPDATE_AOFS_TF     updateaofstf
#    define FORT_UPDATE_AOFS_TF_GP  updateaofstfgp
#    define FORT_SUM_TF_GP_VISC     sumtfgpvisc
#    define FORT_SUM_TF_DIVU_VISC   sumtfdivuvisc
#    define FORT_SUM_TF_GP          sumtfgp
#    define FORT_SUM_TF_DIVU        sumtfdivu
#    define FORT_SET_PARAMS         setparams
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_TEST_U_RHO         TESTURHO 
#    define FORT_TEST_UMAC_RHO      TESTUMACRHO
#    define FORT_ESTDT              ESTMDT 
#    define FORT_MAXCHNG_VELMAG     MAXCHNGVELMAG
#    define FORT_TRANSVEL           TRANSVEL  
#    define FORT_ESTATE             ESTATE    
#    define FORT_ESTATE_FPU         ESTATEFPU
#    define FORT_ESTATE_BDS         ESTATEBDS
#    define FORT_ADV_FORCING        ADVFORCING    
#    define FORT_SYNC_ADV_FORCING   SYNCADVFORCING   
#    define FORT_CONSSCALMINMAX     CONSSCALMINMAX
#    define FORT_CONVSCALMINMAX     CONVSCALMINMAX
#    define FORT_UPDATE_TF          UPDATETF 
#    define FORT_CORRECT_TF         CORRECTTF 
#    define FORT_UPDATE_AOFS_TF     UPDATEAOFSTF 
#    define FORT_UPDATE_AOFS_TF_GP  UPDATEAOFSTFGP 
#    define FORT_SUM_TF_GP_VISC     SUMTFGPVISC 
#    define FORT_SUM_TF_DIVU_VISC   SUMTFDIVUVISC 
#    define FORT_SUM_TF_GP          SUMTFGP 
#    define FORT_SUM_TF_DIVU        SUMTFDIVU 
#    define FORT_SET_PARAMS         SETPARAMS
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_TEST_U_RHO         testurho
#    define FORT_TEST_UMAC_RHO      testumacrho
#    define FORT_ESTDT              estmdt
#    define FORT_MAXCHNG_VELMAG     maxchngvelmax
#    define FORT_TRANSVEL           transvel
#    define FORT_ESTATE             estate
#    define FORT_ESTATE_FPU         estatefpu
#    define FORT_ESTATE_BDS         estatebds
#    define FORT_ADV_FORCING        advforcing
#    define FORT_SYNC_ADV_FORCING   syncadvforcing
#    define FORT_CONSSCALMINMAX     consscalminmax
#    define FORT_CONVSCALMINMAX     convscalminmax
#    define FORT_UPDATE_TF          updatetf
#    define FORT_CORRECT_TF         correcttf
#    define FORT_UPDATE_AOFS_TF     updateaofstf
#    define FORT_UPDATE_AOFS_TF_GP  updateaofstfgp
#    define FORT_SUM_TF_GP_VISC     sumtfgpvisc
#    define FORT_SUM_TF_DIVU_VISC   sumtfdivuvisc
#    define FORT_SUM_TF_GP          sumtfgp
#    define FORT_SUM_TF_DIVU        sumtfdivu
#    define FORT_SET_PARAMS         setparams
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_TEST_U_RHO         testurho_
#    define FORT_TEST_UMAC_RHO      testumacrho_
#    define FORT_ESTDT              estmdt_
#    define FORT_MAXCHNG_VELMAG     maxchngvelmag_
#    define FORT_TRANSVEL           transvel_
#    define FORT_ESTATE             estate_
#    define FORT_ESTATE_FPU         estatefpu_
#    define FORT_ESTATE_BDS         estatebds_
#    define FORT_ADV_FORCING        advforcing_
#    define FORT_SYNC_ADV_FORCING   syncadvforcing_
#    define FORT_CONSSCALMINMAX     consscalminmax_
#    define FORT_CONVSCALMINMAX     convscalminmax_
#    define FORT_UPDATE_TF          updatetf_
#    define FORT_CORRECT_TF         correcttf_
#    define FORT_UPDATE_AOFS_TF     updateaofstf_
#    define FORT_UPDATE_AOFS_TF_GP  updateaofstfgp_
#    define FORT_SUM_TF_GP_VISC     sumtfgpvisc_
#    define FORT_SUM_TF_DIVU_VISC   sumtfdivuvisc_
#    define FORT_SUM_TF_GP          sumtfgp_
#    define FORT_SUM_TF_DIVU        sumtfdivu_
#    define FORT_SET_PARAMS         setparams_
#  endif

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C"
{
#endif
   void FORT_TEST_U_RHO(const amrex::Real* u,  ARLIM_P(vlo), ARLIM_P(vhi),
			const amrex::Real* v,  ARLIM_P(vlo), ARLIM_P(vhi),
#if (BL_SPACEDIM == 3)                          
			const amrex::Real* w,  ARLIM_P(vlo), ARLIM_P(vhi),
#endif
			const amrex::Real* rh, ARLIM_P(rlo), ARLIM_P(rhi),
			const int* lo, const int* hi, const amrex::Real* dt, 
			const amrex::Real* dx, amrex::Real* cflmax, 
			const amrex::Real* u_max, int* verbose);
   
   void FORT_TEST_UMAC_RHO(const amrex::Real* um, ARLIM_P(ulo), ARLIM_P(uhi),
			   const amrex::Real* vm, ARLIM_P(vlo), ARLIM_P(vhi),
#if (BL_SPACEDIM == 3)                            
			   const amrex::Real* wm, ARLIM_P(wlo), ARLIM_P(whi),
#endif                                              
			   const amrex::Real* rh, ARLIM_P(rlo), ARLIM_P(rhi),
			   const int* lo, const int* hi, const amrex::Real* dt, 
			   const amrex::Real* dx, amrex::Real* cfl, const amrex::Real* u_max);
   
   void FORT_ESTDT(const amrex::Real* Udat,  ARLIM_P(vlo), ARLIM_P(vhi),
		   const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
		   const amrex::Real* rdat,  ARLIM_P(rlo), ARLIM_P(rhi),
		   const int* lo, const int* hi, amrex::Real* dt, 
		   const amrex::Real* dx, amrex::Real* cfl, amrex::Real* u_max);
   
   void FORT_MAXCHNG_VELMAG(const amrex::Real* Uodat,  ARLIM_P(uo_lo), ARLIM_P(uo_hi),
		   const amrex::Real* Undat,  ARLIM_P(un_lo), ARLIM_P(un_hi),
		   const int* lo, const int* hi, amrex::Real* max_change);
   
   void FORT_TRANSVEL( const amrex::Real* u_dat, const amrex::Real* uad_dat, 
		       const amrex::Real* xhi_dat, const amrex::Real* slx_dat, 
		       const int* ubc, amrex::Real* slxscr, amrex::Real* Imx, amrex::Real* Ipx, 
                       amrex::Real* sedgex, ARLIM_P(sedgexlo), ARLIM_P(sedgexhi),

		       const amrex::Real* v_dat, const amrex::Real* vad_dat, 
		       const amrex::Real* yhi_dat, const amrex::Real* sly_dat, 
		       const int* vbc, amrex::Real*slyscr, amrex::Real* Imy, amrex::Real* Ipy,
                       amrex::Real* sedgey, ARLIM_P(sedgeylo), ARLIM_P(sedgeyhi),
#if (BL_SPACEDIM == 3)
		       const amrex::Real* w_dat, const amrex::Real* wad_dat, 
		       const amrex::Real* zhi_dat, const amrex::Real* slz_dat, 
		       const int* wbc, amrex::Real* slzscr, amrex::Real* Imz, amrex::Real* Ipz,
                       amrex::Real* sedgez, ARLIM_P(sedgezlo), ARLIM_P(sedgezhi),
#endif    
		       ARLIM_P(u_lo),ARLIM_P(u_hi),
		       ARLIM_P(w_lo),ARLIM_P(w_hi),
		       ARLIM_P(I_lo),ARLIM_P(I_hi),
		       amrex::Real* dsvl, ARLIM_P(dsvl_lo),ARLIM_P(dsvl_hi),
		       amrex::Real* sm, amrex::Real* sp, ARLIM_P(smp_lo),ARLIM_P(smp_hi),
		       const int* lo, const int* hi, amrex::Real* dt, 
		       const amrex::Real* dx, int* use_forces_in_trans, 
		       const amrex::Real* tforcedat, int* ppm_type);
   
   void FORT_ESTATE(const amrex::Real* s_dat, const amrex::Real* tfr_dat, 
		    ARLIM_P(s_lo), ARLIM_P(s_hi),
		    const amrex::Real* u_dat, ARLIM_P(u_lo),ARLIM_P(u_hi),
		    const amrex::Real* xlo_dat, const amrex::Real* xhi_dat,
		    const amrex::Real* slx_dat, const amrex::Real* uad_dat,
		    amrex::Real* slxscr, amrex::Real* stxlo, amrex::Real* stxhi,
		    const amrex::Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
                    const amrex::Real* stx,       ARLIM_P(stx_lo),ARLIM_P(stx_hi), amrex::Real* Imx, amrex::Real* Ipx,
                    amrex::Real* sedgex, ARLIM_P(sedgexlo), ARLIM_P(sedgexhi),
		    
		    const amrex::Real* v_dat, const amrex::Real* ylo_dat, 
		    const amrex::Real* yhi_dat, const amrex::Real* sly_dat, 
		    const amrex::Real* vad_dat,
		    amrex::Real* slyscr, amrex::Real* stylo, amrex::Real* styhi,
		    const amrex::Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
                    const amrex::Real* sty,       ARLIM_P(sty_lo),ARLIM_P(sty_hi), amrex::Real* Imy, amrex::Real* Ipy,
                    amrex::Real* sedgey, ARLIM_P(sedgeylo), ARLIM_P(sedgeyhi),
#if (BL_SPACEDIM == 3)
		    const amrex::Real* w_dat, const amrex::Real* zlo_dat, 
		    const amrex::Real* zhi_dat, const amrex::Real* slz_dat, 
		    const amrex::Real* wad_dat,
		    amrex::Real* slzscr, amrex::Real* stzlo, amrex::Real* stzhi,
		    const amrex::Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
                    const amrex::Real* stz,       ARLIM_P(stz_lo),ARLIM_P(stz_hi), amrex::Real* Imz, amrex::Real* Ipz,
                    amrex::Real* sedgez, ARLIM_P(sedgezlo), ARLIM_P(sedgezhi),

		    const amrex::Real* xedge_dat, const amrex::Real* yedge_dat, const amrex::Real* zedge_dat,
		    const amrex::Real* xylo_dat,  const amrex::Real* xzlo_dat,  const amrex::Real* yxlo_dat,
		    const amrex::Real* yzlo_dat,  const amrex::Real* zxlo_dat,  const amrex::Real* zylo_dat,
		    const amrex::Real* xyhi_dat,  const amrex::Real* xzhi_dat,  const amrex::Real* yxhi_dat,
		    const amrex::Real* yzhi_dat,  const amrex::Real* zxhi_dat,  const amrex::Real* zyhi_dat,

		    int* corner_couple,
#endif
		    ARLIM_P(ww_lo),ARLIM_P(ww_hi),
		    ARLIM_P(I_lo),ARLIM_P(I_hi),
                    amrex::Real* dsvl, ARLIM_P(dsvl_lo),ARLIM_P(dsvl_hi),
                    amrex::Real* sm, amrex::Real* sp, ARLIM_P(smp_lo),ARLIM_P(smp_hi),
		    const int* bc, const int* lo, const int* hi, 
		    amrex::Real* dt, const amrex::Real* dx, int* fort_ind, int* velpred, 
		    int* use_forces_in_trans, int* ppm_type);
   
   void FORT_ESTATE_FPU(
                    const amrex::Real* s_dat,    ARLIM_P(s_lo),ARLIM_P(s_hi),
		    const amrex::Real* tfr_dat,  ARLIM_P(t_lo),ARLIM_P(t_hi),
		    const amrex::Real* divu_dat, ARLIM_P(d_lo),ARLIM_P(d_hi),
		    const amrex::Real* xlo_dat, 
		    const amrex::Real* xhi_dat, const amrex::Real* slx_dat, 
		    amrex::Real* slxscr, amrex::Real* stxlo, amrex::Real* stxhi,
		    const amrex::Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
                    const amrex::Real* stx,       ARLIM_P(stx_lo),ARLIM_P(stx_hi), amrex::Real* Imx, amrex::Real* Ipx,
                    amrex::Real* sedgex, ARLIM_P(sedgexlo), ARLIM_P(sedgexhi),

		    const amrex::Real* ylo_dat, 
		    const amrex::Real* yhi_dat, const amrex::Real* sly_dat, 
		    amrex::Real* slyscr, amrex::Real* stylo, amrex::Real* styhi,
		    const amrex::Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
                    const amrex::Real* sty,       ARLIM_P(sty_lo),ARLIM_P(sty_hi), amrex::Real* Imy, amrex::Real* Ipy,
                    amrex::Real* sedgey, ARLIM_P(sedgeylo), ARLIM_P(sedgeyhi),

#if (BL_SPACEDIM == 3)
		    const amrex::Real* zlo_dat, 
		    const amrex::Real* zhi_dat, const amrex::Real* slz_dat, 
		    amrex::Real* slzscr, amrex::Real* stzlo, amrex::Real* stzhi,
		    const amrex::Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
                    const amrex::Real* stz,       ARLIM_P(stz_lo),ARLIM_P(stz_hi), amrex::Real* Imz, amrex::Real* Ipz,
                    amrex::Real* sedgez, ARLIM_P(sedgezlo), ARLIM_P(sedgezhi),

		    const amrex::Real* xedge_dat, const amrex::Real* yedge_dat, const amrex::Real* zedge_dat,
		    const amrex::Real* xylo_dat,  const amrex::Real* xzlo_dat,  const amrex::Real* yxlo_dat,
		    const amrex::Real* yzlo_dat,  const amrex::Real* zxlo_dat,  const amrex::Real* zylo_dat,
		    const amrex::Real* xyhi_dat,  const amrex::Real* xzhi_dat,  const amrex::Real* yxhi_dat,
		    const amrex::Real* yzhi_dat,  const amrex::Real* zxhi_dat,  const amrex::Real* zyhi_dat,

       	            int* corner_couple,
#endif
		    ARLIM_P(ww_lo),ARLIM_P(ww_hi),
		    ARLIM_P(I_lo),ARLIM_P(I_hi),
                    amrex::Real* dsvl, ARLIM_P(dsvl_lo),ARLIM_P(dsvl_hi),
                    amrex::Real* sm, amrex::Real* sp, ARLIM_P(smp_lo),ARLIM_P(smp_hi),
		    const int* bc, const int* lo, const int* hi, 
		    amrex::Real* dt, const amrex::Real* dx, int* fort_ind, 
		    int* use_forces_in_trans, int* iconserv, int* ppm_type);

   void FORT_ESTATE_BDS(const amrex::Real* s_dat, const amrex::Real* tfr_dat, const amrex::Real* divu_dat,
		    ARLIM_P(s_lo),ARLIM_P(s_hi),
		    const amrex::Real* xlo_dat, 
		    const amrex::Real* xhi_dat, const amrex::Real* slx_dat, 
		    amrex::Real* slxscr, amrex::Real* stxlo, amrex::Real* stxhi,
		    const amrex::Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
                    const amrex::Real* stx,       ARLIM_P(stx_lo),ARLIM_P(stx_hi),
		    const amrex::Real* ylo_dat, 
		    const amrex::Real* yhi_dat, const amrex::Real* sly_dat, 
		    amrex::Real* slyscr, amrex::Real* stylo, amrex::Real* styhi,
		    const amrex::Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
                    const amrex::Real* sty,       ARLIM_P(sty_lo),ARLIM_P(sty_hi),
#if (BL_SPACEDIM == 3)
		    const amrex::Real* zlo_dat, 
		    const amrex::Real* zhi_dat, const amrex::Real* slz_dat, 
		    amrex::Real* slzscr, amrex::Real* stzlo, amrex::Real* stzhi,
		    const amrex::Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
                    const amrex::Real* stz,       ARLIM_P(stz_lo),ARLIM_P(stz_hi),
#endif
		    ARLIM_P(ww_lo),ARLIM_P(ww_hi),
		    const int* bc, const int* lo, const int* hi, 
		    amrex::Real* dt, const amrex::Real* dx, int* fort_ind, 
		    int* use_forces_in_trans, int* iconserv);
   
   void FORT_ADV_FORCING(const amrex::Real* aofs_dat, ARLIM_P(a_lo), ARLIM_P(a_hi),

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
   
   void FORT_SYNC_ADV_FORCING(const amrex::Real* sync_dat,ARLIM_P(s_lo),ARLIM_P(s_hi),
			      
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
   
   void FORT_CONVSCALMINMAX (const amrex::Real* Sold_dat, 
                             ARLIM_P(slo), ARLIM_P(shi),
                             const amrex::Real* Snew_dat, 
                             ARLIM_P(snlo), ARLIM_P(snhi),
#if (BL_SPACEDIM == 3)
			     const amrex::Real* smin_dat, const amrex::Real* smax_dat,
			     ARLIM_P(lo),ARLIM_P(hi),
#endif
			     const int* lo, const int* hi, const int* bc);

   void FORT_CONSSCALMINMAX (const amrex::Real* Sold_dat, const amrex::Real* Snew_dat, 
			     const amrex::Real*  Rho_dat, const amrex::Real* Rhon_dat,
                             ARLIM_P(slo), ARLIM_P(shi),
#if (BL_SPACEDIM == 3)
			     const amrex::Real* smin_dat, const amrex::Real* smax_dat,
			     ARLIM_P(lo),ARLIM_P(hi),
#endif
			     const int* lo, const int* hi, const int* bc);
   
   void  FORT_UPDATE_TF(const amrex::Real* SOdat, ARLIM_P(slo), ARLIM_P(shi), 
			const amrex::Real* SNdat, ARLIM_P(slo), ARLIM_P(shi),
			const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			const int* lo, const int* hi, amrex::Real* dt, int* num_comp);
   
   void FORT_CORRECT_TF(const amrex::Real* SSdat, const amrex::Real* SPdat, ARLIM_P(slo), ARLIM_P(shi),
			const amrex::Real* TSdat, const amrex::Real* TNdat, ARLIM_P(tlo), ARLIM_P(thi),
			const int* lo, const int* hi, amrex::Real* dt, int* num_comp);
   
   void FORT_UPDATE_AOFS_TF(const amrex::Real* SOdat, ARLIM_P(slo), ARLIM_P(shi), 
			    const amrex::Real* SNdat, ARLIM_P(slo), ARLIM_P(shi),
			    const amrex::Real* AOdat, ARLIM_P(alo), ARLIM_P(ahi),
			    const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			    const int* lo, const int* hi, 
			    amrex::Real* dt, int* num_comp);
   
   void FORT_UPDATE_AOFS_TF_GP(const amrex::Real* UOdat, ARLIM_P(ulo), ARLIM_P(uhi),
			       const amrex::Real* UNdat, ARLIM_P(ulo), ARLIM_P(uhi),
			       const amrex::Real* AOdat, ARLIM_P(alo), ARLIM_P(ahi),
			       const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			       const amrex::Real* GPdat, ARLIM_P(glo), ARLIM_P(ghi),
			       const amrex::Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
			       const int* lo, const int* hi, amrex::Real* dt);
   
   void FORT_SUM_TF_GP_VISC(const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			    const amrex::Real* VIdat, ARLIM_P(vlo), ARLIM_P(vhi),
			    const amrex::Real* GPdat, ARLIM_P(glo), ARLIM_P(ghi),
			    const amrex::Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
			    const int* tlo, const int* thi);
   
   void FORT_SUM_TF_DIVU_VISC(const amrex::Real* Sdat,  ARLIM_P(slo), ARLIM_P(shi),
			      const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			      const amrex::Real* DUdat, ARLIM_P(dlo), ARLIM_P(dhi),
			      const amrex::Real* VIdat, ARLIM_P(vlo), ARLIM_P(vhi),
			      const amrex::Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
			      const int* tlo, const int* thi, 
			      int* num_comp, int* iconserv);
   
   void FORT_SUM_TF_GP(const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
		       const amrex::Real* GPdat, ARLIM_P(glo), ARLIM_P(ghi),
		       const amrex::Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
		       const int* tlo, const int* thi);
   
   void FORT_SUM_TF_DIVU(const amrex::Real* Sdat,  ARLIM_P(slo), ARLIM_P(shi),
			 const amrex::Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			 const amrex::Real* DUdat, ARLIM_P(dlo), ARLIM_P(dhi),
			 const amrex::Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
			 const int* tlo, const int* thi, 
			 int* num_comp, int* iconserv);
   
   void FORT_SET_PARAMS(const int& slope_order,
			const int& use_unlimited_slopes);
#ifdef __cplusplus
}
#endif
#endif

#endif /*_GODUNOV_F_H_*/
