#ifndef _GODUNOV_F_H_
#define _GODUNOV_F_H_

/*
** $Id: GODUNOV_F.H,v 1.7 1999-08-25 18:04:04 propp Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_TEST_U_RHO         testurho
#    define FORT_TEST_UMAC_RHO      testumacrho
#    define FORT_ESTDT              estdt
#    define FORT_TRANSVEL           transvel
#    define FORT_ESTATE             estate
#    define FORT_ADV_FORCING        advforcing
#    define FORT_SYNC_ADV_FORCING   syncadvforcing
#    define FORT_SCALMINMAX         scalminmax
#    define FORT_UPDATE_TF          updatetf
#    define FORT_CORRECT_TF         correcttf
#    define FORT_UPDATE_AOFS_TF     updateaofstf
#    define FORT_UPDATE_AOFS_TF_GP  updateaofstfgp
#    define FORT_SUM_TF_GP_VISC     sumtfgpvisc
#    define FORT_SUM_TF_DIVU_VISC   sumtfdivuvisc
#    define FORT_SUM_TF_GP          sumtfgp
#    define FORT_SUM_TF_DIVU        sumtfdivu
#    define FORT_SET_PARAMS         setparams
#endif

#ifdef BL_LANG_CC
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_TEST_U_RHO         TESTURHO 
#    define FORT_TEST_UMAC_RHO      TESTUMACRHO
#    define FORT_ESTDT              ESTDT 
#    define FORT_TRANSVEL           TRANSVEL  
#    define FORT_ESTATE             ESTATE    
#    define FORT_ADV_FORCING        ADVFORCING    
#    define FORT_SYNC_ADV_FORCING   SYNCADVFORCING   
#    define FORT_SCALMINMAX         SCALMINMAX
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
#    define FORT_ESTDT              estdt
#    define FORT_TRANSVEL           transvel
#    define FORT_ESTATE             estate
#    define FORT_ADV_FORCING        advforcing
#    define FORT_SYNC_ADV_FORCING   syncadvforcing
#    define FORT_SCALMINMAX         scalminmax
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
#    define FORT_ESTDT              estdt_
#    define FORT_TRANSVEL           transvel_
#    define FORT_ESTATE             estate_
#    define FORT_ADV_FORCING        advforcing_
#    define FORT_SYNC_ADV_FORCING   syncadvforcing_
#    define FORT_SCALMINMAX         scalminmax_
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

#include <ArrayLim.H>

extern "C"
{
   void FORT_TEST_U_RHO(const Real* u,  ARLIM_P(vlo), ARLIM_P(vhi),
			const Real* v,  ARLIM_P(vlo), ARLIM_P(vhi),
#if (BL_SPACEDIM == 3)                          
			const Real* w,  ARLIM_P(vlo), ARLIM_P(vhi),
#endif
			const Real* rh, ARLIM_P(rlo), ARLIM_P(rhi),
			const int* lo, const int* hi, const Real* dt, 
			const Real* dx, Real* cflmax, 
			const Real* u_max, int* verbose);
   
   void FORT_TEST_UMAC_RHO(const Real* um, ARLIM_P(ulo), ARLIM_P(uhi),
			   const Real* vm, ARLIM_P(vlo), ARLIM_P(vhi),
#if (BL_SPACEDIM == 3)                            
			   const Real* wm, ARLIM_P(wlo), ARLIM_P(whi),
#endif                                              
			   const Real* rh, ARLIM_P(rlo), ARLIM_P(rhi),
			   const int* lo, const int* hi, const Real* dt, 
			   const Real* dx, Real* cfl, const Real* u_max);
   
   void FORT_ESTDT(const Real* Udat,  ARLIM_P(vlo), ARLIM_P(vhi),
		   const Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
		   const Real* rdat,  ARLIM_P(rlo), ARLIM_P(rhi),
		   const int* lo, const int* hi, Real* dt, 
		   const Real* dx, Real* cfl, Real* u_max);
   
   void FORT_TRANSVEL( const Real* u_dat, const Real* uad_dat, 
		       const Real* xhi_dat, const Real* slx_dat, 
		       const int* ubc, Real* slxscr, 
		       const Real* v_dat, const Real* vad_dat, 
		       const Real* yhi_dat, const Real* sly_dat, 
		       const int* vbc, Real*slyscr, 
#if (BL_SPACEDIM == 3)
		       const Real* w_dat, const Real* wad_dat, 
		       const Real* zhi_dat, const Real* slz_dat, 
		       const int* wbc, Real* slzscr,
#endif    
		       ARLIM_P(u_lo),ARLIM_P(u_hi),
		       ARLIM_P(w_lo),ARLIM_P(w_hi),
		       const int* lo, const int* hi, Real* dt, 
		       const Real* dx, int* use_forces_in_trans, 
		       const Real* tforcedat);
   
   void FORT_ESTATE(const Real* s_dat, const Real* tfr_dat, 
		    ARLIM_P(s_lo),ARLIM_P(s_hi),
		    const Real* u_dat, const Real* xlo_dat, 
		    const Real* xhi_dat, const Real* slx_dat, 
		    const Real* uad_dat,
		    Real* slxscr, Real* stxlo, Real* stxhi,
		    const Real* uedge_dat, const Real* stx_dat, 
		    ARLIM_P(u_lo),ARLIM_P(u_hi), 
		    
		    const Real* v_dat, const Real* ylo_dat, 
		    const Real* yhi_dat, const Real* sly_dat, 
		    const Real* vad_dat,
		    Real* slyscr, Real* stylo, Real* styhi,
		    const Real* vedge_dat, const Real* sty_dat, 
		    ARLIM_P(v_lo),ARLIM_P(v_hi),
#if (BL_SPACEDIM == 3)
		    const Real* w_dat, const Real* zlo_dat, 
		    const Real* zhi_dat, const Real* slz_dat, 
		    const Real* wad_dat,
		    Real* slzscr, Real* stzlo, Real* stzhi,
		    const Real* wedge_dat, const Real* stz_dat,
		    ARLIM_P(w_lo),ARLIM_P(w_hi), 
#endif
		    ARLIM_P(ww_lo),ARLIM_P(ww_hi),
		    const int* bc, const int* lo, const int* hi, 
		    Real* dt, const Real* dx, int* fort_ind, int* velpred, 
		    int* use_forces_in_trans);
   
   void FORT_ADV_FORCING(const Real* aofs_dat, ARLIM_P(a_lo), ARLIM_P(a_hi),
			 const Real* xflux_dat, const Real* uedge_dat, 
			 ARLIM_P(xflux_lo),ARLIM_P(xflux_hi),
			 const Real* areax_dat, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
			 const Real* yflux_dat, const Real* vedge_dat, 
			 ARLIM_P(yflux_lo),ARLIM_P(yflux_hi),
			 const Real* areay_dat, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
#if (BL_SPACEDIM == 3)                                                    
			 const Real* zflux_dat, const Real* wedge_dat, 
			 ARLIM_P(zflux_lo),ARLIM_P(zflux_hi),
			 const Real* areaz_dat, ARLIM_P(az_lo), ARLIM_P(az_hi),
#endif
			 const Real* vol_dat, ARLIM_P(v_lo),ARLIM_P(v_hi),
			 const int* lo, const int* hi, int* iconserv);
   
   void FORT_SYNC_ADV_FORCING(const Real* sync_dat,ARLIM_P(s_lo),ARLIM_P(s_hi),
			      
			      const Real* xflux_dat, const Real* ucorr_dat, 
			      ARLIM_P(xflux_lo),ARLIM_P(xflux_hi),
			      const Real* areax_dat, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
			      
			      const Real* yflux_dat, const Real* vcorr_dat, 
			      ARLIM_P(yflux_lo),ARLIM_P(yflux_hi),
			      const Real* areay_dat, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
#if (BL_SPACEDIM == 3)                                             
			      const Real* zflux_dat, const Real* wcorr_dat, 
			      ARLIM_P(zflux_lo),ARLIM_P(zflux_hi),
			      const Real* areaz_dat, ARLIM_P(az_lo), ARLIM_P(az_hi),
#endif
			      const Real* vol_dat, ARLIM_P(v_lo),ARLIM_P(v_hi),
			      const int* lo, const int* hi, int* iconserv);
   
   void FORT_SCALMINMAX (const Real* Sold_dat, ARLIM_P(slo), ARLIM_P(shi),
			 const Real* Snew_dat, ARLIM_P(slo), ARLIM_P(shi),
#if (BL_SPACEDIM == 3)
			 const Real* smin_dat, const Real* smax_dat,
			 ARLIM_P(lo),ARLIM_P(hi),
#endif
			 const int* lo, const int* hi, const int* bc);
   
   void  FORT_UPDATE_TF(const Real* SOdat, ARLIM_P(slo), ARLIM_P(shi), 
			const Real* SNdat, ARLIM_P(slo), ARLIM_P(shi),
			const Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			const int* lo, const int* hi, Real* dt, int* num_comp);
   
   void FORT_CORRECT_TF(const Real* SSdat, const Real* SPdat, ARLIM_P(slo), ARLIM_P(shi),
			const Real* TSdat, const Real* TNdat, ARLIM_P(tlo), ARLIM_P(thi),
			const int* lo, const int* hi, Real* dt, int* num_comp);
   
   void FORT_UPDATE_AOFS_TF(const Real* SOdat, ARLIM_P(slo), ARLIM_P(shi), 
			    const Real* SNdat, ARLIM_P(slo), ARLIM_P(shi),
			    const Real* AOdat, ARLIM_P(alo), ARLIM_P(ahi),
			    const Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			    const int* lo, const int* hi, 
			    Real* dt, int* num_comp);
   
   void FORT_UPDATE_AOFS_TF_GP(const Real* UOdat, ARLIM_P(ulo), ARLIM_P(uhi),
			       const Real* UNdat, ARLIM_P(ulo), ARLIM_P(uhi),
			       const Real* AOdat, ARLIM_P(alo), ARLIM_P(ahi),
			       const Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			       const Real* GPdat, ARLIM_P(glo), ARLIM_P(ghi),
			       const Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
			       const int* lo, const int* hi, Real* dt);
   
   void FORT_SUM_TF_GP_VISC(const Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			    const Real* VIdat, ARLIM_P(vlo), ARLIM_P(vhi),
			    const Real* GPdat, ARLIM_P(glo), ARLIM_P(ghi),
			    const Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
			    const int* tlo, const int* thi);
   
   void FORT_SUM_TF_DIVU_VISC(const Real* Sdat,  ARLIM_P(slo), ARLIM_P(shi),
			      const Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			      const Real* DUdat, ARLIM_P(dlo), ARLIM_P(dhi),
			      const Real* VIdat, ARLIM_P(vlo), ARLIM_P(vhi),
			      const Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
			      const int* tlo, const int* thi, 
			      int* num_comp, int* iconserv);
   
   void FORT_SUM_TF_GP(const Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
		       const Real* GPdat, ARLIM_P(glo), ARLIM_P(ghi),
		       const Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
		       const int* tlo, const int* thi);
   
   void FORT_SUM_TF_DIVU(const Real* Sdat,  ARLIM_P(slo), ARLIM_P(shi),
			 const Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			 const Real* DUdat, ARLIM_P(dlo), ARLIM_P(dhi),
			 const Real* RHdat, ARLIM_P(rlo), ARLIM_P(rhi),
			 const int* tlo, const int* thi, 
			 int* num_comp, int* iconserv);
   
   void FORT_SET_PARAMS(const int& slope_order,
			const int& use_unlimited_slopes);
}
#endif

#endif /*_GODUNOV_F_H_*/
