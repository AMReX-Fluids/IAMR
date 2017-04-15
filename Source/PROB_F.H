
#ifndef _PROB_NS_F_H_
#define _PROB_NS_F_H_

#ifdef BL_LANG_FORT
#    define FORT_INITDATA    initdata
#    define FORT_DENFILL     denfill
#    define FORT_ADVFILL     advfill
#    define FORT_ADV2FILL    adv2fill
#    define FORT_XVELFILL    xvelfill
#    define FORT_YVELFILL    yvelfill
#    define FORT_ZVELFILL    zvelfill
#    define FORT_PRESFILL    presfill
#    define FORT_TEMPFILL    tempfill
#    define FORT_DIVUFILL    divufill
#    define FORT_DSDTFILL    dsdtfill
#if defined(GENGETFORCE) || defined(MOREGENGETFORCE)
#    define FORT_AVERAGE_EDGE_STATES      average_edge_states
#    define FORT_MAKEFORCE      makeforce
#endif
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_INITDATA    INITDATA
#    define FORT_DENFILL     DENFILL
#    define FORT_ADVFILL     ADVFILL
#    define FORT_ADV2FILL    ADV2FILL
#    define FORT_XVELFILL    XVELFILL
#    define FORT_YVELFILL    YVELFILL
#    define FORT_ZVELFILL    ZVELFILL
#    define FORT_PRESFILL    PRESFILL
#    define FORT_TEMPFILL    TEMPFILL
#    define FORT_DIVUFILL    DIVUFILL
#    define FORT_DSDTFILL    DSDTFILL
#if defined(GENGETFORCE) || defined(MOREGENGETFORCE)
#    define FORT_AVERAGE_EDGE_STATES      AVERAGE_EDGE_STATES
#    define FORT_MAKEFORCE      MAKEFORCE
#endif
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_INITDATA    initdata
#    define FORT_DENFILL     denfill
#    define FORT_ADVFILL     advfill
#    define FORT_ADV2FILL    adv2fill
#    define FORT_XVELFILL    xvelfill
#    define FORT_YVELFILL    yvelfill
#    define FORT_ZVELFILL    zvelfill
#    define FORT_PRESFILL    presfill
#    define FORT_TEMPFILL    tempfill
#    define FORT_DIVUFILL    divufill
#    define FORT_DSDTFILL    dsdtfill
#if defined(GENGETFORCE) || defined(MOREGENGETFORCE)
#    define FORT_AVERAGE_EDGE_STATES      average_edge_states
#    define FORT_MAKEFORCE      makeforce
#endif
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_INITDATA    initdata_
#    define FORT_DENFILL     denfill_
#    define FORT_ADVFILL     advfill_
#    define FORT_ADV2FILL    adv2fill_
#    define FORT_XVELFILL    xvelfill_
#    define FORT_YVELFILL    yvelfill_
#    define FORT_ZVELFILL    zvelfill_
#    define FORT_PRESFILL    presfill_
#    define FORT_TEMPFILL    tempfill_
#    define FORT_DIVUFILL    divufill_
#    define FORT_DSDTFILL    dsdtfill_
#if defined(GENGETFORCE) || defined(MOREGENGETFORCE)
#    define FORT_AVERAGE_EDGE_STATES      average_edge_states_
#    define FORT_MAKEFORCE      makeforce_
#endif
#  endif

#include <AMReX_ArrayLim.H>

  extern "C" {

     void FORT_INITDATA (const int* level, const amrex::Real* time,
                         const int* lo, const int* hi,
                         const int* nscal,
                         amrex::Real* vel, amrex::Real* scal, 
                         ARLIM_P(s_lo), ARLIM_P(s_hi),
                         amrex::Real* press,
                         ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const amrex::Real* dx,
                         const amrex::Real* xlo, const amrex::Real* xhi);

     void FORT_DENFILL  (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_ADVFILL  (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi), 
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_ADV2FILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi), 
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_XVELFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_YVELFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

#if (BL_SPACEDIM == 3)
     void FORT_ZVELFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);
#endif
     void FORT_PRESFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_TEMPFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi), 
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_DIVUFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_DSDTFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

#ifdef BOUSSINESQ
     void FORT_MAKEFORCE (amrex::Real* force,
			  const amrex::Real* scal,
                          ARLIM_P(f_lo), ARLIM_P(f_hi),
                          ARLIM_P(s_lo), ARLIM_P(s_hi),
                          const amrex::Real* dx,
                          const amrex::Real* xlo, const amrex::Real* xhi,
			  const int* scomp, const int* ncomp, const int* nscal);
#else
#ifdef MOREGENGETFORCE

     void FORT_AVERAGE_EDGE_STATES(const amrex::Real* Vel,
				   const amrex::Real* umacx,
				   const amrex::Real* umacy,
#if (BL_SPACEDIM==3)
				   const amrex::Real* umacz,
#endif
				   ARLIM_P(vel_lo),  ARLIM_P(vel_hi),
				   ARLIM_P(umacx_lo), ARLIM_P(umacx_hi),
				   ARLIM_P(umacy_lo), ARLIM_P(umacy_hi),
#if (BL_SPACEDIM==3)
				   ARLIM_P(umacz_lo), ARLIM_P(umacz_hi),
#endif
				   const int* getForceVerbose);

     void FORT_MAKEFORCE (const amrex::Real* time,
			  amrex::Real* force,
			  const amrex::Real* vel,
			  const amrex::Real* scal,
                          ARLIM_P(f_lo), ARLIM_P(f_hi),
                          ARLIM_P(v_lo), ARLIM_P(v_hi),
                          ARLIM_P(s_lo), ARLIM_P(s_hi),
                          const amrex::Real* dx,
                          const amrex::Real* xlo, const amrex::Real* xhi,
			  const amrex::Real* gravity,
			  const int* scomp, const int* ncomp, const int* nscal,
			  const int* getForceVerbose);
#endif
#ifdef GENGETFORCE
     void FORT_MAKEFORCE (const amrex::Real* time, amrex::Real* force, const amrex::Real* rho,
                          ARLIM_P(i_lo), ARLIM_P(i_hi),
                          ARLIM_P(s_lo), ARLIM_P(s_hi),
                          const amrex::Real* dx,
                          const amrex::Real* xlo, const amrex::Real* xhi,
			  const amrex::Real* gravity,
			  const int* scomp, const int* ncomp);
#endif
#endif
  }
#endif
#endif /*_PROB_F_H_*/
