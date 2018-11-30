#ifndef _LS_F_H_
#define _LS_F_H_

#ifdef BL_LANG_FORT
#    define FORT_PHIUPD      phiupd
#    define FORT_LSCFL       lscfl
#    define FORT_NARROWBAND  narrowband
#    define FORT_RETYPIFY    retypify
#    define FORT_FINDINTRFCE findintrfce
#    define FORT_FASTMARCH   fastmarch
#    define FORT_NBANDNUMIFY nbandnumify
#    define FORT_MINE        setmine
#    define FORT_FASTMARCH2  fastmarch2

#if (BL_SPACEDIM == 3)
#    define FORT_SUMTURB     sumturb
#ifdef SUMJET
#    define FORT_SUMJET      sumjet
#endif
#endif
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_PHIUPD      PHIUPD
#    define FORT_INITPHI     INITPHI
#    define FORT_LSCFL       LSCFL
#    define FORT_NARROWBAND  NARROWBAND
#    define FORT_RETYPIFY    RETYPIFY
#    define FORT_FINDINTRFCE FINDINTRFCE
#    define FORT_FASTMARCH   FASTMARCH
#    define FORT_NBANDNUMIFY NBANDNUMIFY
#    define FORT_MINE        SETMINE
#    define FORT_FASTMARCH2  FASTMARCH2
#if (BL_SPACEDIM == 3)
#    define FORT_SUMTURB     SUMTURB
#ifdef SUMJET
#    define FORT_SUMJET      SUMJET
#endif
#endif
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_PHIUPD      phiupd
#    define FORT_INITPHI     initphi
#    define FORT_LSCFL       lscfl
#    define FORT_NARROWBAND  narrowband
#    define FORT_RETYPIFY    retypify
#    define FORT_FINDINTRFCE findintrfce
#    define FORT_FASTMARCH   fastmarch
#    define FORT_NBANDNUMIFY nbandnumify
#    define FORT_MINE        setmine
#    define FORT_FASTMARCH2  fastmarch2
#if (BL_SPACEDIM == 3)
#    define FORT_SUMTURB     sumturb
#ifdef SUMJET
#    define FORT_SUMJET      sumjet
#endif
#endif
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_PHIUPD      phiupd_
#    define FORT_INITPHI     initphi_
#    define FORT_LSCFL       lscfl_
#    define FORT_NARROWBAND  narrowband_
#    define FORT_RETYPIFY    retypify_
#    define FORT_FINDINTRFCE findintrfce_
#    define FORT_FASTMARCH   fastmarch_
#    define FORT_NBANDNUMIFY nbandnumify_
#    define FORT_MINE        setmine_
#    define FORT_FASTMARCH2  fastmarch2_
#if (BL_SPACEDIM == 3)
#    define FORT_SUMTURB     sumturb_
#ifdef SUMJET
#    define FORT_SUMJET      sumjet_
#endif
#endif
#  endif

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" {
#endif

    int FORT_PHIUPD(const amrex::Real * phi,ARLIM_P(phi_lo),ARLIM_P(phi_hi),
                     amrex::Real * phin,ARLIM_P(phin_lo),ARLIM_P(phin_hi),
                     const amrex::Real * uadv,ARLIM_P(u_lo),ARLIM_P(u_hi),
                     const amrex::Real * vadv,ARLIM_P(v_lo),ARLIM_P(v_hi),
#if (BL_SPACEDIM == 3)
                     const amrex::Real * wadv,ARLIM_P(w_lo),ARLIM_P(w_hi),
#endif
                     int * nband, int * nbandsize,
                     int * mine, int * minesize,
                     const int* lo, const int* hi, amrex::Real* dt, const amrex::Real* dx,
                     int * type, ARLIM_P(type_lo), ARLIM_P(type_hi));

    amrex::Real FORT_LSCFL(const amrex::Real * phi, ARLIM_P(phi_lo),ARLIM_P(phi_hi),
                    const amrex::Real * uadv,ARLIM_P(u_lo),ARLIM_P(u_hi),
                    const amrex::Real * vadv,ARLIM_P(v_lo),ARLIM_P(v_hi),
#if (BL_SPACEDIM == 3)
                    const amrex::Real * wadv,ARLIM_P(w_lo),ARLIM_P(w_hi),
#endif
                    int * nband, int * nbandsize,
                    int * mine, int * minesize,
                    const int* lo, const int* hi,amrex::Real* phit,const amrex::Real* dx,
                    int * type, ARLIM_P(type_lo), ARLIM_P(type_hi));

    void FORT_NARROWBAND(int * type, ARLIM_P(type_lo), ARLIM_P(type_hi),
     		             int * nband, int * nbandsize,
     		             int * mine, int * minesize,
      		             const int * lo, const int * hi);

    void FORT_FINDINTRFCE(const amrex::Real * phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
    					  amrex::Real * phin,ARLIM_P(phin_lo), ARLIM_P(phin_hi),
	                      int * type, ARLIM_P(type_lo), ARLIM_P(type_hi),
	                      const int* lo, const int* hi,const amrex::Real *  dx,
	                      int * intfacenump, int * intfacenumn, int * intfacep, int * intfacen,
                          int * nband, int * nbandsize, int * intfacesize);

    void FORT_RETYPIFY(int * type, ARLIM_P(type_lo), ARLIM_P(type_hi),
                       int * nband, int * nbandsize);

    void  FORT_FASTMARCH(amrex::Real * phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
                         int * type, ARLIM_P(type_lo), ARLIM_P(type_hi),
                         const int * lo, const int * hi, const amrex::Real * dx,
                         int * intfacenum, int* intface,
                         int * nband, int* nbandsize, int* nbandnum,
                         int* mine, int* sgn, int* intfacesize, int* heap,
                         int * heaploc);


    void FORT_NBANDNUMIFY(int * nband,int* nbandsize,int* nbandnum );


    void FORT_MINE(int* type, ARLIM_P(type_lo), ARLIM_P(type_hi),
	       		   int* nband,int* nbandsize,
	               int* mine,int* minesize,
     		       const int *lo, const int * hi);

    int FORT_FASTMARCH2(amrex::Real * phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
                        int * type, ARLIM_P(type_lo), ARLIM_P(type_hi),
                        const int* lo, const int* hi,const amrex::Real * dx,
                        int *    nband, int *nbandsize, int *nbandnum,
                        int * sgn, int* heaploc);
#ifdef __cplusplus
}
#endif

#endif
#endif 


