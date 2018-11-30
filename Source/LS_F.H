#ifndef _LS_F_H_
#define _LS_F_H_

#ifndef BL_LANG_FORT

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


