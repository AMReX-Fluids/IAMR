#ifndef _MACOUTFLOWBC_F_H
#define _MACOUTFLOWBC_F_H

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" {
#endif

#if (BL_SPACEDIM == 2)
    void macphibc (const amrex::Real* phi, int* length,
                        const amrex::Real* divuEptr, const amrex::Real* rhoEptr,
                        const amrex::Real* redge,
                        const amrex::Real* hx,
                        const int* isPeriodic);

    void macfill_oned(int* lenx, int* leny, int* length,
                           int * faces, int * numOutFlowFaces,
                           const amrex::Real* cc0, const amrex::Real* cc1,
                           const amrex::Real* cc2, const amrex::Real* cc3,
                           const amrex::Real*  r0, const amrex::Real*  r1,
                           const amrex::Real*  r2, const amrex::Real*  r3,
                           amrex::Real * cc_conn, amrex::Real* redge_conn);

    void macallphi_from_x(int* lenx, int* leny, int* length,
                            int * faces, int * numOutFlowFaces,
                            const amrex::Real* phi0, const amrex::Real* phi1,
                            const amrex::Real* phi2, const amrex::Real* phi3,
                            amrex::Real * x);

#else

    void  compute_maccoeff(ARLIM_P(rhslo),ARLIM_P(rhshi),const amrex::Real* rhsPtr,
				ARLIM_P(beta0lo),ARLIM_P(beta0hi),const amrex::Real* beta0Ptr,
				ARLIM_P(beta1lo),ARLIM_P(beta1hi),const amrex::Real* beta1Ptr,
				ARLIM_P(divuElo),ARLIM_P(divuEhi), 
				const amrex::Real* divuEPtr,
				ARLIM_P(rhoElo),ARLIM_P(rhoEhi),const amrex::Real* rhoEPtr,
				const int* faceLo, const int* faceHi,
				const amrex::Real* dx, const int* isPeriodicPerp);

    void macfill_twod(int* lenx, int* leny, int* length, int* width,
                           int * faces, int * numOutFlowFaces,
                           const amrex::Real* cc0, const amrex::Real* cc1, const amrex::Real* cc2,
                           const amrex::Real* cc3, const amrex::Real* cc4, const amrex::Real* cc5,
                           amrex::Real * cc_conn);


    void macallphi_from_x(int* lenx, int* leny, int* lenz, int* length, int* width,
                            int * faces, int * numOutFlowFaces,
                            const amrex::Real* phi0, const amrex::Real* phi1, const amrex::Real* phi2,
                            const amrex::Real* phi3, const amrex::Real* phi4, const amrex::Real* phi5,
                            amrex::Real * x, ARLIM_P(x_lo), ARLIM_P(x_hi));

#endif

    void macphi_from_x(ARLIM_P(phi_lo), ARLIM_P(phi_hi), const amrex::Real* phi,
                         int * length, amrex::Real * x);

    void extrap_mac (ARLIM_P(u0_lo), ARLIM_P(u0_hi), const amrex::Real* u0,
			  ARLIM_P(u1_lo), ARLIM_P(u1_hi), const amrex::Real* u1,
#if (BL_SPACEDIM == 3)
			  ARLIM_P(u2_lo), ARLIM_P(u2_hi), const amrex::Real* u2,
#endif
			  ARLIM_P(d_lo), ARLIM_P(d_hi), const amrex::Real* divcor,
			  ARLIM_P(rho_lo), ARLIM_P(rho_hi), const amrex::Real* rho,
#if (BL_SPACEDIM == 2)
			  int* r_len, const amrex::Real* redge,
#endif
			  ARLIM_P(dE_lo), ARLIM_P(dE_hi), const amrex::Real* dE,
			  ARLIM_P(rhoE_lo), ARLIM_P(rhoE_hi), const amrex::Real* rhoE,
			  const amrex::Real* dx,
			  const int* faceLo, const int* faceHi,
			  int* face, int* per, int* zeroIt, amrex::Real* small_udiff);

    void macrelax (ARLIM_P(rhslo), ARLIM_P(rhshi), const amrex::Real* rhs,
			ARLIM_P(beta0lo), ARLIM_P(beta0hi), const amrex::Real* beta0,
#if (BL_SPACEDIM == 3)
			ARLIM_P(beta1lo), ARLIM_P(beta1hi), const amrex::Real* beta1,
#endif
			ARLIM_P(philo), ARLIM_P(phihi), const amrex::Real* phi,
			const int* faceLo, const int* faceHi,
			const amrex::Real* h, const int* isPeriodic, 
			const int* niter);

    void mac_reshift_phi(
                           ARLIM_P(outlo),ARLIM_P(outhi),const amrex::Real* outPtr,
			   ARLIM_P(inlo),ARLIM_P(inhi),const amrex::Real* inPtr,
			   const int* face);

   void mac_shift_phi(
                           ARLIM_P(outlo),ARLIM_P(outhi),const amrex::Real* outPtr,
                           ARLIM_P(inlo),ARLIM_P(inhi),const amrex::Real* inPtr,
                           const int* face);


    void macresid(ARLIM_P(rhslo),ARLIM_P(rhshi),const amrex::Real* rhsPtr,
		       ARLIM_P(beta0),ARLIM_P(beta0),const amrex::Real* beta0,
#if (BL_SPACEDIM == 3)
		       ARLIM_P(beta1),ARLIM_P(beta1),const amrex::Real* beta1,
#endif
		       ARLIM_P(phi), ARLIM_P(phi),const amrex::Real* phi,
		       ARLIM_P(resid),ARLIM_P(resid), const amrex::Real* resid,
		       const int* lo, const int* hi, const amrex::Real* h,
		       const int* isPeriodic, const amrex::Real* maxnorm);


    void macsubtractavgphi(ARLIM_P(philo),ARLIM_P(phihi),const amrex::Real* phi,
#if (BL_SPACEDIM == 2)
		       int* r_lo, int* r_hi, const amrex::Real* r,
#endif
		       const int* faceLo, const int* faceHi,
				const int* isPeriodic);

  void solvemac(amrex::Real * dest, ARLIM_P(destlo),ARLIM_P(desthi),
		     amrex::Real * dest0, ARLIM_P(dest0lo),ARLIM_P(dest0hi),
		     amrex::Real *source, ARLIM_P(sourcelo),ARLIM_P(sourcehi),
                     amrex::Real * sigmax,  ARLIM_P(sigmaxlo),ARLIM_P(simgaxhi),
#if (BL_SPACEDIM == 3)
		     amrex::Real * sigmay, ARLIM_P(sigmaylo),ARLIM_P(sigmayhi),
#endif
		     amrex::Real * sum, ARLIM_P(sumlo),ARLIM_P(sumhi), 
                     amrex::Real * r, ARLIM_P(rlo),ARLIM_P(rhi), 
		     amrex::Real * w, ARLIM_P(wlo),ARLIM_P(whi), 
		     amrex::Real * z, ARLIM_P(zlo),ARLIM_P(zhi), 
		     amrex::Real * work,ARLIM_P(worklo),ARLIM_P(workhi), 
                     const int* lo, const int* hi,
                     amrex::Real* dxFiltered, int* isPeriodicFiltered,
                     int* maxIter, amrex::Real* tol, 
		     amrex::Real* cg_abs_tol, amrex::Real* max_jump,amrex::Real * norm);

  void coarsigma(amrex::Real * sigmax , ARLIM_P(sigmax_lo),ARLIM_P(sigmax_hi), 
#if (BL_SPACEDIM == 3)
		      amrex::Real * sigmay, ARLIM_P(sigmay_lo),ARLIM_P(sigmay_hi), 
#endif
                      amrex::Real * sigmaxc, ARLIM_P(sigmaxc_lo),ARLIM_P(sigmaxc_hi), 
#if (BL_SPACEDIM == 3)
		      amrex::Real * sigmayc,ARLIM_P(sigmayc_lo),ARLIM_P(sigmayc_hi), 
#endif
                      const int* lo, const int* hi,
		      const int* loc, const int* hic);

  void outflowbc_restrict(amrex::Real * res, ARLIM_P(res_lo),ARLIM_P(res_hi), 
		     amrex::Real * resc, ARLIM_P(resc_lo),ARLIM_P(resc_hi), 
                     const int* lo, const int* hi,
		     const int* loc, const int* hic);

  void fort_interpolate(amrex::Real * phi, ARLIM_P(phi_lo),ARLIM_P(phi_hi), 
			amrex::Real * deltac,ARLIM_P(deltac_lo),ARLIM_P(deltac_hi), 
			const int* lo, const int* hi,
			const int* loc, const int* hic);


#ifdef __cplusplus
}
#endif

#endif
#endif 
