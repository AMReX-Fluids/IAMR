#ifndef _PROJOUTFLOWBC_F_H_
#define _PROJOUTFLOWBC_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_OLDHGPHIBC    oldhgphibc
#    define FORT_HGPHIBC       hgphibc
#    define FORT_PHI_FROM_X    phifromx
#    define FORT_ALLPHI_FROM_X allphifromx
#    define FORT_FILL_ONED     filloned
#    define FORT_FILL_TWOD     filltwod
#    define FORT_RHOGBC        rhogbc
#    define FORT_EXTRAP_PROJ   extrapproj
#    define FORT_COMPUTE_COEFF computecoeff
#    define FORT_HGRELAX       hgrelax
#    define FORT_HG_RESHIFT_PHI   hgreshiftphi
#    define FORT_HG_SHIFT_PHI   hgshiftphi
#    define FORT_HGRESID       hgresid
#    define FORT_HGSUBTRACTAVGPHI hgsubtractavgphi
#    define FORT_SOLVEHG        solvehg
#  define FORT_COARSIG     coarsigb
#  define FORT_RESTRICT    coarsresb
#  define FORT_INTERP      interphg
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_HGPHIBC       HGPHIBC
#    define FORT_OLDHGPHIBC    OLDHGPHIBC
#    define FORT_RHOGBC        RHOGBC
#    define FORT_PHI_FROM_X    PHIFROMX
#    define FORT_ALLPHI_FROM_X ALLPHIFROMX
#    define FORT_FILL_ONED     FILLONED
#    define FORT_FILL_TWOD     FILLTWOD
#    define FORT_EXTRAP_PROJ   EXTRAPPROJ
#    define FORT_COMPUTE_COEFF COMPUTECOEFF
#    define FORT_HGRELAX       HGRELAX
#    define FORT_HG_SHIFT_PHI   HGSHIFTPHI
#    define FORT_HG_RESHIFT_PHI   HGRESHIFTPHI
#    define FORT_HGRESID       HGRESID
#    define FORT_HGSUBTRACTAVGPHI HGSUBTRACTAVGPHI
#    define FORT_SOLVEHG        SOLVEHG
#  define FORT_COARSIG     COARSIGB
#  define FORT_RESTRICT    COARSRESB
#  define FORT_INTERP      INTERPHG
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_HGPHIBC       hgphibc
#    define FORT_OLDHGPHIBC    oldhgphibc
#    define FORT_PHI_FROM_X    phifromx
#    define FORT_ALLPHI_FROM_X allphifromx
#    define FORT_FILL_ONED     filloned
#    define FORT_FILL_TWOD     filltwod
#    define FORT_RHOGBC        rhogbc
#    define FORT_EXTRAP_PROJ   extrapproj
#    define FORT_COMPUTE_COEFF computecoeff
#    define FORT_HGRELAX       hgrelax
#    define FORT_HG_RESHIFT_PHI  hgreshiftphi
#    define FORT_HG_SHIFT_PHI    hgshiftphi
#    define FORT_HGRESID       hgresid
#    define FORT_HGSUBTRACTAVGPHI hgsubtractavgphi
#    define FORT_SOLVEHG        solvehg
#  define FORT_COARSIG     coarsigb
#  define FORT_RESTRICT    coarsresb
#  define FORT_INTERP      interphg
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_HGPHIBC       hgphibc_
#    define FORT_OLDHGPHIBC    oldhgphibc_
#    define FORT_PHI_FROM_X    phifromx_
#    define FORT_ALLPHI_FROM_X allphifromx_
#    define FORT_FILL_ONED     filloned_
#    define FORT_FILL_TWOD     filltwod_
#    define FORT_RHOGBC        rhogbc_
#    define FORT_EXTRAP_PROJ   extrapproj_
#    define FORT_COMPUTE_COEFF computecoeff_
#    define FORT_HGRELAX       hgrelax_
#    define FORT_HG_SHIFT_PHI   hgshiftphi_
#    define FORT_HG_RESHIFT_PHI hgreshiftphi_
#    define FORT_HGRESID       hgresid_
#    define FORT_HGSUBTRACTAVGPHI hgsubtractavgphi_
#    define FORT_SOLVEHG        solvehg_
#  define FORT_COARSIG     coarsigb_
#  define FORT_RESTRICT    coarsresb_
#  define FORT_INTERP      interphg_
#  endif

#include <ArrayLim.H>

extern "C" {
#if (BL_SPACEDIM == 2)
    void FORT_HGPHIBC(const Real* hx,   const Real* sigE,
		      const Real* s, const Real* x, 
                      const int* length, const int* isPeriodicInX);

    void FORT_OLDHGPHIBC(const Real* hx,   const Real* rcen,
                         const Real* uE,   const Real* divE,
		         const Real* rhoE, const Real* x, 
		         const Real* flag,
                         const int* length, const int* isPeriodicInX);

    void FORT_FILL_ONED(int* lenx, int* leny, int* length, 
                        int * faces, int * numOutFlowFaces, 
                        const Real* cc0, const Real* cc1, 
                        const Real* cc2, const Real* cc3,
                        const Real*  r0, const Real*  r1, 
                        const Real*  r2, const Real*  r3,
                        Real * cc_conn, Real* s, const int* isPeriodic,
                        const Real* hx, const Real* hy);

    void FORT_ALLPHI_FROM_X(int* lenx, int* leny, int* length, 
                            int * faces, int * numOutFlowFaces,
                            const Real* phi0, const Real* phi1, 
                            const Real* phi2, const Real* phi3,
                            Real * x);
#elif (BL_SPACEDIM == 3)

    void FORT_FILL_TWOD(int* lenx, int* leny, int* lenz, int* length, int* width,
                        int * faces, int * numOutFlowFaces, 
                        const Real* cc0, const Real* cc1, const Real* cc2,
                        const Real* cc3, const Real* cc4, const Real* cc5,
                        Real * cc_conn);

    void FORT_ALLPHI_FROM_X(int* lenx, int* leny, int* lenz, int* length, int* width,
                            int * faces, int * numOutFlowFaces,
                            const Real* phi0, const Real* phi1, const Real* phi2,
                            const Real* phi3, const Real* phi4, const Real* phi5,
                            Real * x, ARLIM_P(x_lo), ARLIM_P(x_hi));

#endif

    void FORT_PHI_FROM_X(ARLIM_P(phi_lo), ARLIM_P(phi_hi), const Real* phi,
			 int * length, Real * x, int * isPeriodicInX);

    void FORT_RHOGBC(const Real* rho, ARLIM_P(rho_lo), ARLIM_P(rho_hi),
                     const Real* phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi), 
		     int * face, const Real* gravity, const Real* dx,
                     const int* domlo, const int* domhi,
                     const int* lo_bc, const int* hi_bc);

    void  FORT_COMPUTE_COEFF(ARLIM_P(rhslo),ARLIM_P(rhshi),const Real* rhsPtr,
			ARLIM_P(betalo),ARLIM_P(betahi),const Real* betaPtr,
			ARLIM_P(uElo), ARLIM_P(uEhi), const Real* uExtPtr,
			ARLIM_P(divuElo),ARLIM_P(divuEhi), const Real* divuEPtr,
			ARLIM_P(rhoElo),ARLIM_P(rhoEhi),const Real* rhoEPtr,
#if (BL_SPACEDIM == 2)
			int* r_lo,int* r_hi,Real* rcen,
#endif
			  const int* faceLo, const int* faceHi,
			const Real* dx,const int* isPeriodicPerp);

    void FORT_EXTRAP_PROJ(ARLIM_P(u_lo), ARLIM_P(u_hi), const Real* u,
			  ARLIM_P(d_lo), ARLIM_P(d_hi), const Real* div,
			  ARLIM_P(rho_lo), ARLIM_P(rho_hi), const Real* rho,
#if (BL_SPACEDIM == 2)
			  int* r_len, const Real* redge,
#endif
			  ARLIM_P(uE_lo), ARLIM_P(uE_hi), const Real* uE,
			  ARLIM_P(dE_lo), ARLIM_P(dE_hi), const Real* divE,
			  ARLIM_P(rhoE_lo), ARLIM_P(rhoE_hi), const Real* rhoE,
			  const int* faceLo, const int* faceHi,
			  int* face, int* zeroIt, const Real* hx);

    void FORT_HGRELAX (ARLIM_P(rhslo), ARLIM_P(rhshi), const Real* rhs,
                       ARLIM_P(betalo), ARLIM_P(betahi), const Real* beta,
                       ARLIM_P(philo), ARLIM_P(phihi), const Real* phi,
                       ARLIM_P(dgphilo), ARLIM_P(dgphihi), const Real* dgphi,
		       const int* faceLo, const int* faceHi,
		       const Real* h, const int* isPeriodic,
		       const int* niter);

    void FORT_HG_RESHIFT_PHI(ARLIM_P(outlo),ARLIM_P(outhi),const Real* outPtr,
			   ARLIM_P(inlo),ARLIM_P(inhi),const Real* inPtr,
			   const int* face);

    void FORT_HG_SHIFT_PHI(ARLIM_P(outlo),ARLIM_P(outhi),const Real* outPtr,
                           ARLIM_P(inlo),ARLIM_P(inhi),const Real* inPtr,
                           const int* face);


    void FORT_HGRESID (ARLIM_P(rhslo), ARLIM_P(rhshi), const Real* rhs,
                       ARLIM_P(betalo), ARLIM_P(betahi), const Real* beta,
                       ARLIM_P(philo), ARLIM_P(phihi), const Real* phi,
                       ARLIM_P(residlo), ARLIM_P(residhi), const Real* resid,
                       ARLIM_P(dgphilo), ARLIM_P(dgphihi), const Real* dgphi,
		       const int* faceLo, const int* faceHi,
		       const Real* h, const int* isPeriodic, 
		       const Real* maxnorm);
    
    void FORT_HGSUBTRACTAVGPHI(ARLIM_P(philo),ARLIM_P(phihi),const Real* phi,
#if (BL_SPACEDIM == 2)
			       int* r_lo, int* r_hi, const Real* r,
#endif
			       const int* faceLo, const int* faceHi,
			       const int* isPeriodic);

  void FORT_SOLVEHG(Real * dest, ARLIM_P(dest_lo),ARLIM_P(dest_hi), 
                    Real * dest0, ARLIM_P(dest0_lo),ARLIM_P(dest0_hi), 
                    Real * source, ARLIM_P(source_lo),ARLIM_P(source_hi), 
                    Real * sigma,  ARLIM_P(sigma_lo),ARLIM_P(sigma_hi), 
                    Real * cen, ARLIM_P(cen_lo),ARLIM_P(cen_hi), 
                    Real * r, ARLIM_P(r_lo),ARLIM_P(r_hi), 
                    Real * w, ARLIM_P(w_lo),ARLIM_P(w_hi), 
                    Real * z, ARLIM_P(z_lo),ARLIM_P(z_hi), 
                    Real * x, ARLIM_P(x_lo),ARLIM_P(x_hi), 
                    const int* lo, const int* hi, Real * dxFiltered,
                    int * isPeriodicFiltered, int * maxiter,
                    Real* tol, Real* abs_tol, Real* max_jump,Real * norm);

  void FORT_RESTRICT(Real * resid,ARLIM_P(resid_lo),ARLIM_P(resid_hi), 
                     Real * resid_crse, ARLIM_P(resc_lo),ARLIM_P(resc_hi),
                     const int* lo, const int* hi,
                     const int* loc, const int* hic,
                     int * isPeriodic);

  void FORT_COARSIG(Real * sigma, ARLIM_P(si_lo),ARLIM_P(si_hi),
                    Real * sigmac,ARLIM_P(ns_lo),ARLIM_P(ns_hi),
                    const int* lo, const int* hi,
                    const int* loc, const int* hic,
                    int * isPeriodic);

  void FORT_INTERP(Real* phi, ARLIM_P(phi_lo),ARLIM_P(phi_hi),
		   Real * temp, ARLIM_P(temp_lo),ARLIM_P(temp_hi), 
                   Real * phi_crse, ARLIM_P(phic_lo),ARLIM_P(phic_hi), 
                   Real * sig, ARLIM_P(sig_lo),ARLIM_P(sig_hi), 
                   const int* lo, const int* hi,
                   const int* loc, const int* hic,
		   const int* isPeriodic);


}
#endif
#endif
