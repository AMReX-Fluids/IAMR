
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

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" {
#endif
#if (BL_SPACEDIM == 2)
    void FORT_HGPHIBC(const amrex::Real* hx,   const amrex::Real* sigE,
		      const amrex::Real* s, const amrex::Real* x, 
                      const int* length, const int* isPeriodicInX);

    void FORT_OLDHGPHIBC(const amrex::Real* hx,   const amrex::Real* rcen,
                         const amrex::Real* uE,   const amrex::Real* divE,
		         const amrex::Real* rhoE, const amrex::Real* x, 
		         const amrex::Real* flag,
                         const int* length, const int* isPeriodicInX);

    void FORT_FILL_ONED(int* lenx, int* leny, int* length, 
                        int * faces, int * numOutFlowFaces, 
                        const amrex::Real* cc0, const amrex::Real* cc1, 
                        const amrex::Real* cc2, const amrex::Real* cc3,
                        const amrex::Real*  r0, const amrex::Real*  r1, 
                        const amrex::Real*  r2, const amrex::Real*  r3,
                        amrex::Real * cc_conn, amrex::Real* s, const int* isPeriodic,
                        const amrex::Real* hx, const amrex::Real* hy);

    void FORT_ALLPHI_FROM_X(int* lenx, int* leny, int* length, 
                            int * faces, int * numOutFlowFaces,
                            const amrex::Real* phi0, const amrex::Real* phi1, 
                            const amrex::Real* phi2, const amrex::Real* phi3,
                            amrex::Real * x);
#elif (BL_SPACEDIM == 3)

    void FORT_FILL_TWOD(int* lenx, int* leny, int* lenz, int* length, int* width,
                        int * faces, int * numOutFlowFaces, 
                        const amrex::Real* cc0, const amrex::Real* cc1, const amrex::Real* cc2,
                        const amrex::Real* cc3, const amrex::Real* cc4, const amrex::Real* cc5,
                        amrex::Real * cc_conn);

    void FORT_ALLPHI_FROM_X(int* lenx, int* leny, int* lenz, int* length, int* width,
                            int * faces, int * numOutFlowFaces,
                            const amrex::Real* phi0, const amrex::Real* phi1, const amrex::Real* phi2,
                            const amrex::Real* phi3, const amrex::Real* phi4, const amrex::Real* phi5,
                            amrex::Real * x, ARLIM_P(x_lo), ARLIM_P(x_hi));

#endif

    void FORT_PHI_FROM_X(ARLIM_P(phi_lo), ARLIM_P(phi_hi), const amrex::Real* phi,
			 int * length, amrex::Real * x, int * isPeriodicInX);

    void FORT_RHOGBC(const amrex::Real* rho, ARLIM_P(rho_lo), ARLIM_P(rho_hi),
                     const amrex::Real* phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi), 
		     int * face, const amrex::Real* gravity, const amrex::Real* dx,
                     const int* domlo, const int* domhi,
                     const int* lo_bc, const int* hi_bc);

    void  FORT_COMPUTE_COEFF(ARLIM_P(rhslo),ARLIM_P(rhshi),const amrex::Real* rhsPtr,
			ARLIM_P(betalo),ARLIM_P(betahi),const amrex::Real* betaPtr,
			ARLIM_P(uElo), ARLIM_P(uEhi), const amrex::Real* uExtPtr,
			ARLIM_P(divuElo),ARLIM_P(divuEhi), const amrex::Real* divuEPtr,
			ARLIM_P(rhoElo),ARLIM_P(rhoEhi),const amrex::Real* rhoEPtr,
#if (BL_SPACEDIM == 2)
			int* r_lo,int* r_hi,amrex::Real* rcen,
#endif
			  const int* faceLo, const int* faceHi,
			const amrex::Real* dx,const int* isPeriodicPerp);

    void FORT_EXTRAP_PROJ(ARLIM_P(u_lo), ARLIM_P(u_hi), const amrex::Real* u,
			  ARLIM_P(d_lo), ARLIM_P(d_hi), const amrex::Real* div,
			  ARLIM_P(rho_lo), ARLIM_P(rho_hi), const amrex::Real* rho,
#if (BL_SPACEDIM == 2)
			  int* r_len, const amrex::Real* redge,
#endif
			  ARLIM_P(uE_lo), ARLIM_P(uE_hi), const amrex::Real* uE,
			  ARLIM_P(dE_lo), ARLIM_P(dE_hi), const amrex::Real* divE,
			  ARLIM_P(rhoE_lo), ARLIM_P(rhoE_hi), const amrex::Real* rhoE,
			  const int* faceLo, const int* faceHi,
			  int* face, int* zeroIt, const amrex::Real* hx);

    void FORT_HGRELAX (ARLIM_P(rhslo), ARLIM_P(rhshi), const amrex::Real* rhs,
                       ARLIM_P(betalo), ARLIM_P(betahi), const amrex::Real* beta,
                       ARLIM_P(philo), ARLIM_P(phihi), const amrex::Real* phi,
                       ARLIM_P(dgphilo), ARLIM_P(dgphihi), const amrex::Real* dgphi,
		       const int* faceLo, const int* faceHi,
		       const amrex::Real* h, const int* isPeriodic,
		       const int* niter);

    void FORT_HG_RESHIFT_PHI(ARLIM_P(outlo),ARLIM_P(outhi),const amrex::Real* outPtr,
			   ARLIM_P(inlo),ARLIM_P(inhi),const amrex::Real* inPtr,
			   const int* face);

    void FORT_HG_SHIFT_PHI(ARLIM_P(outlo),ARLIM_P(outhi),const amrex::Real* outPtr,
                           ARLIM_P(inlo),ARLIM_P(inhi),const amrex::Real* inPtr,
                           const int* face);


    void FORT_HGRESID (ARLIM_P(rhslo), ARLIM_P(rhshi), const amrex::Real* rhs,
                       ARLIM_P(betalo), ARLIM_P(betahi), const amrex::Real* beta,
                       ARLIM_P(philo), ARLIM_P(phihi), const amrex::Real* phi,
                       ARLIM_P(residlo), ARLIM_P(residhi), const amrex::Real* resid,
                       ARLIM_P(dgphilo), ARLIM_P(dgphihi), const amrex::Real* dgphi,
		       const int* faceLo, const int* faceHi,
		       const amrex::Real* h, const int* isPeriodic, 
		       const amrex::Real* maxnorm);
    
    void FORT_HGSUBTRACTAVGPHI(ARLIM_P(philo),ARLIM_P(phihi),const amrex::Real* phi,
#if (BL_SPACEDIM == 2)
			       int* r_lo, int* r_hi, const amrex::Real* r,
#endif
			       const int* faceLo, const int* faceHi,
			       const int* isPeriodic);

  void FORT_SOLVEHG(amrex::Real * dest, ARLIM_P(dest_lo),ARLIM_P(dest_hi), 
                    amrex::Real * dest0, ARLIM_P(dest0_lo),ARLIM_P(dest0_hi), 
                    amrex::Real * source, ARLIM_P(source_lo),ARLIM_P(source_hi), 
                    amrex::Real * sigma,  ARLIM_P(sigma_lo),ARLIM_P(sigma_hi), 
                    amrex::Real * cen, ARLIM_P(cen_lo),ARLIM_P(cen_hi), 
                    amrex::Real * r, ARLIM_P(r_lo),ARLIM_P(r_hi), 
                    amrex::Real * w, ARLIM_P(w_lo),ARLIM_P(w_hi), 
                    amrex::Real * z, ARLIM_P(z_lo),ARLIM_P(z_hi), 
                    amrex::Real * x, ARLIM_P(x_lo),ARLIM_P(x_hi), 
                    const int* lo, const int* hi, amrex::Real * dxFiltered,
                    int * isPeriodicFiltered, int * maxiter,
                    amrex::Real* tol, amrex::Real* abs_tol, amrex::Real* max_jump,amrex::Real * norm);

  void FORT_RESTRICT(amrex::Real * resid,ARLIM_P(resid_lo),ARLIM_P(resid_hi), 
                     amrex::Real * resid_crse, ARLIM_P(resc_lo),ARLIM_P(resc_hi),
                     const int* lo, const int* hi,
                     const int* loc, const int* hic,
                     int * isPeriodic);

  void FORT_COARSIG(amrex::Real * sigma, ARLIM_P(si_lo),ARLIM_P(si_hi),
                    amrex::Real * sigmac,ARLIM_P(ns_lo),ARLIM_P(ns_hi),
                    const int* lo, const int* hi,
                    const int* loc, const int* hic,
                    int * isPeriodic);

  void FORT_INTERP(amrex::Real* phi, ARLIM_P(phi_lo),ARLIM_P(phi_hi),
		   amrex::Real * temp, ARLIM_P(temp_lo),ARLIM_P(temp_hi), 
                   amrex::Real * phi_crse, ARLIM_P(phic_lo),ARLIM_P(phic_hi), 
                   amrex::Real * sig, ARLIM_P(sig_lo),ARLIM_P(sig_hi), 
                   const int* lo, const int* hi,
                   const int* loc, const int* hic,
		   const int* isPeriodic);


#ifdef __cplusplus
}
#endif
#endif
#endif
