
#ifndef _PROJOUTFLOWBC_F_H_
#define _PROJOUTFLOWBC_F_H_ 

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" {
#endif
#if (BL_SPACEDIM == 2)
    void hgphibc(const amrex::Real* hx,   const amrex::Real* sigE,
		      const amrex::Real* s, const amrex::Real* x, 
                      const int* length, const int* isPeriodicInX);

    void oldhgphibc(const amrex::Real* hx,   const amrex::Real* rcen,
                         const amrex::Real* uE,   const amrex::Real* divE,
		         const amrex::Real* rhoE, const amrex::Real* x, 
		         const amrex::Real* flag,
                         const int* length, const int* isPeriodicInX);

    void fill_oned(int* lenx, int* leny, int* length, 
                        int * faces, int * numOutFlowFaces, 
                        const amrex::Real* cc0, const amrex::Real* cc1, 
                        const amrex::Real* cc2, const amrex::Real* cc3,
                        const amrex::Real*  r0, const amrex::Real*  r1, 
                        const amrex::Real*  r2, const amrex::Real*  r3,
                        amrex::Real * cc_conn, amrex::Real* s, const int* isPeriodic,
                        const amrex::Real* hx, const amrex::Real* hy);

    void allphi_from_x(int* lenx, int* leny, int* length, 
                            int * faces, int * numOutFlowFaces,
                            const amrex::Real* phi0, const amrex::Real* phi1, 
                            const amrex::Real* phi2, const amrex::Real* phi3,
                            amrex::Real * x);
#elif (BL_SPACEDIM == 3)

    void fill_twod(int* lenx, int* leny, int* lenz, int* length, int* width,
                        int * faces, int * numOutFlowFaces, 
                        const amrex::Real* cc0, const amrex::Real* cc1, const amrex::Real* cc2,
                        const amrex::Real* cc3, const amrex::Real* cc4, const amrex::Real* cc5,
                        amrex::Real * cc_conn);

    void allphi_from_x(int* lenx, int* leny, int* lenz, int* length, int* width,
                            int * faces, int * numOutFlowFaces,
                            const amrex::Real* phi0, const amrex::Real* phi1, const amrex::Real* phi2,
                            const amrex::Real* phi3, const amrex::Real* phi4, const amrex::Real* phi5,
                            amrex::Real * x, ARLIM_P(x_lo), ARLIM_P(x_hi));

#endif

    void phi_from_x(ARLIM_P(phi_lo), ARLIM_P(phi_hi), const amrex::Real* phi,
			 int * length, amrex::Real * x, int * isPeriodicInX);

    void rhogbc(const amrex::Real* rho, ARLIM_P(rho_lo), ARLIM_P(rho_hi),
                     const amrex::Real* phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi), 
		     int * face, const amrex::Real* gravity, const amrex::Real* dx,
                     const int* domlo, const int* domhi,
                     const int* lo_bc, const int* hi_bc);

    void  compute_coeff(ARLIM_P(rhslo),ARLIM_P(rhshi),const amrex::Real* rhsPtr,
			ARLIM_P(betalo),ARLIM_P(betahi),const amrex::Real* betaPtr,
			ARLIM_P(uElo), ARLIM_P(uEhi), const amrex::Real* uExtPtr,
			ARLIM_P(divuElo),ARLIM_P(divuEhi), const amrex::Real* divuEPtr,
			ARLIM_P(rhoElo),ARLIM_P(rhoEhi),const amrex::Real* rhoEPtr,
#if (BL_SPACEDIM == 2)
			int* r_lo,int* r_hi,amrex::Real* rcen,
#endif
			  const int* faceLo, const int* faceHi,
			const amrex::Real* dx,const int* isPeriodicPerp);

    void extrap_proj(ARLIM_P(u_lo), ARLIM_P(u_hi), const amrex::Real* u,
			  ARLIM_P(d_lo), ARLIM_P(d_hi), const amrex::Real* div,
			  ARLIM_P(rho_lo), ARLIM_P(rho_hi), const amrex::Real* rho,
#if (BL_SPACEDIM == 2)
			  int* r_len, const amrex::Real* redge,
#endif
			  ARLIM_P(uE_lo), ARLIM_P(uE_hi), const amrex::Real* uE,
			  ARLIM_P(dE_lo), ARLIM_P(dE_hi), const amrex::Real* divE,
			  ARLIM_P(rhoE_lo), ARLIM_P(rhoE_hi), const amrex::Real* rhoE,
			  const int* faceLo, const int* faceHi,
			  int* face, int* zeroIt);

    void hgrelax (ARLIM_P(rhslo), ARLIM_P(rhshi), const amrex::Real* rhs,
                       ARLIM_P(betalo), ARLIM_P(betahi), const amrex::Real* beta,
                       ARLIM_P(philo), ARLIM_P(phihi), const amrex::Real* phi,
                       ARLIM_P(dgphilo), ARLIM_P(dgphihi), const amrex::Real* dgphi,
		       const int* faceLo, const int* faceHi,
		       const amrex::Real* h, const int* isPeriodic,
		       const int* niter);

    void hg_reshift_phi(ARLIM_P(outlo),ARLIM_P(outhi),const amrex::Real* outPtr,
			   ARLIM_P(inlo),ARLIM_P(inhi),const amrex::Real* inPtr,
			   const int* face);

    void hg_shift_phi(ARLIM_P(outlo),ARLIM_P(outhi),const amrex::Real* outPtr,
                           ARLIM_P(inlo),ARLIM_P(inhi),const amrex::Real* inPtr,
                           const int* face);


    void hgresid (ARLIM_P(rhslo), ARLIM_P(rhshi), const amrex::Real* rhs,
                       ARLIM_P(betalo), ARLIM_P(betahi), const amrex::Real* beta,
                       ARLIM_P(philo), ARLIM_P(phihi), const amrex::Real* phi,
                       ARLIM_P(residlo), ARLIM_P(residhi), const amrex::Real* resid,
                       ARLIM_P(dgphilo), ARLIM_P(dgphihi), const amrex::Real* dgphi,
		       const int* faceLo, const int* faceHi,
		       const amrex::Real* h, const int* isPeriodic, 
		       const amrex::Real* maxnorm);
    
    void hgsubtractavgphi(ARLIM_P(philo),ARLIM_P(phihi),const amrex::Real* phi,
#if (BL_SPACEDIM == 2)
			       int* r_lo, int* r_hi, const amrex::Real* r,
#endif
			       const int* faceLo, const int* faceHi,
			       const int* isPeriodic);

  void solvehg(amrex::Real * dest, ARLIM_P(dest_lo),ARLIM_P(dest_hi), 
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

  void fort_restrict(amrex::Real * resid,ARLIM_P(resid_lo),ARLIM_P(resid_hi), 
                     amrex::Real * resid_crse, ARLIM_P(resc_lo),ARLIM_P(resc_hi),
                     const int* lo, const int* hi,
                     const int* loc, const int* hic,
                     int * isPeriodic);

  void coarsig(amrex::Real * sigma, ARLIM_P(si_lo),ARLIM_P(si_hi),
                    amrex::Real * sigmac,ARLIM_P(ns_lo),ARLIM_P(ns_hi),
                    const int* lo, const int* hi,
                    const int* loc, const int* hic,
                    int * isPeriodic);

  void interp(amrex::Real* phi, ARLIM_P(phi_lo),ARLIM_P(phi_hi),
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
