#ifndef _AABBEC_F_H_
#define _AABBEC_F_H_ "$Id"

#include <REAL.H>

#if	defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 2)
#define FORT_GSRB	  gsrb2daabbec
#define FORT_CGPRECND	  cgprecnd2daabbec
#define FORT_ADOTX	  adotx2daabbec
#endif

#if (BL_SPACEDIM == 3)
#define FORT_GSRB	  gsrb3daabbec
#define FORT_CGPRECND	  cgprecnd3daabbec
#define FORT_ADOTX	  adotx3daabbec
#endif

#endif


#if	defined(BL_LANG_CC)

#if (BL_SPACEDIM == 2)

#ifdef BL_FORT_USE_UPPERCASE
#define FORT_GSRB	  GSRB2DAABBEC
#define FORT_CGPRECND  CGPRECND2DAABBEC
#define FORT_ADOTX  ADOTX2DAABBEC
#endif

#ifdef BL_FORT_USE_UNDERSCORE
#define FORT_GSRB	  gsrb2daabbec_
#define FORT_CGPRECND  cgprecnd2daabbec_
#define FORT_ADOTX  adotx2daabbec_
#endif

#endif


#if (BL_SPACEDIM == 3)

#ifdef BL_FORT_USE_UPPERCASE
#define FORT_GSRB	  GSRB3DAABBEC
#define FORT_CGPRECND  CGPRECND3DAABBEC
#define FORT_ADOTX  ADOTX3DAABBEC
#endif

#ifdef BL_FORT_USE_UNDERSCORE
#define FORT_GSRB	  gsrb3daabbec_
#define FORT_CGPRECND  cgprecnd3daabbec_
#define FORT_ADOTX  adotx3daabbec_
#endif

#endif

#include <ArrayLim.H>

extern "C"
{
#if (BL_SPACEDIM == 2)
    void FORT_GSRB (
	REAL* phi       , ARLIM_P(phi_lo), ARLIM_P(phi_hi),
	const REAL* rhs , ARLIM_P(rhs_lo), ARLIM_P(phi_hi),
	const REAL* alpha, const REAL* beta,
	const REAL* a   , ARLIM_P(a_lo),   ARLIM_P(a_hi),
	const REAL* bX  , ARLIM_P(bX_lo),  ARLIM_P(bX_hi),
	const REAL* bY  , ARLIM_P(bY_lo),  ARLIM_P(bY_hi),
	const REAL* den0, ARLIM_P(den0_lo),ARLIM_P(den0_hi),
	const int* m0   , ARLIM_P(m0_lo),  ARLIM_P(m0_hi),
	const REAL* den1, ARLIM_P(den1_lo),ARLIM_P(den1_hi),
	const int* m1   , ARLIM_P(m1_lo),  ARLIM_P(m1_hi),
	const REAL* den2, ARLIM_P(den2_lo),ARLIM_P(den2_hi),
	const int* m2   , ARLIM_P(m2_lo),  ARLIM_P(m2_hi),
	const REAL* den3, ARLIM_P(den3_lo),ARLIM_P(den3_hi),
	const int* m3   , ARLIM_P(m3_lo),  ARLIM_P(m3_hi),
	const int* lo, const int* hi, const int *nc,
	const REAL *h, const  int* redblack
	);

    void FORT_CGPRECND(
	REAL *tmpz      , ARLIM_P(tmpz_lo), ARLIM_P(tmpz_hi),
	REAL *rho,
	const REAL *tmpr, ARLIM_P(tmpr_lo), ARLIM_P(tmpr_hi),
	const REAL* alpha, const REAL* beta,
	const REAL* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
	const REAL* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
	const REAL* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
	const int *lo, const int *hi, const int *nc,
	const REAL *h
	);

    void FORT_ADOTX(
	REAL *y      , ARLIM_P(y_lo), ARLIM_P(y_hi),
	const REAL *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
	const REAL* alpha, const REAL* beta,
	const REAL* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
	const REAL* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
	const REAL* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
	const int *lo, const int *hi, const int *nc,
	const REAL *h
	);
#endif    

#if (BL_SPACEDIM == 3)
    void FORT_GSRB (
	REAL* phi,       ARLIM_P(phi_lo), ARLIM_P(phi_hi),
	const REAL* rhs, ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
	const REAL* alpha, const REAL* beta,
	const REAL* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
	const REAL* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
	const REAL* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
	const REAL* bZ, ARLIM_P(bZ_lo), ARLIM_P(bZ_hi),
	const REAL* den0, ARLIM_P(den0_lo), ARLIM_P(den0_hi),
	const int* m0   , ARLIM_P(m0_lo),   ARLIM_P(m0_hi),
	const REAL* den1, ARLIM_P(den1_lo), ARLIM_P(den1_hi),
	const int* m1   , ARLIM_P(m1_lo),   ARLIM_P(m1_hi),
	const REAL* den2, ARLIM_P(den2_lo), ARLIM_P(den2_hi),
	const int* m2   , ARLIM_P(m2_lo),   ARLIM_P(m2_hi),
	const REAL* den3, ARLIM_P(den3_lo), ARLIM_P(den3_hi),
	const int* m3   , ARLIM_P(m3_lo),   ARLIM_P(m3_hi),
	const REAL* den4, ARLIM_P(den4_lo), ARLIM_P(den4_hi),
	const int* m4   , ARLIM_P(m4_lo),   ARLIM_P(m4_hi),
	const REAL* den5, ARLIM_P(den5_lo), ARLIM_P(den5_hi),
	const int* m5   , ARLIM_P(m5_lo),   ARLIM_P(m5_hi),
	const int* lo, const int* hi, const int *nc,
	const REAL *h, const  int* redblack
	);

    void FORT_CGPRECND(
	REAL *tmpz      , ARLIM_P(tmpz_lo), ARLIM_P(tmpz_hi),
	REAL *rho,
	const REAL *tmpr, ARLIM_P(tmpr_lo), ARLIM_P(tmpr_hi),
	const REAL* alpha, const REAL* beta,
	const REAL* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
	const REAL* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
	const REAL* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
	const REAL* bZ, ARLIM_P(bZ_lo), ARLIM_P(bZ_hi),
	const int *lo, const int *hi, const int *nc,
	const REAL *h
	);

    void FORT_ADOTX(
	REAL *y      , ARLIM_P(y_lo), ARLIM_P(y_hi),
	const REAL *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
	const REAL* alpha, const REAL* beta,
	const REAL* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
	const REAL* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
	const REAL* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
	const REAL* bZ, ARLIM_P(bZ_lo), ARLIM_P(bZ_hi),
	const int *lo, const int *hi, const int *nc,
	const REAL *h
	);
#endif
}
#endif

#endif
