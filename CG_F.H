#ifndef _CG_F_H_
#define _CG_F_H_ "$Id: CG_F.H,v 1.1 1997-07-08 23:08:03 vince Exp $"

#include <REAL.H>

#if	defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 2)
#define FORT_CGUPDATE  cgupdate2dgen
#define FORT_CGADVCP   cgadvcp2dgen
#define FORT_CGXDOTY   cgxdoty2dgen
#endif

#if (BL_SPACEDIM == 3)
#define FORT_CGUPDATE  cgupdate3dgen
#define FORT_CGADVCP   cgadvcp3dgen
#define FORT_CGXDOTY   cgxdoty3dgen
#endif

#endif

#if	defined(BL_LANG_CC)

#ifdef BL_FORT_USE_UPPERCASE

#if (BL_SPACEDIM == 2)
#define FORT_CGUPDATE  CGUPDATE2DGEN
#define FORT_CGADVCP   CGADVCP2DGEN
#define FORT_CGXDOTY   CGXDOTY2DGEN
#endif

#if (BL_SPACEDIM == 3)
#define FORT_CGUPDATE  CGUPDATE3DGEN
#define FORT_CGADVCP   CGADVCP3DGEN
#define FORT_CGXDOTY   CGXDOTY3DGEN
#endif

#endif
#ifdef BL_FORT_USE_UNDERSCORE

#if (BL_SPACEDIM == 2)
#define FORT_CGUPDATE  cgupdate2dgen_
#define FORT_CGADVCP   cgadvcp2dgen_
#define FORT_CGXDOTY   cgxdoty2dgen_
#endif

#if (BL_SPACEDIM == 3)
#define FORT_CGUPDATE  cgupdate3dgen_
#define FORT_CGADVCP   cgadvcp3dgen_
#define FORT_CGXDOTY   cgxdoty3dgen_
#endif

#endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_CGUPDATE (
	REAL* phi,  ARLIM_P(phi_lo),  ARLIM_P(phi_hi),
	REAL* tmpr, ARLIM_P(tmpr_lo), ARLIM_P(tmpr_hi),
	const REAL *h,
	const REAL* tmpw, ARLIM_P(tmpw_lo), ARLIM_P(tmpw_hi),
	const REAL* tmpp, ARLIM_P(tmpp_lo), ARLIM_P(tmpp_hi),
	const int* lo, const int* hi,
	const int* nc
	);

    void FORT_CGADVCP (
	REAL* tmpp,       ARLIM_P(tmpp_lo), ARLIM_P(tmpp_hi),
	const REAL* tmpr, ARLIM_P(tmpr_lo), ARLIM_P(tmpr_hi),
	const REAL *b,
	const int* lo, const int* hi,
	const int* nc
	);

    void FORT_CGXDOTY (
	REAL *pw,
	const REAL* phi,  ARLIM_P(phi_lo),  ARLIM_P(phi_hi),
	const REAL* tmpw, ARLIM_P(tmpw_lo), ARLIM_P(tmpw_hi),
	const int* lo, const int* hi,
	const int* nc
	);

}
#endif

#endif
