#ifndef _mainIBD_F_H_
#define _mainIBD_F_H_

#include <ArrayLim.H>

#if	defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 2)
#define FORT_FILLCRSE  fillcrse2d
#define FORT_FILLFINE  fillfine2d
#endif

#if (BL_SPACEDIM == 3)
#define FORT_FILLCRSE  fillcrse3d
#define FORT_FILLFINE  fillfine3d
#endif

#endif


#if	defined(BL_LANG_CC)

#if (BL_SPACEDIM == 2)
#ifdef BL_FORT_USE_UPPERCASE
#define FORT_FILLCRSE  FILLCRSE2D
#define FORT_FILLFINE  FILLFINE2D
#endif
#ifdef BL_FORT_USE_UNDERSCORE
#define FORT_FILLCRSE  fillcrse2d_
#define FORT_FILLFINE  fillfine2d_
#endif
#endif

#if (BL_SPACEDIM == 3)
#ifdef BL_FORT_USE_UPPERCASE
#define FORT_FILLCRSE  FILLCRSE3D
#define FORT_FILLFINE  FILLFINE3D
#endif
#ifdef BL_FORT_USE_UNDERSCORE
#define FORT_FILLCRSE  fillcrse3d_
#define FORT_FILLFINE  fillfine3d_
#endif

#endif

extern "C"
{
    void FORT_FILLCRSE (
	REAL* crse, ARLIM_P(crse_lo), ARLIM_P(crse_hi),
	const REAL* h, const int* nc
	);
    void FORT_FILLFINE (
	REAL* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
	const REAL* h, const int* nc
	);
};

#endif
#endif
