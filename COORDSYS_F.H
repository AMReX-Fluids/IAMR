#ifndef _COORDSYS_F_H_
#define _COORDSYS_F_H_  

#ifdef BL_LANG_FORT
#    define FORT_SETVOL      setvol
#    define FORT_SETAREA     setarea
#    define FORT_SETDLOGA    setdloga
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_SETVOL      SETVOL
#    define FORT_SETAREA     SETAREA
#    define FORT_SETDLOGA    SETDLOGA
#  else
#    define FORT_SETVOL      setvol_
#    define FORT_SETAREA     setarea_
#    define FORT_SETDLOGA    setdloga_
#  endif

#include "ArrayLim.H"

  extern "C" {
     void FORT_SETVOL (REAL* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       const REAL* offset,
                       const REAL* dx, const int* coord);

#if (BL_SPACEDIM == 2)
     void FORT_SETDLOGA (REAL* vol, ARLIM_P(vlo), ARLIM_P(vhi),
			 const REAL* offset,
			 const REAL* dx, const int* dir,
			 const int* coord);
#endif

     void FORT_SETAREA (REAL* area,  ARLIM_P(alo), ARLIM_P(ahi),
                        const REAL* offset,
                        const REAL* dx, const int* dir,
                        const int* coord);
  }

#endif

#endif
