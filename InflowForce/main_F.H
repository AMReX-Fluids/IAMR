#ifndef _main_F_H_
#define _main_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_READ_SWIRL_HEADER  read_swirl_header
#    define FORT_FILL fill
#    define FORT_WRITE writef
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_READ_SWIRL_HEADER  READ_SWIRL_HEADER
#    define FORT_FILL FILL
#    define FORT_WRITE WRITEF
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_READ_SWIRL_HEADER  read_swirl_header
#    define FORT_FILL fill
#    define FORT_WRITE writef
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_READ_SWIRL_HEADER  read_swirl_header_
#    define FORT_FILL fill_
#    define FORT_WRITE writef_
#  endif

  extern "C"
  {
      void FORT_READ_SWIRL_HEADER(int*       fname,
                                  int*       flen,
                                  int*       imax,
                                  int*       jmax,
                                  int*       kmax,
                                  int*       ncomp,
                                  Real*      scalex,
                                  Real*      scaley,
                                  Real*      scalez,
                                  Real*      times,
                                  const int* NN);

      void FORT_FILL(Real*       dat,
                     const int*  imax,
                     const int*  jmax,
                     const int*  kmax,
                     const int*  ncomp,
                     const Real* time,
                     const int*  putidx,
                     const int*  getidx,
                     int*        fname,
                     int*        flen,
                     int*        count);

      void FORT_WRITE(Real*       dat,
                      const int*  imax,
                      const int*  jmax,
                      const int*  kmax,
                      const int*  ncomp,
                      const Real* times,
                      const Real* scalex,
                      const Real* scaley,
                      const Real* scalez,
                      int*        fname,
                      int*        flen);
  }
#endif

#endif
