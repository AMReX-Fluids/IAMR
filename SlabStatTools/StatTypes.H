#ifndef _STATS_FORT_INTERFACE_
#define _STATS_FORT_INTERFACE_

#ifdef BL_LANG_FORT
c
c  -------------------------
c  ::: Fortran Interface :::
c  -------------------------
c
#  define   FORT_VEL_RA_RND       vel_ra_rnd
#  define   FORT_VEL_FA_RND       vel_fa_rnd
#else
//
// ---------------------------------------
// ::: C Interface to Fortran Routines :::
// ---------------------------------------
//

#if defined(BL_FORT_USE_UNDERSCORE)
#  define   FORT_VEL_RA_RND       vel_ra_rnd_
#  define   FORT_VEL_FA_RND       vel_fa_rnd_
#elif defined(BL_FORT_USE_LOWERCASE)
#  define   FORT_VEL_RA_RND       vel_ra_rnd
#  define   FORT_VEL_FA_RND       vel_fa_rnd
#elif defined(BL_FORT_USE_UPPERCASE)
#  define   FORT_VEL_RA_RND       VEL_RA_RND
#  define   FORT_VEL_FA_RND       VEL_FA_RND
#endif

//
// Use ARLIM in the following.
//

extern "C"
{
  void FORT_VEL_RA_RND(const Real *ssData, const int *nSSComp,
                         ARLIM_P(lo), ARLIM_P(hi), 
                       const int *nStats, const int *nStns,
                             Real *stats,
                             Real *physStn,
                       const int *vblo, const int *vbhi,
                       const Real *dx, 
                       const Real *probLo, const Real *probHi,
                       const int *axialDir,
                             int *nActualStations);

  void FORT_VEL_FA_RND(const Real *ssData, const int *nSSComp,
                         ARLIM_P(lo), ARLIM_P(hi), 
                       const int *nStats, const int *nStns,
                             Real *stats,
                             Real *physStn,
                       const int *vblo, const int *vbhi,
                       const Real *dx, 
                       const Real *probLo, const Real *probHi,
                       const int *axialDir,
                             int *nActualStations);


  //
  // Typedef needed for defining a pointer to the FORTRAN functions
  // used in calculating the statistics from DATABOX state variables.
  //
  typedef void (*STATS_FUNC) (const Real *ssData, const int *nSSComp,
                                ARLIM_P(lo), ARLIM_P(hi), 
                              const int *nStats, const int *nStns,
                                    Real *stats,
                                    Real *physStn,
                              const int *vblo, const int *vbhi,
                              const Real *dx, 
                              const Real *probLo, const Real *probHi,
                              const int *axialDir,
                                    int *nActualStations);
}



//
// Set Up the Structure of Available Statistics Types
//
struct StatsRecord
{
      aString    name;  // Name of statistics type
      aString    oFileBase;            // Base name for output file
      STATS_FUNC func;                 // Pointer to function that computes 
                                       //   these statistics
      int        nResStats;            // Number of resulting statistic values
};

//
// The Types
//

#if (BL_SPACEDIM == 2)
extern const struct StatsRecord StatsTypes[] =
{
    { "vel_ra_rnd",      "Vel_ra_",    FORT_VEL_RA_RND,     6},
    { "vel_fa_rnd",      "Vel_fa_",    FORT_VEL_FA_RND,     6}
};
#endif
#if (BL_SPACEDIM == 3)
extern const struct StatsRecord StatsTypes[] =
{
    { "vel_ra_rnd",      "Vel_ra_",    FORT_VEL_RA_RND,    10},
    { "vel_fa_rnd",      "Vel_fa_",    FORT_VEL_FA_RND,    10}
};
#endif

static const int NumStatsTypes = sizeof(StatsTypes)/sizeof(StatsRecord);

#endif
#endif
