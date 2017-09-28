#ifndef _NS_error_F_H_
#define _NS_error_F_H_
#include <AMReX_BLFort.H>

BL_FORT_PROC_DECL(FORT_DENERROR, fort_denerror)
     (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
      const int* tagval, const int* clearval,
      amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
      const int* lo, const int* hi, const int* nvar,
      const int* domain_lo, const int* domain_hi,
      const amrex::Real* dx, const amrex::Real* xlo,
      const amrex::Real* prob_lo, const amrex::Real* time,
      const int* level);

BL_FORT_PROC_DECL(FORT_ADVERROR, fort_adverror)
     (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
      const int* tagval, const int* clearval,
      amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
      const int* lo, const int* hi, const int* nvar,
      const int* domain_lo, const int* domain_hi,
      const amrex::Real* dx, const amrex::Real* xlo,
      const amrex::Real* prob_lo, const amrex::Real* time,
      const int* level);

BL_FORT_PROC_DECL(FORT_ADV2ERROR, fort_adv2error)
     (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
      const int* tagval, const int* clearval,
      amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
      const int* lo, const int* hi, const int* nvar,
      const int* domain_lo, const int* domain_hi,
      const amrex::Real* dx, const amrex::Real* xlo,
      const amrex::Real* prob_lo, const amrex::Real* time,
      const int* level);

BL_FORT_PROC_DECL(FORT_TEMPERROR, fort_temperror)
     (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
      const int* tagval, const int* clearval,
      amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
      const int* lo, const int* hi, const int* nvar,
      const int* domain_lo, const int* domain_hi,
      const amrex::Real* dx, const amrex::Real* xlo,
      const amrex::Real* prob_lo, const amrex::Real* time,
      const int* level);

BL_FORT_PROC_DECL(FORT_MVERROR, fort_mverror)
     (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
      const int* tagval, const int* clearval,
      amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
      const int* lo, const int* hi, const int* nvar,
      const int* domain_lo, const int* domain_hi,
      const amrex::Real* dx, const amrex::Real* xlo,
      const amrex::Real* prob_lo, const amrex::Real* time,
      const int* level);

BL_FORT_PROC_DECL(FORT_STRSERROR, fort_strserror)
     (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
      const int* tagval, const int* clearval,
      amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
      const int* lo, const int* hi, const int* nvar,
      const int* domain_lo, const int* domain_hi,
      const amrex::Real* dx, const amrex::Real* xlo,
      const amrex::Real* prob_lo, const amrex::Real* time,
      const int* level);
#endif
