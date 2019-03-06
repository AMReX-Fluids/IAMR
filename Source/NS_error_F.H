#ifndef _NS_error_F_H_
#define _NS_error_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
  extern "C" {
#endif


void FORT_DENERROR
     (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
      const int* tagval, const int* clearval,
      amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
      const int* lo, const int* hi, const int* nvar,
      const int* domain_lo, const int* domain_hi,
      const amrex::Real* dx, const amrex::Real* xlo,
      const amrex::Real* prob_lo, const amrex::Real* time,
      const int* level);

void FORT_ADVERROR
     (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
      const int* tagval, const int* clearval,
      amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
      const int* lo, const int* hi, const int* nvar,
      const int* domain_lo, const int* domain_hi,
      const amrex::Real* dx, const amrex::Real* xlo,
      const amrex::Real* prob_lo, const amrex::Real* time,
      const int* level);

void FORT_ADV2ERROR
     (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
      const int* tagval, const int* clearval,
      amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
      const int* lo, const int* hi, const int* nvar,
      const int* domain_lo, const int* domain_hi,
      const amrex::Real* dx, const amrex::Real* xlo,
      const amrex::Real* prob_lo, const amrex::Real* time,
      const int* level);

void FORT_TEMPERROR
     (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
      const int* tagval, const int* clearval,
      amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
      const int* lo, const int* hi, const int* nvar,
      const int* domain_lo, const int* domain_hi,
      const amrex::Real* dx, const amrex::Real* xlo,
      const amrex::Real* prob_lo, const amrex::Real* time,
      const int* level);

void FORT_MVERROR
     (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
      const int* tagval, const int* clearval,
      amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
      const int* lo, const int* hi, const int* nvar,
      const int* domain_lo, const int* domain_hi,
      const amrex::Real* dx, const amrex::Real* xlo,
      const amrex::Real* prob_lo, const amrex::Real* time,
      const int* level);
#ifdef __cplusplus
  }
#endif
#endif


