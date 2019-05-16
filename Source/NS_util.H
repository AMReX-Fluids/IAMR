#ifndef NS_util_H 
#define NS_util_H

#include <AMReX_MultiFab.H>

namespace amrex
{
  amrex::Vector<amrex::Real>
  VectorMax(const amrex::Vector<const amrex::MultiFab *>& mf,
            const amrex::IntVect&                         tilesize,
            int                                           sComp,
            int                                           nComp,
            int                                           nGrow);

  amrex::Vector<amrex::Real>
  VectorMin(const amrex::Vector<const amrex::MultiFab *>& mf,
            const amrex::IntVect&                         tilesize,
            int                                           sComp,
            int                                           nComp,
            int                                           nGrow);

  amrex::Vector<amrex::Real>
  VectorMaxAbs(const amrex::Vector<const amrex::MultiFab *>& mf,
               const amrex::IntVect&                         tilesize,
               int                                           sComp,
               int                                           nComp,
               int                                           nGrow);
}

#endif
