
#ifndef _PackArrayArray__
#define _PackArrayArray_

//
// Includes required for the interface definition
//
#include <Array.H>
#include <REAL.H>


// ================================================================
// ===== Routine:  PackArrayArray(int outerDim, int innerDim, 
//                                Array< Array<Real> > data)
// ================================================================
//
// This routine pack an Array<Array<Real>> into an Array<Real> so it
// can be passed into Fortran.  The ordering of the pack is such that the 
// inner index in the Array<Array<>> is the first index in fortran
// (the slowest varying).
//
Array<Real>  PackArrayArray(int outerDim, 
                            int innerDim, 
                            Array< Array<Real> > data);


// ================================================================
// ===== Routine:  UnPackArrayArray(int outerDim, int innerDim, 
//                                  Array<Real> data)
// ================================================================
//
// This routine unpacks an Array<Real> into an Array<Array<Real>> so it
// can be passed into Fortran.  The ordering of the pack is such that the 
// inner index in the Array<Array<>> is the first index in fortran
// (the slowest varying).  This is the unpack to go with the pack above.
//
Array< Array<Real> >  UnPackArrayArray(int outerDim, 
                                       int innerDim, 
                                       Array<Real> data);
#endif

