
#include <PackArrayArray.H>

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
Array<Real>
PackArrayArray(int outerDim, int innerDim, Array< Array<Real> > data)
{
    Array<Real> dataRet(outerDim*innerDim);

    for (int no = 0; no < outerDim; no++)
    {
        for (int ni = 0; ni < innerDim; ni++)
            dataRet[ni*outerDim + no] = (data[no])[ni];
    }

    return dataRet;
}


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
Array< Array<Real> >
UnPackArrayArray(int outerDim, int innerDim, Array<Real> data)
{
    Array< Array<Real> > dataRet(outerDim);
    for (int no = 0; no < outerDim; no++)
      dataRet[no].resize(innerDim);

    for (int no = 0; no < outerDim; no++)
    {
        for (int ni = 0; ni < innerDim; ni++)
            (dataRet[no])[ni] = data[ni*outerDim + no];
    }

    return dataRet;
}

