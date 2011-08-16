
#ifndef _ValidType_
#define _ValidType_

#include <iostream>

#include "REAL.H"
#include "aString.H"
#include "ArrayLim.H"

#include "StatTypes.H"

// ============================================================
// ===== Routine:  ValidType(aString reqType) 
// ============================================================
//
// This routine checks the type specified in reqType against those defined
// in the structure StatsTypes to make sure that reqType is valid.  It then
// returns the integer index in StatsTypes[] to the defined type corresponding 
// to reqType.  If reqType is not valid, the code is aborted.
//

int
ValidType (aString reqType)
{
  int whichType = -1;
  bool foundType = false;

  for (int nDefined = 0; nDefined < NumStatsTypes; ++nDefined) 
  {
    if (reqType == StatsTypes[nDefined].name)
    {
        foundType = true;
        whichType = nDefined;
        break;
    }
  }

  if (!foundType)
  {
    cerr << "Invalid statistics type:  " << reqType << '\n';
    BoxLib::Abort();
  }

  return whichType;
}
#endif
