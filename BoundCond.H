#ifndef _BOUNDCOND_H_
#define _BOUNDCOND_H_

//
// $Id: BoundCond.H,v 1.3 1997-09-26 16:56:53 lijewski Exp $
//

#include <Assert.H>

//@Man:
/*@Memo:
        Maintain an identifier for boundary condition types.
*/	
/*@Doc:
        This is a placeholder for more extensive boundary condition
	implementations, which might include stencils, etc.  Presently,
	boundary conditions are specified via an integer identifier.
	This class maintains that integer.
*/	
class BoundCond
{
public:
//@ManMemo: administrative functions
//@ManDoc: constructor, initializes boundary condition type to default
    BoundCond() : bctype(-1){}
//@ManDoc: constructor, initializes boundary condition type to specified value
    BoundCond(int _bctype) : bctype(_bctype){}

//@ManMemo: access functions
//@ManDoc: return index of boundary condition type
    operator int () const
    {
	assert(bctype != -1);
	return bctype;
    }
//@ManDoc: set index of boundary condition type
    BoundCond& operator = (const int _bctype)
    {
	bctype = _bctype;
	return *this;
    }
private:
// Index of boundary condition type
    int bctype;
};

#endif
