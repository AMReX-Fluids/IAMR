#ifndef _RealBox_H_
#define _RealBox_H_ 

//
// $Id: RealBox.H,v 1.4 1997-10-01 01:03:14 car Exp $
//

#include <Array.H>
#include <REAL.H>
#include <Box.H>

// rules
// ok  if and only if (volume > 0)
// is zero volume intersection valid?

class RealBox
{
protected:
    REAL xlo[BL_SPACEDIM], xhi[BL_SPACEDIM];
    REAL len[BL_SPACEDIM];
    static REAL eps;

public:
#ifdef __GNUG__
     bool operator== (const RealBox& rhs) const { return true; }
#endif 
      // -------------------------------------------
      // ---------  I/O friend functions  ----------
      // -------------------------------------------
    friend ostream& operator << (ostream&, const RealBox& );
    friend istream& operator >> (istream&, RealBox& );

      // -------------------------------------------
      // ------------  constructors  ---------------
      // -------------------------------------------

      // construct invalid region
    RealBox();

      // construct region given diagional points
    RealBox(const REAL *lo, const REAL* hi);
      
      // construct region given index box, cell spacing
      // and physical location of index (0,0,0)
    RealBox(const BOX& bx, const REAL* dx, const REAL* base);

      // read region from input file
    RealBox(istream &is);

      // explicit dimension specific constructors
    RealBox(D_DECL(REAL x0, REAL y0, REAL z0),
	    D_DECL(REAL x1, REAL y1, REAL z1));

      // -------------------------------------------
      // -----------  static functions  ------------
      // -------------------------------------------
      // set epsilon for determining intersections
    static void setEpsilon(REAL epsilon) { eps = epsilon; }
    static REAL epsilon() { return eps; }

      // -------------------------------------------
      // ------------  access functions  -----------
      // -------------------------------------------
    const REAL* lo() const {return xlo;}
    const REAL* hi() const {return xhi;}
    const REAL* length() const {return len;}
    REAL  lo(int dir) const {return xlo[dir];}
    REAL  hi(int dir) const {return xhi[dir];}
    REAL  length(int dir) const {return len[dir];}

      // -------------------------------------------
      // ----------  Modify functions  -------------
      // -------------------------------------------
    void setLo(const REAL* lo);
    void setLo(const Array<REAL> &lo);
    void setHi(const REAL* hi);
    void setHi(const Array<REAL> &hi);
    void setLo(int indx, REAL lo);
    void setHi(int indx, REAL hi);

      // -------------------------------------------
      // ----------  Boolean functions  ------------
      // -------------------------------------------
    int ok() const;
    int contains(const REAL* point);
    int contains(const RealBox &bx);

      // -------------------------------------------
      // ---------  conversion functions  ----------
      // -------------------------------------------
    // BOX  CCBox(const REAL *dx, const REAL *base);
    // BOX  NodeBox(const REAL *dx, const REAL* base);

protected:
    void computeBoxLen() 
    {
	D_EXPR(len[0] = xhi[0]-xlo[0],
	       len[1] = xhi[1]-xlo[1],
	       len[2] = xhi[2]-xlo[2]);
    }
};

#endif
