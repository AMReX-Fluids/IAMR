
#include <aString.H>
#include <Misc.H>
#include <Utility.H>
#include <RealBox.H>

REAL REALBOX::eps = 1.0e-6;

// -----------------------------------------------------------
REALBOX::REALBOX()
{
    D_TERM(xlo[0] , = xlo[1] , = xlo[2] ) = 0.;
    D_TERM(xhi[0] , = xhi[1] , = xhi[2] ) = -1.;
    computeBoxLen();
}

// -----------------------------------------------------------
REALBOX::REALBOX(const REAL* lo, const REAL* hi)
{
    D_EXPR(xlo[0] = lo[0] , xlo[1] = lo[1] , xlo[2] = lo[2]);
    D_EXPR(xhi[0] = hi[0] , xhi[1] = hi[1] , xhi[2] = hi[2]);
    computeBoxLen() ;
}

// -----------------------------------------------------------
REALBOX::REALBOX(const BOX& bx, const REAL* dx, const REAL* base)
{
    const int* lo = bx.loVect();
    const int* hi = bx.hiVect();
    int i;
    for (i = 0; i < BL_SPACEDIM; i++) {
	xlo[i] = base[i] + dx[i]*lo[i];
	int shft = (bx.type(i)==IndexType::CELL ? 1 : 0);
	xhi[i] = base[i] + dx[i]*(hi[i]+ shft);
    }   
    computeBoxLen() ;
}

// -----------------------------------------------------------
REALBOX::REALBOX(D_DECL(REAL x0, REAL y0, REAL z0),
		 D_DECL(REAL x1, REAL y1, REAL z1))
{
    D_EXPR(xlo[0] = x0 , xlo[1] = y0 , xlo[2] = z0);
    D_EXPR(xhi[0] = x1 , xhi[1] = y1 , xhi[2] = z1);
    computeBoxLen() ;
}

// -----------------------------------------------------------
void REALBOX::setLo(const REAL* lo)
{
    D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]);
    computeBoxLen();
}

void REALBOX::setLo(const Array<REAL> &lo)
{
    D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]);
    computeBoxLen();
}

// -----------------------------------------------------------
void REALBOX::setHi(const REAL* hi)
{
    D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]);
    computeBoxLen();
}

void REALBOX::setHi(const Array<REAL>& hi)
{
    D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]);
    computeBoxLen();
}

// -----------------------------------------------------------
void
REALBOX::setLo(int indx, REAL lo)
{
   assert( indx >= 0 && indx < BL_SPACEDIM);
   xlo[indx] = lo;
   computeBoxLen();
}

// -----------------------------------------------------------
void
REALBOX::setHi(int indx, REAL hi)
{
    assert( indx >= 0 && indx < BL_SPACEDIM);
    xhi[indx] = hi;
    computeBoxLen();
}

// -----------------------------------------------------------
int
REALBOX::contains(const REAL* point)
{
    return  (xlo[0]-eps < point[0]) && (point[0] < xhi[0]+eps)
#if (BL_SPACEDIM > 1)   
        && (xlo[1]-eps < point[1]) && (point[1] < xhi[1]+eps)
#endif
#if (BL_SPACEDIM > 2)   
        && (xlo[2]-eps < point[2]) && (point[2] < xhi[2]+eps)
#endif
   ;
}

// -----------------------------------------------------------
int
REALBOX::contains(const REALBOX& rb)
{
    return (contains(rb.xlo) && contains(rb.xhi));
}

// -----------------------------------------------------------
int
REALBOX::ok() const
{
    return (len[0] > eps)
#if (BL_SPACEDIM > 1)
	&& (len[1] > eps)
#endif   
#if (BL_SPACEDIM > 2)
	&& (len[2] > eps)
#endif
   ;
}

#if 0
// -----------------------------------------------------------
int REALBOX::intersects(const REALBOX& /*rb*/) 
{
   // can realboxes intersect in zero volume areas?
   cerr << "REALBOX::intersects not implemented" << endl;
   mpAbort();
   return 0;
}

// -----------------------------------------------------------
REALBOX& REALBOX::operator &= (const REALBOX& /*bx*/)
{
   cerr << "REALBOX::op &= not implemented" << endl;
   mpAbort();
   return *this;
}

// -----------------------------------------------------------
REALBOX
REALBOX::operator & (const REALBOX& /*bx*/)
{
   cerr << "REALBOX::op & not implemented" << endl;
   mpAbort();
   return *this;
}
#endif 

// -----------------------------------------------------------
ostream&
operator << (ostream &os, const REALBOX& b)
{
    os << "(REALBOX ";
    int i;
    for (i = 0; i < BL_SPACEDIM; i++) {
        os << b.xlo[i] << ' ' << b.xhi[i] << ' ';
    }
    os << ')';

    return os;
}

// -----------------------------------------------------------
istream&
operator >> (istream &is, REALBOX& b)
{
    is.ignore(BL_IGNORE_MAX,'(');
    aString s;
    is >> s;
    if(s != "REALBOX") {
        cerr << "unexpected token in RealBox: " << s ;
        abort();
    }
    int i;
    for (i = 0; i < BL_SPACEDIM; i++) {
        is >> b.xlo[i] >> b.xhi[i];
    }
    is.ignore(BL_IGNORE_MAX, ')');
    b.computeBoxLen();

    return is;
}
