//BL_COPYRIGHT_NOTICE

#ifndef _IFRAME_H_
#define _IFRAME_H_	"@(#)IFrame.H	1.2\t10/5/92"

#if BL_SPACEDIM==2
#include "Box.H"

// describes linear transformation between integer coords and real coords
// in 2d.  useful for IGraphTool
class IFrame {
	private:
		int ilo,ihi,jlo,jhi;
		double xlo,xhi,ylo,yhi;
		double delx,dely;
	public:
		IFrame(){}
		IFrame(const Box&,double xmin,double xmax, 
				double ymin, double ymax);
		~IFrame(){}
		
		IntVect toIV( const double* coord) const;
		double toX(const IntVect) const;
		double toY(const IntVect) const;
};
#endif
#endif
