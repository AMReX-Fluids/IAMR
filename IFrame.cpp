//BL_COPYRIGHT_NOTICE

#if (BL_SPACEDIM == 2)
#include "IFrame.H"
#include <BoxLib.H>

IFrame::IFrame(const Box& box, double xmin,double xmax, 
			double ymin, double ymax) : 
			xlo(xmin),xhi(xmax),ylo(ymin),yhi(ymax)		     {
	ilo = box.smallEnd()[0];
	jlo = box.smallEnd()[1];
	ihi = box.bigEnd()[0];
	jhi = box.bigEnd()[1];
	if( (ihi-ilo)*(jhi-jlo) == 0 ) BoxLib::Error("IFrame():zero width frame");
	delx = (xmax-xmin)/(ihi-ilo);
	dely = (ymax-ymin)/(jhi-jlo);
}

IntVect IFrame::toIV( const double* coord ) const {
	double x,y;
	int i,j;
	x = (coord[0]-xlo)/delx;
	if( x > 0 ) i = int(x+0.5)+ilo;
	else i = int(x-0.5)+ilo;
	
	y = (coord[1]-ylo)/dely;
	if( y > 0 ) j = int(y+0.5)+jlo;
	else j = int(y-0.5)+jlo;
	
	return IntVect(i,j);
}

double IFrame::toX(const IntVect p ) const {
	return xlo+(p[0]-ilo)*delx;
}

double IFrame::toY(const IntVect p ) const {
	return ylo+(p[1]-jlo)*dely;
}
#endif
