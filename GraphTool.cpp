//BL_COPYRIGHT_NOTICE

#include "GraphTool.H"

#if (BL_SPACEDIM == 2)
GraphTool& GraphTool::setFrame(const IFrame &fr)
{
   frame = fr;
   return *this;
}
#endif

GTDevice GraphTool::getDevice() const 
{
   return dev;
}

#if (BL_SPACEDIM == 2)
IFrame GraphTool::getFrame() const
{
   return frame;
}
#endif

#ifdef _COMPLEX_IS_OK_
GraphTool& GraphTool::movePen(const complex &z)
{
   return movePen(real(z),imag(z));
}
#endif

#if (BL_SPACEDIM == 2)
GraphTool& GraphTool::movePen(const IntVect &v)
{
   return movePen(frame.toX(v),frame.toY(v));
}
#endif

#ifdef _COMPLEX_IS_OK_
GraphTool& GraphTool::drawLine(const complex &z)
{
   return drawLine(real(z),imag(z));
}
#endif

#if (BL_SPACEDIM == 2)
GraphTool& GraphTool::drawLine(const IntVect &v)
{
   return drawLine(frame.toX(v),frame.toY(v));
}
#endif

#ifdef _COMPLEX_IS_OK_
GraphTool& GraphTool::drawBox(const complex &zlo, const complex &zhi)
{
   return drawBox(real(zlo),imag(zlo),real(zhi),imag(zhi));
}
#endif

#if (BL_SPACEDIM == 2)
GraphTool& GraphTool::drawBox(const Box &b)
{
   double x1 = frame.toX(b.smallEnd());
   double x2 = frame.toX(b.bigEnd());
   double y1 = frame.toY(b.smallEnd());
   double y2 = frame.toY(b.bigEnd());
   return drawBox(x1,y1,x2,y2);
}
#endif

#ifdef _COMPLEX_IS_OK_
GraphTool& GraphTool::putString(const complex &z, const char *str)
{
   return putString(real(z),imag(z),str);
}
#endif

#if (BL_SPACEDIM == 2)
GraphTool& GraphTool::putString(const IntVect &v, const char *str)
{
   return putString(frame.toX(v),frame.toY(v),str);
}
#endif

#define XTOINT(x) int( width* ((x-xlo)/(xhi-xlo)) )
#define YTOINT(y) int( height*((y-ylo)/(yhi-ylo)) )
#define PS_XSCALE(x) ((x)-xlo)/(xhi-xlo)
#define PS_YSCALE(y) ((y)-ylo)/(yhi-ylo)

// private function to init a GraphTool
void GraphTool::initGT(int wid, int high, double x1, double y1,
                  double x2, double y2, const char *str, GTDevice device)
{
   name = new char[strlen(str) + 1];
   strcpy(name,str);
   dev = device;
   if ( (x1 >= x2) || (y1 >= y2) ) {
      cerr << "GraphTool: invalid domain" << endl;
      abort();
   };
   xcL = xlo = x1;     ycL = ylo = y1;
   xcU = xhi = x2;     ycU = yhi = y2;
   x0 = xlo;     y0 = ylo;
   was_out = 0;
   width = wid;
   height = high;

#if (BL_SPACEDIM == 2)
   IntVect v_lo(0,0);
   IntVect v_hi(wid-1,high-1);
   Box     bx(v_lo,v_hi);
   IFrame  fr(bx,xlo,xhi,ylo,yhi);
   frame = fr;
#endif

   win = NULL;
   ps  = NULL;
   if ((dev&xWinDevice) == xWinDevice) {
      win = new XWindow(width,height,str);
   };
   if ((dev&psDevice) == psDevice) {
      ps = new PSfile(width,height,str);
   };
   
}

#if (BL_SPACEDIM == 2)
GraphTool::GraphTool(const Box &bx,
                     const char *str, int wid, int high, GTDevice device)
{
   const int* lo_d = bx.smallEnd().getVect();
   const int* hi_d = bx.bigEnd().getVect();
   double x1 = (double) lo_d[0];
   double x2 = (double) hi_d[0];
   double y1 = (double) lo_d[1];
   double y2 = (double) hi_d[1];
   initGT(wid,high,x1,y1,x2,y2,str,device);
}
#endif

GraphTool::GraphTool(double x1, double y1, double x2, double y2,
                     const char *str, int maxwinsize,
		     GTDevice  device)
{
   double maxlen = Max(x2-x1,y2-y1);
   int    wid    = int( maxwinsize*(x2-x1)/maxlen );
   int    high   = int( maxwinsize*(y2-y1)/maxlen );
   initGT(wid,high,x1,y1,x2,y2,str,device);
}		     

GraphTool::GraphTool(double x1, double y1, double x2, double y2,
                     const char *str, int wid, int high, GTDevice device)
{
   initGT(wid,high,x1,y1,x2,y2,str,device);
}		     

GraphTool::GraphTool(const double *lo_pt, const double *hi_pt,
                     const char *str, int wid, int high, GTDevice device)
{
   initGT(wid,high,lo_pt[0],lo_pt[1],hi_pt[0],hi_pt[1],str,device);
}	  
		
#ifdef _COMPLEX_IS_OK_
GraphTool::GraphTool(const complex &zlo, const complex &zhi,
                     const char *str, int wid, int high, GTDevice device)
{
   initGT(wid,high,real(zlo),imag(zlo),real(zhi),imag(zhi),str,device);
}	  
#endif

GraphTool::GraphTool(const double *lo_pt, const double *hi_pt,
                     const char *str, int maxwinsize,
		     GTDevice  device)
{
   double x1 = lo_pt[0];
   double x2 = hi_pt[0];
   double y1 = lo_pt[1];
   double y2 = hi_pt[1];
   double maxlen = Max(x2-x1,y2-y1);
   int    wid    = int( maxwinsize*(x2-x1)/maxlen );
   int    high   = int( maxwinsize*(y2-y1)/maxlen );
   initGT(wid,high,x1,y1,x2,y2,str,device);
}

#ifdef _COMPLEX_IS_OK_
GraphTool::GraphTool(const complex &zlo, const complex &zhi,
                     const char *str, int maxwinsize,
		     GTDevice  device)
{
   double x1 = real(zlo);
   double x2 = real(zhi);
   double y1 = imag(zlo);
   double y2 = imag(zhi);
   double maxlen = Max(x2-x1,y2-y1);
   int    wid    = int( maxwinsize*(x2-x1)/maxlen );
   int    high   = int( maxwinsize*(y2-y1)/maxlen );
   initGT(wid,high,x1,y1,x2,y2,str,device);
}
#endif

#if (BL_SPACEDIM == 2)
GraphTool::GraphTool(const Box &bx,
                     const char *str, int maxwinsize,
		     GTDevice  device)
{
   const int* lo_d = bx.smallEnd().getVect();
   const int* hi_d = bx.bigEnd().getVect();
   double x1 = (double) lo_d[0];
   double x2 = (double) hi_d[0];
   double y1 = (double) lo_d[1];
   double y2 = (double) hi_d[1];
   double maxlen = Max(x2-x1,y2-y1);
   int    wid    = int( maxwinsize*(x2-x1)/maxlen );
   int    high   = int( maxwinsize*(y2-y1)/maxlen );
   initGT(wid,high,x1,y1,x2,y2,str,device);
   IFrame dummy(bx,x1,x2,y1,y2);
   frame = dummy;
}
#endif

GraphTool::~GraphTool()
{
   delete name;
   delete win;
   delete ps;
}

GraphTool& GraphTool::rmDevice(GTDevice device)
{
  if (dev & device) dev ^= device;
  return *this;
}

GraphTool& GraphTool::setDevice(GTDevice device)
{
   dev = device;
   if ( ((dev&xWinDevice) == xWinDevice) && (win == NULL) ) {
      win = new XWindow(width,height,name);
   };
   if ( ((dev&psDevice) == psDevice) && (ps == NULL) ) {
      ps = new PSfile(width,height,name);
   };
   return *this;
}

GraphTool& GraphTool::addDevice(GTDevice device)
{
   dev |= device;
   if ( ((dev&xWinDevice) == xWinDevice) && (win == NULL) ) {
      win = new XWindow(width,height,name);
   };
   if ( ((dev&psDevice) == psDevice) && (ps == NULL) ) {
      ps = new PSfile(width,height,name);
   };
   return *this;
}

GraphTool& GraphTool::newPage()
{
   if ((dev&xWinDevice) == xWinDevice) {
      win->newPage();
   };
   if ((dev&psDevice) == psDevice) {
      ps->newPage();
   };
   return *this;
}

GraphTool& GraphTool::movePen(double x, double y)
{
   int out = ( (x<xcL) || (x>xcU) || (y<ycL) || (y>ycU) );
   x0 = x;
   y0 = y;
   was_out = out;
   if (!out) {
      if ((dev&xWinDevice) == xWinDevice) {
         int xi = XTOINT(x);
         int yi = YTOINT(y);
         win->movePen(xi,yi);
      };
      if ((dev&psDevice) == psDevice) {
         ps->movePen(PS_XSCALE(x),PS_YSCALE(y));
      };
   };
   return *this;
}

#ifdef _COMPLEX_IS_OK_
GraphTool& GraphTool::setClipRegion(const complex &zlo,
                                           const complex &zhi)
{
   return setClipRegion(real(zlo),imag(zlo),real(zhi),imag(zhi));
}
#endif

#if (BL_SPACEDIM == 2)
GraphTool& GraphTool::setClipRegion(const Box &b)
{
   double x1 = frame.toX(b.smallEnd());
   double x2 = frame.toX(b.bigEnd());
   double y1 = frame.toY(b.smallEnd());
   double y2 = frame.toY(b.bigEnd());
   return setClipRegion(x1,y1,x2,y2);
}
#endif
                                    
GraphTool& GraphTool::drawLine(double x, double y, int lev)
{
   int out = ( (x<xcL) || (x>xcU) || (y<ycL) || (y>ycU) );
   if ( (!out) && (!was_out) ) {
      // old point and new point are in clip region
      if ((dev&xWinDevice) == xWinDevice) {
         int xi = XTOINT(x);
         int yi = YTOINT(y);
         win->drawLine(xi,yi);
      };
      if ((dev&psDevice) == psDevice) {
         ps->drawLine(PS_XSCALE(x),PS_YSCALE(y), lev);
      };
   } else if ((!was_out) && out) {
      // old point was in range, new point out of range
      // find exit point and only draw to there
      double c = clipRatio(x0,y0,x,y);
      double xc = x0 + c*(x-x0);
      double yc = y0 + c*(y-y0);
      if ((dev&xWinDevice) == xWinDevice) {
         int xi = XTOINT(xc);
         int yi = YTOINT(yc); 
         win->drawLine(xi,yi);
      };
      if ((dev&psDevice) == psDevice) {
         ps->drawLine(PS_XSCALE(xc),PS_YSCALE(yc), lev);
      };
      
   } else if (was_out && (!out)) {
      // old point was out of range, new point is in range
      // find entry point and only draw from there
      double c = clipRatio(x,y,x0,y0);
      double xc = x + c*(x0-x);
      double yc = y + c*(y0-y);
      if ((dev&xWinDevice) == xWinDevice) {
         int xi = XTOINT(xc);
         int yi = YTOINT(yc); 
         win->movePen(xi,yi);
         xi = XTOINT(x);
         yi = YTOINT(y); 
	 win->drawLine(xi,yi);
      };
      if ((dev&psDevice) == psDevice) {
         ps->movePen(PS_XSCALE(xc),PS_YSCALE(yc));
         ps->drawLine(PS_XSCALE(x),PS_YSCALE(y), lev);
      };
   } else {
      // both points outside range, does line intersect
      // box at all?  if so, draw on intersection
   };
   x0 = x;
   y0 = y;
   was_out = out;
   return *this;
}

GraphTool& GraphTool::drawBox(double x1, double y1,
                              double x2, double y2, int lev) 
{
   movePen(x1,y1);
   drawLine(x2,y1,lev);
   drawLine(x2,y2,lev);
   drawLine(x1,y2,lev);
   drawLine(x1,y1,lev);
   return *this;
}

GraphTool& GraphTool::setClipRegion(double x1, double y1,
                                    double x2, double y2)
{
   if ( (x1 >= x2) || (y1 >= y2) ) {
      cout << "Invalid Clip region, ignoring..." << endl;
   } else {
      xcL = x1;
      xcU = x2;
      ycL = y1;
      ycU = y2;
   };
   return *this;
}

GraphTool& GraphTool::setLineWidth(int lw)
{
   if ((dev&xWinDevice) == xWinDevice) {
      win->setLineWidth(lw);
   };
   if ((dev&psDevice) == psDevice) {
      ps->setLineWidth(lw);
   };
   return *this;
}

GraphTool& GraphTool::setFont(char *font_name)
{
   if ((dev&xWinDevice) == xWinDevice) {
      win->setFont(font_name);
   };
   if ((dev&psDevice) == psDevice) {
   };
   return *this;
}

GraphTool& GraphTool::putString(double x, double y, const char *str)
{
   if ((dev&xWinDevice) == xWinDevice) {
         int xi = XTOINT(x);
         int yi = YTOINT(y);
         win->putString(xi,yi,str);
   };
   if ((dev&psDevice) == psDevice) {
   };
   return *this;
}

GraphTool& GraphTool::setfgColor(const char* color)
{
   if ((dev&xWinDevice) == xWinDevice) {
      win->setfgColor(color);
   };
   if ((dev&psDevice) == psDevice) {
   };
   return *this;
}

GraphTool& GraphTool::setbgColor(const char* color)
{
   if ((dev&xWinDevice) == xWinDevice) {
      win->setbgColor(color);
   };
   if ((dev&psDevice) == psDevice) {
   };
   return *this;
}

GraphTool& GraphTool::setfgColor(int color)
{
   if ((dev&xWinDevice) == xWinDevice) {
      win->setfgColor(color);
   };
   if ((dev&psDevice) == psDevice) {
   };
   return *this;
}

GraphTool& GraphTool::setbgColor(int color)
{
   if ((dev&xWinDevice) == xWinDevice) {
      win->setbgColor(color);
   };
   if ((dev&psDevice) == psDevice) {
   };
   return *this;
}

GraphTool& GraphTool::defineCmap(unsigned short *red, unsigned short *green,
                                 unsigned short *blue, int num)
{
   if ((dev&xWinDevice) == xWinDevice) {
      win->defineCmap(red, green, blue, num);
   };
   if ((dev&psDevice) == psDevice) {
   };
   return *this;
}

#ifdef _COMPLEX_IS_OK_
int GraphTool::getMouse(complex &z) const
{
   double x,y;
   int but = getMouse(x,y);
   complex c(x,y);
   z = c;
   return but;
}
#endif

#if (BL_SPACEDIM == 2)
int GraphTool::getMouse(IntVect &v) const
{
   double x[2];
   int but = getMouse(x[0],x[1]);
   v = frame.toIV(x);
   return but;
}
#endif

int GraphTool::getMouse(double &x, double &y) const
{
   int button = -1;
   if ((dev&xWinDevice) == xWinDevice) {
      int i,j;
      button = win->getMouse(i,j);
      x = xlo + (xhi-xlo)*( double(i)/double(width) );
      y = ylo + (yhi-ylo)*( double(j)/double(height) );
   };
   return button;
}

// contour plotting
int GraphTool::contour(const double *data, double value,
                       int has_mask, const int *mask,
                       int nx, int /* ny */, int mx, int my,
                       double xlft, double ybot, double xrgt, double ytop)
{
// data     = of data to be contoured
// value    = value to contour
// has_mask = true if mask array available
// mask     = array of mask values.  will not contour in masked off cells
// nx       = dimension of arrays in X direction
// ny       = dimension of arrays in Y direction
// mx       = number of cells in X direction to plot
// my       = number of cells in Y direction to plot
// xlft     = position of left   edge of grid in domain
// xrgt     = position of right  edge of grid in domain
// ybot     = position of bottom edge of grid in domain
// ytop     = position of top    edge of grid in domain
#define VAL(i,j) data[(i)+(j)*nx]
#define MSK(i,j) mask[(i)+(j)*nx]
#define BTWN(a,b,c) ( ((a<=b)&&(b<=c)) || ((a>=b)&&(b>=c)) )
//   double dx = (xrgt-xlft)/(nx-1);
//   double dy = (ytop-ybot)/(ny-1);
   double dx = (xrgt-xlft)/(mx-1);
   double dy = (ytop-ybot)/(my-1);
   bool lft, rgt, bot, top;   // does contour line intersect this side?
   double  xl, yl;               // where contour line intersects lft side
   double  xr, yr;               // where contour line intersects rgt side
   double  xb, yb;               // where contour line intersects bot side
   double  xt, yt;               // where contour line intersects top side
   bool failure_status = false;
   for (int j = 0; j < my-1; j++)
   for (int i = 0; i < mx-1; i++) {
      if (has_mask) {
         int m = MSK(i,j)+MSK(i+1,j)+MSK(i+1,j+1)+MSK(i,j+1);
         if (m > 0) continue;
      };
      double  lb = VAL(i,j);            // left bottom value
      double  lt = VAL(i,j+1);          // left top value
      double  rb = VAL(i+1,j);          // right bottom value
      double  rt = VAL(i+1,j+1);        // right top value
      xl = xlft + dx*(i);
      xr = xl + dx;
      yb = ybot + dy*(j);
      yt = yb + dy;

      // figure out where things intersect the cell
      if (lft = BTWN(lb,value,lt)) {
         if (lb != lt) {
            yl = yb + dy*(value-lb)/(lt-lb);
         } else {
            yl = yb;
            failure_status = true;
         };
      };
      if (rgt = BTWN(rb,value,rt)) {
         if (rb != rt) {
            yr = yb + dy*(value-rb)/(rt-rb);
         } else {
            yr = yb;
            failure_status = true;
         };
      };
      if (bot = BTWN(lb,value,rb)) {
         if (lb != rb) {
            xb = xl + dx*(value-lb)/(rb-lb);
         } else {
            xb = xr;
            failure_status = true;
         };
      };
      if (top = BTWN(lt,value,rt)) {
         if (lt != rt) {
            xt = xl + dx*(value-lt)/(rt-lt);
         } else {
            xt = xr;
            failure_status = true;
         };
      };

      // finally, draw contour line
      if (lft && rgt && bot && top) {
         // intersects all sides, generate saddle point
         movePen(xl,yl);
         drawLine(xr,yr);
         movePen(xt,yt);
         drawLine(xb,yb);
      } else if (top && bot) {
         // only intersects top and bottom sides
         movePen(xt,yt);
         drawLine(xb,yb);
      } else if (lft) {
         movePen(xl,yl);
         if (rgt) {
            drawLine(xr,yr);
         } else if (top) {
            drawLine(xt,yt);
         } else {
            drawLine(xb,yb);
         };
      } else if (rgt) {
         movePen(xr,yr);
         if (top) {
            drawLine(xt,yt);
         } else {
            drawLine(xb,yb);
         };
      };

   }; // for I,J

   return failure_status;
#undef VAL
#undef MSK
#undef BTWN
}

#if (BL_SPACEDIM == 2)
// contour plotting
int GraphTool::contour(const double *data, double value,
                       int has_mask, const int *mask,
	               const Box &dims, const Box &subrange,
		       double xlft, double ybot, double xrgt, double ytop)
{
// same as above except for:
// dims     = dimension of arrays (ilo,ihi,jlo,jhi)
// subrange = subrange of array to plot (irlo,irhi,jrlo,jrhi)
   const int *lo_n = dims.smallEnd().getVect();
   const int *hi_n = dims.bigEnd().getVect();
   const int ilo = lo_n[0];
   const int ihi = hi_n[0];
   const int jlo = lo_n[1];
   const int jhi = hi_n[1];
   const int nx = ihi - ilo + 1;
   const int ny = jhi - jlo + 1;
   const int *lo_r = subrange.smallEnd().getVect();
   const int *hi_r = subrange.bigEnd().getVect();
   const int irlo = lo_r[0];
   const int irhi = hi_r[0];
   const int jrlo = lo_r[1];
   const int jrhi = hi_r[1];
   const int mx   = irhi - irlo + 1;
   const int my   = jrhi - jrlo + 1;
   const double *dstart = data + (irlo-ilo) + nx*(jrlo-jlo);
   const int    *mstart = mask + (irlo-ilo) + nx*(jrlo-jlo);
   return contour(dstart,value,has_mask,mstart,nx,ny,mx,my,
                  xlft,ybot,xrgt,ytop);
}
#endif
#if (BL_SPACEDIM == 2)
int GraphTool::contour(const double* data, double value,
               int has_mask, const int *mask,
	       const Box &dims, const Box &subrange,
	       const Box &position)
{
   const int *lo_n = dims.smallEnd().getVect();
   const int *hi_n = dims.bigEnd().getVect();
   const int ilo = lo_n[0];
   const int ihi = hi_n[0];
   const int jlo = lo_n[1];
   const int jhi = hi_n[1];
   const int nx = ihi - ilo + 1;
   const int ny = jhi - jlo + 1;
   const int *lo_r = subrange.smallEnd().getVect();
   const int *hi_r = subrange.bigEnd().getVect();
   const int irlo = lo_r[0];
   const int irhi = hi_r[0];
   const int jrlo = lo_r[1];
   const int jrhi = hi_r[1];
   const int mx   = irhi - irlo + 1;
   const int my   = jrhi - jrlo + 1;
   double xlft = frame.toX(position.smallEnd());
   double xrgt = frame.toX(position.bigEnd());
   double ybot = frame.toY(position.smallEnd());
   double ytop = frame.toY(position.bigEnd());
   const double  *dstart = data + (irlo-ilo) + nx*(jrlo-jlo);
   const int     *mstart = mask + (irlo-ilo) + nx*(jrlo-jlo);
   return contour(dstart,value,has_mask,mstart,nx,ny,mx,my,
                  xlft,ybot,xrgt,ytop);
}
#endif

// contour plotting.  SINGLE PRECISION
int GraphTool::contour(const float *data, double value,
                        int has_mask, const int *mask,
                        int nx, int ny, int mx, int my,
                        double xlft, double ybot, double xrgt, double ytop)
{
// data     = of data to be contoured
// value    = value to contour
// has_mask = true if mask array available
// mask     = array of mask values.  will not contour in masked off cells
// nx       = dimension of arrays in X direction
// ny       = dimension of arrays in Y direction
// mx       = number of cells in X direction to plot
// my       = number of cells in Y direction to plot
// xlft     = position of left   edge of grid in domain
// xrgt     = position of right  edge of grid in domain
// ybot     = position of bottom edge of grid in domain
// ytop     = position of top    edge of grid in domain
#define VAL(i,j) data[(i)+(j)*nx]
#define MSK(i,j) mask[(i)+(j)*nx]
#define BTWN(a,b,c) ( ((a<=b)&&(b<=c)) || ((a>=b)&&(b>=c)) )
   double dx = (xrgt-xlft)/(nx-1);
   double dy = (ytop-ybot)/(ny-1);
   bool lft, rgt, bot, top;   // does contour line intersect this side?
   double  xl, yl;               // where contour line intersects lft side
   double  xr, yr;               // where contour line intersects rgt side
   double  xb, yb;               // where contour line intersects bot side
   double  xt, yt;               // where contour line intersects top side
   bool failure_status = false;
   for (int j = 0; j < my-1; j++)
   for (int i = 0; i < mx-1; i++) {
      if (has_mask) {
         int m = MSK(i,j)+MSK(i+1,j)+MSK(i+1,j+1)+MSK(i,j+1);
         if (m > 0) continue;
      };
      double  lb = VAL(i,j);            // left bottom value
      double  lt = VAL(i,j+1);          // left top value
      double  rb = VAL(i+1,j);          // right bottom value
      double  rt = VAL(i+1,j+1);        // right top value
      xl = xlft + dx*(i);
      xr = xl + dx;
      yb = ybot + dy*(j);
      yt = yb + dy;

      // figure out where things intersect the cell
      if (lft = BTWN(lb,value,lt)) {
         if (lb != lt) {
            yl = yb + dy*(value-lb)/(lt-lb);
         } else {
            yl = yb;
            failure_status = true;
         };
      };
      if (rgt = BTWN(rb,value,rt)) {
         if (rb != rt) {
            yr = yb + dy*(value-rb)/(rt-rb);
         } else {
            yr = yb;
            failure_status = true;
         };
      };
      if (bot = BTWN(lb,value,rb)) {
         if (lb != rb) {
            xb = xl + dx*(value-lb)/(rb-lb);
         } else {
            xb = xr;
            failure_status = true;
         };
      };
      if (top = BTWN(lt,value,rt)) {
         if (lt != rt) {
            xt = xl + dx*(value-lt)/(rt-lt);
         } else {
            xt = xr;
            failure_status = true;
         };
      };

      // finally, draw contour line
      if (lft && rgt && bot && top) {
         // intersects all sides, generate saddle point
         movePen(xl,yl);
         drawLine(xr,yr);
         movePen(xt,yt);
         drawLine(xb,yb);
      } else if (top && bot) {
         // only intersects top and bottom sides
         movePen(xt,yt);
         drawLine(xb,yb);
      } else if (lft) {
         movePen(xl,yl);
         if (rgt) {
            drawLine(xr,yr);
         } else if (top) {
            drawLine(xt,yt);
         } else {
            drawLine(xb,yb);
         };
      } else if (rgt) {
         movePen(xr,yr);
         if (top) {
            drawLine(xt,yt);
         } else {
            drawLine(xb,yb);
         };
      };

   }; // for I,J

   return failure_status;
#undef VAL
#undef MSK
#undef BTWN
}

#if (BL_SPACEDIM == 2)
// contour plotting
int GraphTool::contour(const float *data, double value,
                       int has_mask, const int *mask,
	               const Box &dims, const Box &subrange,
		       double xlft, double ybot, double xrgt, double ytop)
{
// same as above except for:
// dims     = dimension of arrays (ilo,ihi,jlo,jhi)
// subrange = subrange of array to plot (irlo,irhi,jrlo,jrhi)
   const int *lo_n = dims.smallEnd().getVect();
   const int *hi_n = dims.bigEnd().getVect();
   const int ilo = lo_n[0];
   const int ihi = hi_n[0];
   const int jlo = lo_n[1];
   const int jhi = hi_n[1];
   const int nx = ihi - ilo + 1;
   const int ny = jhi - jlo + 1;
   const int *lo_r = subrange.smallEnd().getVect();
   const int *hi_r = subrange.bigEnd().getVect();
   const int irlo = lo_r[0];
   const int irhi = hi_r[0];
   const int jrlo = lo_r[1];
   const int jrhi = hi_r[1];
   const int mx   = irhi - irlo + 1;
   const int my   = jrhi - jrlo + 1;
   const float  *dstart = data + (irlo-ilo) + nx*(jrlo-jlo);
   const int    *mstart = mask + (irlo-ilo) + nx*(jrlo-jlo);
   return contour(dstart,value,has_mask,mstart,nx,ny,mx,my,
                  xlft,ybot,xrgt,ytop);
}
#endif
#if (BL_SPACEDIM == 2)
int GraphTool::contour(const float* data, double value,
               int has_mask, const int *mask,
	       const Box &dims, const Box &subrange,
	       const Box &position)
{
   const int *lo_n = dims.smallEnd().getVect();
   const int *hi_n = dims.bigEnd().getVect();
   const int ilo = lo_n[0];
   const int ihi = hi_n[0];
   const int jlo = lo_n[1];
   const int jhi = hi_n[1];
   const int nx = ihi - ilo + 1;
   const int ny = jhi - jlo + 1;
   const int *lo_r = subrange.smallEnd().getVect();
   const int *hi_r = subrange.bigEnd().getVect();
   const int irlo = lo_r[0];
   const int irhi = hi_r[0];
   const int jrlo = lo_r[1];
   const int jrhi = hi_r[1];
   const int mx   = irhi - irlo + 1;
   const int my   = jrhi - jrlo + 1;
   double xlft = frame.toX(position.smallEnd());
   double xrgt = frame.toX(position.bigEnd());
   double ybot = frame.toY(position.smallEnd());
   double ytop = frame.toY(position.bigEnd());
   const float *dstart = data + (irlo-ilo) + nx*(jrlo-jlo);
   const int   *mstart = mask + (irlo-ilo) + nx*(jrlo-jlo);
   return contour(dstart,value,has_mask,mstart,nx,ny,mx,my,
                  xlft,ybot,xrgt,ytop);
}
#endif

//  Private function: Determines fraction of distance from
//  (x1,y1) to (x2,y2) at which boundary is intersected
double GraphTool::clipRatio(double x1, double y1, double x2, double y2)
{
   double r = 1.0;
   if (x2<xcL) r = Min(r,(x1-xcL)/(x1-x2));
   if (y2<ycL) r = Min(r,(y1-ycL)/(y1-y2));
   if (x2>xcU) r = Min(r,(xcU-x1)/(x2-x1));
   if (y2>ycU) r = Min(r,(ycU-y1)/(y2-y1));
   return 0.9999*r;
}

#undef  XTOINT
#undef  YTOINT
#undef  PS_XSCALE
#undef  PS_YSCALE

