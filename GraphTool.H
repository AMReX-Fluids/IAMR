//BL_COPYRIGHT_NOTICE

#ifndef _GraphTool_H_
#define _GraphTool_H_  "@(#)GraphTool.H	1.3\t1/8/93"

#include <iostream.h>
#include <string.h>
#include <math.h>
#ifdef _COMPLEX_IS_OK_
#include <complex.h>
#endif
#include <Misc.H>
#include "Boolean.H"
#include "XWindow.H"
#include "PSfile.H"

#if (BL_SPACEDIM == 2)
#include "IntVect.H"
#include "Box.H"
#include "IFrame.H"
#endif

typedef unsigned short GTDevice;
const GTDevice noDevice = 0;
const GTDevice xWinDevice = 1;
const GTDevice psDevice = 2;

// This class is a C++ version of Louis Howell's GraphTools library.
// It has the ability to do simple line drawing graphics, including
// contour plotting.
// The devices supported are X-windows and postscript files and 
// an object of this class can have either one or both active at the
// same time.
// Many of the operations, such as querring the mouse and setting color
// and fonts are not implemented for the postscript device and
// become no-ops.
// The user can change the currently active devices at any time
// through a variety of member functions and all subsequent
// graphics operations will only be directed to those devices.
// One sets the active devices by ORing together any combination
// of the above constants (noDevice, XWinDevice, psDevice).
class GraphTool {
   public:
      GraphTool(double x1, double y1, double x2, double y2, 
                const char *str, int wid, int high,
		GTDevice device = xWinDevice);
      // create a graphtool:
      //   (x1,y1) X (x2,y2) physical coordinates of window
      //   str               window/postscript file name
      //   (wid,high)        width and height (in pixels) of the window
      //   device            selected devices           

      GraphTool(const double *lo_pt, const double *hi_pt,
                const char *str, int wid, int high,
		GTDevice device = xWinDevice);
      // create a graphtool:
      //   lo_pt, hi_pt      arrays specifying physical coordinates of window
      //   str               window/postscript file name
      //   (wid,high)        width and height (in pixels) of the window
      //   device            selected devices           

#ifdef _COMPLEX_IS_OK_
      GraphTool(const complex &zlo, const complex &zhi,
                const char *str, int wid, int high,
		GTDevice device = xWinDevice);
      // create a graphtool:
      //   zlo, zhi	     physical coords specified in complex form
      //   str               window/postscript file name
      //   (wid,high)        width and height (in pixels) of the window
      //   device            selected devices           
#endif

#if (BL_SPACEDIM == 2)
      GraphTool(const Box &bx,
                const char *str, int wid, int high,
		GTDevice device = xWinDevice);
      // create a graphtool:
      //   bx	             physical coordinates specified via box
      //   str               window/postscript file name
      //   (wid,high)        width and height (in pixels) of the window
      //   device            selected devices           
#endif

      GraphTool(double x1, double y1, double x2, double y2,
                const char *str, int maxwinsize = 600,
		GTDevice  device = xWinDevice);
      // greate a graphtool with correct aspec ratio and
      // maximum length having maxwinsize pixels.

      GraphTool(const double *lo_pt, const double *hi_pt,
                const char *str, int maxwinsize = 600,
		GTDevice  device = xWinDevice);
      // greate a graphtool with correct aspec ratio and
      // maximum length having maxwinsize pixels.

#ifdef _COMPLEX_IS_OK_
      GraphTool(const complex &zlo, const complex &zhi,
                const char *str, int maxwinsize = 600,
		GTDevice  device = xWinDevice);
      // greate a graphtool with correct aspec ratio and
      // maximum length having maxwinsize pixels.
#endif

#if (BL_SPACEDIM == 2)
      GraphTool(const Box &bx,
                const char *str, int maxwinsize = 600,
		GTDevice  device = xWinDevice);
      // greate a graphtool with correct aspec ratio and
      // maximum length having maxwinsize pixels.
#endif

      ~GraphTool();
      // distructor

      GraphTool& setDevice(GTDevice device);
      // set these as the active devices

#if (BL_SPACEDIM == 2)
      GraphTool& setFrame(const IFrame &fr);
      // set the current int to double mapping
#endif

      GTDevice getDevice() const;
      // returns currently active devices

#if (BL_SPACEDIM == 2)
      IFrame getFrame() const;
      // returns currently active int to double mapping
#endif

      GraphTool& addDevice(GTDevice device);
      // add these devices to the currently active devices

      GraphTool& rmDevice(GTDevice device);
      // remove these devices from the currently active devices

      GraphTool& newPage();
      // clear XWindow if active
      // start new page of PostScript file if active

      GraphTool& movePen(double x, double y);
      // move to given location without drawing line
      
#ifdef _COMPLEX_IS_OK_
      GraphTool& movePen(const complex &z);
      // move to given location without drawing line
#endif

#if (BL_SPACEDIM == 2)
      GraphTool& movePen(const IntVect &v);
      // move to given location without drawing line
#endif

      GraphTool& drawLine(double x, double y, int lev = 0);
      // draw line from current location to this location

#ifdef _COMPLEX_IS_OK_
      GraphTool& drawLine(const complex &z);
      // draw line from current location to this location
#endif

#if (BL_SPACEDIM == 2)
      GraphTool& drawLine(const IntVect &v);
      // draw line from current location to this location
#endif

      GraphTool& drawBox(double xlo, double ylo, double xhi, double yhi,
		int lev=0);
      // draw a rectangle with the given diagional locations

#ifdef _COMPLEX_IS_OK_
      GraphTool& drawBox(const complex &zlo, const complex &zhi);
      // draw a rectangle with the given diagional locations
#endif

#if (BL_SPACEDIM == 2)
      GraphTool& drawBox(const Box &bx);
      // draw a rectangle with the given diagional locations
#endif

      GraphTool& setClipRegion(double xlo, double ylo,
                               double xhi, double yhi);
      // set the clipping region so that nothing is drawn outside
      // this this box.

#ifdef _COMPLEX_IS_OK_
      GraphTool& setClipRegion(const complex &zlo, const complex &zhi);
      // set the clipping region so that nothing is drawn outside
      // this this box.
#endif

#if (BL_SPACEDIM == 2)
      GraphTool& setClipRegion(const Box &bx);
      // set the clipping region so that nothing is drawn outside
      // this this box.
#endif

      GraphTool& setLineWidth(int lw);
      // change the width of the lines drawn.
      
      GraphTool& setFont(char *);
      // change the font to the given XWindow font.
      // NoOp for PostScript device.
      
      GraphTool& putString(double x, double y, const char *str);
      // output the given string at the given location
      // NoOp for PostScript device.

#ifdef _COMPLEX_IS_OK_
      GraphTool& putString(const complex &z, const char *str);
      // output the given string at the given location
      // NoOp for PostScript device.
#endif

#if (BL_SPACEDIM == 2)
      GraphTool& putString(const IntVect &v, const char *str);
      // output the given string at the given location
      // NoOp for PostScript device.
#endif

      GraphTool& setfgColor(const char *colorname);
      // change the XWindow foreground color to the named color.
      // NoOp for PostScript device.

      GraphTool& setbgColor(const char *colorname);
      // change the XWindow background color to the named color.
      // NOTE: will also blank the screen!
      // NoOp for PostScript device.

      GraphTool& defineCmap( unsigned short *red, unsigned short *green,
                             unsigned short *blue, int num);
      // NoOp for PostScript device.

      GraphTool& setfgColor(int colorTableEntry);
      // change the XWindow foreground color to that specified
      // by the colorTableEntry.
      // NoOp for PostScript device.

      GraphTool& setbgColor(int colorTableEntry);
      // change the XWindow background color to that specified
      // by the colorTableEntry.
      // NOTE: will also blank the screen!
      // NoOp for PostScript device.

      int        getMouse(double &x, double &y) const;
      // query the mouse and return the number of the button
      // pushed and the (x,y) coordinates (in physical coordinates)
      // of where the button was pressed.
      // NoOp for PostScript device.

#ifdef _COMPLEX_IS_OK_
      int        getMouse(complex &z) const;
      // query the mouse and return the number of the button
      // pushed and the complex coordinate (in physical coordinates)
      // of where the button was pressed.
      // NoOp for PostScript device.
#endif

#if (BL_SPACEDIM == 2)
      int        getMouse(IntVect &v) const;
      // query the mouse and return the number of the button
      // pushed and the position (mapped double to IntVect)
      // of where the button was pressed.
      // NoOp for PostScript device.
#endif

      int contour(const double *data, double value,
                  int has_mask, const int *mask,
                  int nx, int ny, int mx, int my,
                  double xlft, double ybot, double xrgt, double ytop);
      // draw contour line of value given data array
      // returns true if success, false if problems.
      // data     = of data to be contoured
      // value    = value to contour
      // has_mask = true if mask array available
      // mask     = array of mask values.  will not contour in masked off cells
      // dims     = specifies dimension of data
      // subrange = specifies subregion of data to be drawn
      // xlft     = position of left   edge of grid in domain
      // ybot     = position of bottom edge of grid in domain
      // xrgt     = position of right  edge of grid in domain
      // ytop     = position of top    edge of grid in domain

#if (BL_SPACEDIM == 2)
      int contour(const double* data, double value,
                  int has_mask, const int *mask,
		  const Box &dims, const Box &subrange,
		  double xlft, double ybot, double xrgt, double ytop);
      // same as above except that:
      // dims     = dimension of arrays (ilo,ihi,jlo,jhi)
      // subrange = subrange of array to plot (irlo,irhi,jrlo,jrhi)

      int contour(const double* data, double value,
                  int has_mask, const int *mask,
		  const Box &dims, const Box &subrange,
		  const Box &position);
      // same as above except that:
      // position = specifies position via current map
#endif

      int contour(const float *data, double value,
                   int has_mask, const int *mask,
                   int nx, int ny, int mx, int my,
                   double xlft, double ybot, double xrgt, double ytop);
      // draw contour line of value given data array. SINGLE PRECISION
      // returns true if success, false if problems.
      // data     = of data to be contoured
      // value    = value to contour
      // has_mask = true if mask array available
      // mask     = array of mask values.  will not contour in masked off cells
      // dims     = specifies dimension of data
      // subrange = specifies subregion of data to be drawn
      // xlft     = position of left   edge of grid in domain
      // ybot     = position of bottom edge of grid in domain
      // xrgt     = position of right  edge of grid in domain
      // ytop     = position of top    edge of grid in domain

#if (BL_SPACEDIM == 2)
      int contour(const float* data, double value,
                   int has_mask, const int *mask,
		   const Box &dims, const Box &subrange,
		   double xlft, double ybot, double xrgt, double ytop);
      // same as above except that:
      // dims     = dimension of arrays (ilo,ihi,jlo,jhi)
      // subrange = subrange of array to plot (irlo,irhi,jrlo,jrhi)

      int contour(const float* data, double value,
                   int has_mask, const int *mask,
		   const Box &dims, const Box &subrange,
		   const Box &position);
      // same as above except that:
      // position = specifies position via current map
#endif

   protected:
      int       width, height;      // width and height of window
      double    xlo,ylo,xhi,yhi;    // physical domain
      double    x0, y0;             // current position
      double    xcL,xcU,ycL,ycU;    // clipping coordinates
      int       was_out;            // was last point out of range?
      GTDevice  dev;                // encoding of active devices
      char      *name;		    // name of window/file
      XWindow   *win;               // pointer to XWindow object
      PSfile    *ps;		    // pointer to postscript object

#if (BL_SPACEDIM == 2)
      IFrame  frame;		    // current int to double mapping
#endif

      double    clipRatio(double,double,double,double);
      void      initGT(int,int,double,double,double,double,
                       const char*,GTDevice);

};

#endif
