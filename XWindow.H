//BL_COPYRIGHT_NOTICE

#ifndef _XWindow_H_
#define _XWindow_H_ "@(#)XWindow.H	1.2\t10/5/92"

#include "XWinExterns.H"

const long maxIntensity = 65535;

class XWindow {
   public:
      XWindow(int,int,const char*,int iscale = 1);
      ~XWindow();

      void newPage();
      // clear screen

      void movePen(int x, int y);
      // move current location

      void drawLine(int x, int y);
      // draw line from current location to this point
      // set current location to this point

      void drawBox(int x1, int y1, int x2, int y2);
      // draw the given rectangle

      int  getMouse(int&, int&);
      // return the button number and position of mouse

      void setLineWidth(int lw);
      // change the linewidth

      void setFont(char*);
      // change the current font

      void putString(int x, int y, const char *str);
      // write a character string at the given location

      void setfgColor(const char*);
      // set the foreground color to the named color

      void setbgColor(const char*);
      // set the background color to the named color

      void setfgColor(int);
      // set the foreground color to this color table entry

      void setbgColor(int);
      // set the background color to this color table entry

      void defineCmap(unsigned short *red, unsigned short *green,
                       unsigned short *blue, int num);
      // add these colors to the current color map

      void flush();
      // flus the XWindow buffer

   private:
      Window          win;
      Display         *display;
      GC              gc;
      XFontStruct     *the_font;
      int             screen;
      int             width, height, scale;
      char            *name;
      int             x0, y0;
      int             nplanes;
      Colormap        cmap;
      unsigned long   *cmap_colors;
      int             num_cmap_colors;
};

#endif
