//BL_COPYRIGHT_NOTICE

#ifndef _PSfile_H_
#define _PSfile_H_  "@(#)PSfile.H	1.2\t10/5/92"

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Misc.H>
#include "Boolean.H"
 
class PSfile {
   public:
      PSfile(int wid, int high, const char *name);
      ~PSfile();

      void newPage();
      void movePen(double x, double y);
      void drawLine(double x, double y, int lev=0);
      void setLineWidth(int lw);
      void flush();

   private:
      int      npts;
      int      xfactor;
      int      yfactor;
      int      xcur;
      int      ycur;
      ofstream file;
      char     *filename;

}; 

#endif

