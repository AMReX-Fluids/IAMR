#ifndef  _Contour_H_
#define  _Contour_H_ 

#if    (HAS_WINDOWS)

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <string.h>
#include <math.h>

#include <REAL.H>
#include <aString.H>
#include <List.H>
#include <Box.H>
#include <GraphTool.H>

static const int WINSIZE = 600;

class Amr;

class ContFrame {
   public:
      friend class Contour;
      friend ostream& operator << (ostream &os, Contour &cont);

      ContFrame(const aString& file_nm, const aString& var_nm,
		int dev, int freq, int ncont, int lev,
                int draw_grids, 
#if (BL_SPACEDIM == 3)
                double Xslice, double Yslice, double Zslice, 
#endif
                int sz=WINSIZE);

     ~ContFrame();
      
   private:
      aString    file_name;   // root name of variable
      aString    var_name;    // name of quantity to draw
      GTDevice   device;      // which devices to use
      int        interval;    // how often to plot
      int        level;       // level of plot
      int        num_cont;    // number of contour lines to draw
      int	 show_grids;  // show grid locations?
      int        size;	      // number of pixels in window
      GraphTool  *win[3];     // window that we will draw in
#if (BL_SPACEDIM == 3)
      double     xslice,yslice,zslice;  // for 3d slice based on xmin=0
#endif
};

class Contour {
   friend ostream& operator << (ostream &os, Contour &cont);
   public:
      Contour(Amr& amrsys);

      ~Contour();

      // manually clear Contour space
      void clear();

      // draw frames, return number of frames drawn
      int draw(REAL time, int nstep, int force_draw = 0);

      int numXWin() {return num_X_windows; };  

   private:
      List<ContFrame*>  frames;
      Amr             *amrptr;
      int              num_X_windows;

      static int verbose;
      static REAL tol;

#if (BL_SPACEDIM == 2)
      void draw_vector_field(const ContFrame& r, REAL time);

#elif (BL_SPACEDIM == 3)
      void draw_vector_field(const ContFrame& r, REAL time, int slice,
                             double slice_val,
                             double data_min,
                             double data_max);
#endif

};

#endif

#endif
