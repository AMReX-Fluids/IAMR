#if (BL_SPACEDIM==2) || (BL_SPACEDIM==3)

#ifndef  _Raster_H_
#define  _Raster_H_ 

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <string.h>
#include <math.h>

#include <REAL.H>
#include <aString.H>
#include <List.H>
#include <Box.H>
#include <FArrayBox.H>

// forward ref
class Amr;

class RastFrame {
public:
    RastFrame(aString &file_nm, aString& var_nm, int out_typ, int freq,
              int lev, BOX& graph_box, Amr& amrsys);
      
    ~RastFrame();

private:
    friend class Raster;
    friend ostream& operator << (ostream &os, Raster &r);
      
    int        interval;    // how often to plot
    int        o_type;      // output type
    int        level;       // level of plot
    BOX        box;         // box defining plot region at given level
    aString    file_name;   // root name of file
    aString    var_name;    // name of variable to plot
};

class Raster {
    friend ostream& operator << (ostream &os, Raster &r);

public:

      Raster(Amr& sys);

      ~Raster();

      // manually clear raster space
      void clear();

      void addVar(char* file_nm, char* var_nm, char* out_typ,
                  int freq, int lev, BOX& graph_box);

      // dump plot files, return number of frames drawn
      int draw(REAL time, int nstep, int force_draw = 0);

private:
    List<RastFrame*>   frames;
    Amr       *amrsys;
    int        time_plot;    // 0 = step interval, 1 = time interval
    int	       plot_cnt;     // plotfile number
    REAL       plot_dt;      // plot time interval

    int        dump_locations;
    void       gridLocations(int);
    static int verbose;

};

#endif

#endif
