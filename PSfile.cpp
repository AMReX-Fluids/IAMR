//BL_COPYRIGHT_NOTICE

#include "PSfile.H"

void PSfile::flush() {file << endl; }

// this is a width of a linewidth=1 line in points
static const double WIDTH = 0.24; 

// this is the X dimension of units of WIDTH points
static const int    FACTOR = 3196;

// these margins place the entire image within the printable
// region for the laser writer.
static const int    XMARGIN = 52;
static const int    YMARGIN = 76;

#define DRAWPATH if (npts>0) {file << "S\n";}; npts = 0;

PSfile::PSfile(int wid, int high, const char *name)
{
   // create postscript file
   filename = new char[strlen(name) + 5];
   strcpy(filename,name);
   strcat(filename,".ps");
   file.open(filename,ios::out);
   if (!file) {
      cerr << "cant open postscript file: " << filename << endl;
      abort();
   };

   // init scaling factors
   double dmax = double( Max(wid,high) );
   double xlen = double(wid)/dmax;
   double ylen = double(high)/dmax;
   xfactor = int( xlen * FACTOR );
   yfactor = int( ylen * FACTOR );
   xcur = ycur = 0;
   npts = 0;

   // write out header
   file << "%!PS\n";
   file << "%  " << wid << "   " << high << '\n';
   file << WIDTH << ' ' << WIDTH << " scale\n";
   file << "-90 rotate\n";
   file << -FACTOR-XMARGIN << ' ' << YMARGIN << " translate\n";
   file << "1 setlinecap\n";
   file << "1 setlinejoin\n";
   file << "/M {moveto} def\n";
   file << "/L {lineto} def\n";
   file << "/S {stroke} def\n";
   file << endl;
}

PSfile::~PSfile()
{
   DRAWPATH;
   file << "showpage\n";
   file << endl;
   file.close();
   delete filename;
}

void PSfile::newPage()
{
   DRAWPATH;
   file << "copypage\n";
   file << "erasepage\n";
}

void PSfile::movePen(double x, double y)
{
   DRAWPATH;
   xcur = int( x*xfactor );
   ycur = int( y*yfactor );
   npts++;
}

void PSfile::drawLine(double x, double y, int lev)
{
   if (npts == 1) {
      file << xcur << ' ' << ycur << " M\n";
   };
   xcur = int( x*xfactor );
   ycur = int( y*yfactor );
   if (npts > 0) {
      file << xcur << ' ' << ycur << " L\n";
      if(lev == 0) {
        file << "[]" << " 0 setdash " << "\n";
      } else {
//        file << "[" << 9*(4-lev) << "]" << " 0 setdash " << "\n";
        if(lev == 1) {
           file << "[" << 2*(3-lev) << " " << 5*(3-lev) << "]"
                << " 0 setdash " << "\n";
        } else {
          if(lev == 2) {
              file << "[" << 1*(3-lev) << " " << 5*(3-lev) << "]"
                << " 0 setdash " << "\n";
          } else {
            cout << " drawLine only handles 2 levels of refinement... "
                 << " You lose... " << "\n";
            exit(0);
          }
	}
      }
   } else {
      file << xcur << ' ' << ycur << " M\n";
   };
   npts++;
}

void PSfile::setLineWidth(int lw)
{
   if (npts > 0) {
      DRAWPATH;
      npts = 1;
   };
   file << lw << " setlinewidth\n";
}

#undef DRAWPATH

