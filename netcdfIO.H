/******************** INTRODUCTION *******************************
  This header file processes the IO of diagnostic soundings and scalars
  ****************************************************************/

#ifndef NETCDFIO_H
#define NETCDFIO_H


// =================== IO classes ================

class DimCDF { // create a dimension in a netcdf file
private:
    aString title;
    int cdfid,varid;
    REAL *ptr;
    REAL delta;
    int NX;
    int imin,imax,cell;
public:
    int dimid;
    DimCDF();
    DimCDF( int, const char *, int, int, int, REAL );
    void Write();
    ~DimCDF();
};


class netcdfIO { // write a multifab to a netcdf file
private:
    aString title;
    int cdfid;
    int nioper;
    int level;
    int ngrids;
    int ncomps;
public:
    // empty constructor
    netcdfIO();

    // fancy constructor for more information
    netcdfIO( const char *,
              REAL cur_time,
              const REAL *hx,
              Array<aString> &names,
              Array<aString> &units,
              Array<int>     &ngrow,
              int num_comp,
              const BoxArray &grids );

    // constructor for bare bones IO data only
    netcdfIO( const char *title, const BOX &bounds,
              int num_comp, int ngrow );
    netcdfIO( const char *title, const BoxArray &grids,
              int num_comp, int ngrow );

    // access functions
    void WriteBox(    int, char *, BOX& ); 
    void WriteReal(   int, char *, REAL ); 
    void WriteInt(    int, char *, int  ); 
    void WriteString( int, char *, const char * );
    void WriteIArray( int, char *, int *, int ); 
    void WriteRArray( int, char *, REAL *, int ); 
    void OpenAtts();
    void CloseAtts();
    
    // write out a multifab
    void Write( const MultiFab &data, int it,
                int dest_comp,
                int src_comp, int num_comp );

    // write out a farraybox
    void Write( const FARRAYBOX &data, int it,
                int gridno, int dest_comp,
                int src_comp, int num_comp );

    // the destructor
    ~netcdfIO();
};



// ==============================================================
// debug functions for dumping out boxlib data
// ==============================================================

extern "C" void mfab2cdf( const char *prefix, MultiFab *mfab  );
extern "C" void fab2cdf(  const char *prefix, FARRAYBOX *fab );

extern "C" void cdf2mfab( const char *prefix, MultiFab *mfab  );
extern "C" void cdf2fab(  const char *prefix, FARRAYBOX *fab, int grd, int comp );

#if ( BL_SPACEDIM == 2 )
extern "C" void edge2cdf(  const char *, FARRAYBOX *, FARRAYBOX * );
extern "C" void medge2cdf( const char *, MultiFab *,  MultiFab * );
#else
extern "C" void edge2cdf(  const char *, FARRAYBOX *, FARRAYBOX *, FARRAYBOX * );
extern "C" void medge2cdf( const char *, MultiFab *,  MultiFab *,  MultiFab *  );
#endif


#endif





