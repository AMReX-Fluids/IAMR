
/*************************************************************************
  This is David Stevens' netcdf debug stuff
  ************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <strstream.h>
#include <Misc.H>
#include <CoordSys.H>
#include <Geometry.H>
#include <BoxDomain.H>
#include <ParmParse.H>
#include <netcdf.h>
#include <netcdfIO.H>
#include <sys/time.h>

#define LINE       200
#define MAXVDIMS   10
#define BADLIM     1.e20
#define MAXFIELDS  2000
#define MAXDSIZES  8000

#ifdef BL_USE_FLOAT
//#define NC_REAL_REP NC_FLOAT
#define NC_REAL_REP NC_DOUBLE
#else
#define NC_REAL_REP NC_DOUBLE
#endif



// a format function which lets aString use standard c formatting
aString aStringFormat( const char *fmt , ... )
{
    va_list ap;
    static char buf[512];

    va_start( ap , fmt );
    vsprintf( buf , fmt , ap );
    int len  = strlen(buf);
    assert( len < 512 );

    return aString(buf);
}


// ==================================== IO OBJECTS FOLLOW ================

// ====================== netCDF field IO ======================
//
//  DimCDF is a server class that saves a dimension to a file
//
//  FieldIO is a class that creates a netCDF representation of the data
//  in a multifab.
//

DimCDF::DimCDF()
{
    cdfid = 0;
    dimid = 0;
    varid = 0;
    ptr   = NULL;
    delta = 0;
    NX    = 0;
    imin  = 0;
    imax  = 0;
    cell  = 0;
}

DimCDF::DimCDF( int cdfid, const char *name,
                int imin, int imax, int node, REAL delta ) :
        cdfid(cdfid), title( name ),
        imin(imin), imax(imax), cell(cell), delta(delta)

{
    // create an array for the DimCDF
    assert( imin <= imax );
    int i,check;
    REAL vr[2];
    char *valid_range = "valid_range";
    NX  = imax-imin+1;
    ptr = new REAL[NX];

    for ( i = imin ; i <= imax ; i++ ) {
	ptr[i-imin] = delta*(i+0.5*(REAL)(1-node));
    }
    vr[0] = ptr[0];
    vr[1] = ptr[NX-1];

    // make the netcdf calls that define the DimCDF
    dimid = ncdimdef(cdfid,name,NX);
    varid = ncvardef(cdfid,name,NC_REAL_REP,1,&dimid);
    check   = ncattput(cdfid,varid,valid_range,NC_REAL_REP,2,vr);
    assert( dimid != -1 );
    assert( varid != -1 );
    assert( check   != -1 );
}

void DimCDF::Write()
{
    int check;
    long start,count;
    start = 0;
    count = NX;
    check   = ncvarput(cdfid,varid,&start,&count,ptr);
    assert( check != -1 );
}

DimCDF::~DimCDF()
{
    delete ptr;
}


// the empty constructor initializes to bad values
netcdfIO::netcdfIO()
{
    level  = -1;
    ngrids = -1;
    ncomps = -1;
    nioper = -1;
    cdfid  = -1;
}


// fancy constructor which probides more info
netcdfIO::netcdfIO( const char *title,
                    REAL cur_time,
                    const REAL *hx,
                    Array<aString> &names,
                    Array<aString> &units,
                    Array<int>     &ngrow,
                    int num_comp,
                    const BoxArray &grids ) :
        title( title ), nioper(1), level(0),
        ncomps(num_comp), ngrids(grids.length())
{
    int check,i,j,k,ind,n;
    
    // create the netcdf file  
    aString filename = aStringFormat( "%s.cdf", title );
    cdfid = nccreate(filename.c_str(),NC_CLOBBER);
    assert( cdfid != -1 );
    
    // create the record dimension
    aString recname = aStringFormat( "time");
    long recsize = NC_UNLIMITED;
    int recdimid = ncdimdef(cdfid,recname.c_str(),recsize);
    int recvarid = ncvardef(cdfid,recname.c_str(),NC_REAL_REP,1,&recdimid);
    assert( recdimid != -1 && recvarid != -1 ); 
    
    // create the dimensions for each grid and grow factor
    Box b;
    IndexType id;
    aString vartitle;
    int max_ngrow = 3;
    DimCDF **xdim = new DimCDF *[ngrids*max_ngrow];
    DimCDF **ydim = new DimCDF *[ngrids*max_ngrow];
#if (BL_SPACEDIM == 3 )
    DimCDF **zdim = new DimCDF *[ngrids*max_ngrow];
#endif
    ind = 0;
    for ( n = 0 ; n < ngrids ; n++ ) {
        for ( i = 0 ; i < max_ngrow ; i++ ) {
            b  = grow(grids[n],i);
            id = b.ixType();
            const int *lo = b.loVect();
            const int *hi = b.hiVect();
            
            vartitle = aStringFormat( "x_%d_%d", i, n );
            xdim[ind]  = new DimCDF( cdfid, vartitle.c_str(), lo[0], hi[0],
                                   id.test(0), hx[0] );
            
            vartitle = aStringFormat( "y_%d_%d", i, n );
            ydim[ind]  = new DimCDF( cdfid, vartitle.c_str(), lo[1], hi[1],
                                   id.test(1), hx[1] );
#if (BL_SPACEDIM == 3 )
            vartitle = aStringFormat( "z_%d_%d", i, n );
            zdim[ind]  = new DimCDF( cdfid, vartitle.c_str(), lo[2], hi[2],
                                   id.test(2), hx[2] );
#endif
            ind++;
        }
    }
    
    // create the variables for each grid
    char *valid_range = "valid_range";
    char *dim_const   = "dim_const";
    char *boundaries  = "bounds";
    REAL vr[2];   vr[0] = -BADLIM;  vr[1] = BADLIM;
    int vdimid[MAXVDIMS];
    int varid,vdims,vsize;
    nc_type vtype;   vtype = NC_REAL_REP;
    for ( n = 0 ; n < ngrids ; n++ ) {
        for ( i = 0 ; i < ncomps ; i++ ) {

            assert( ngrow[i] < max_ngrow );
            ind   = ngrow[i] + n*max_ngrow;
            vdims = 0;
            vdimid[vdims++] = recdimid;
#if (BL_SPACEDIM == 3)
            vdimid[vdims++] = zdim[ind]->dimid;
#endif
            vdimid[vdims++] = ydim[ind]->dimid;
            vdimid[vdims++] = xdim[ind]->dimid;
            
            // define the variable
            // its name is psi_component_ grid number
            vartitle = aStringFormat( "%s_%d", names[i].c_str(), n );
            assert( vdims == BL_SPACEDIM+1 );
            varid = ncvardef( cdfid, vartitle.c_str(), vtype, vdims, vdimid);
            assert( varid != -1 );
            
            // define the valid range
            WriteRArray( varid, valid_range, vr, 2 );
            
            // define the units
            WriteString( varid, "units", units[i].c_str() );

            // and its bounds
            b  = grow(grids[n],ngrow[n]);
            WriteBox( varid, "bounds", b );
        }
    }
	
    // create the global attributes
    WriteReal(  NC_GLOBAL, "cur_time", cur_time );
    WriteReal(  NC_GLOBAL, "dx",       hx[0] );
    WriteReal(  NC_GLOBAL, "dy",       hx[1] );
#if (BL_SPACEDIM == 3)
    WriteReal(  NC_GLOBAL, "dz",       hx[2] );
#endif
    WriteInt(    NC_GLOBAL, "tdomains",  ngrids );
    WriteInt(    NC_GLOBAL, "ngrids",    ngrids );
    WriteInt(    NC_GLOBAL, "ncomps",    ncomps );
    WriteInt(    NC_GLOBAL, "nio_per",   nioper );  
    WriteInt(    NC_GLOBAL, "ratio",     2      ); // a hack
    WriteInt(    NC_GLOBAL, "MultiFab",  1      );
    WriteInt(    NC_GLOBAL, "spacedim",  BL_SPACEDIM );
    
    // leave define mode
    ncendef(cdfid);

    // write out the time 
    
    // write out the dimension variables
    for ( i = 0 ; i < ngrids*max_ngrow; i++ ) {
	xdim[i]->Write();  delete xdim[i];
    }
    delete xdim;

    for ( j = 0 ; j < ngrids*max_ngrow; j++ ) {
	ydim[j]->Write();  delete ydim[j];
    }
    delete ydim;
#if ( BL_SPACEDIM == 3 )
    for ( k = 0 ; k < ngrids*max_ngrow; k++ ) {
	zdim[k]->Write();  delete zdim[k];
    }
    delete zdim;
#endif

    // sync the data file
    ncsync( cdfid );
}


// barebones constructor which is data only
netcdfIO::netcdfIO( const char *title,
                    const BoxArray &grids, int num_comp, int ngrow ) :
        title( title ), nioper(1), level(0),
        ngrids(grids.length()), ncomps(num_comp)
{
    int check,i,j,k,n;
    
    // create the netcdf file  
    aString filename = aStringFormat( "%s.cdf", title );
    cdfid = nccreate(filename.c_str(),NC_CLOBBER);
    assert( cdfid != -1 );
    
    // create the record dimension
    aString recname = aStringFormat( "time");
    long recsize = NC_UNLIMITED;
    int recdimid = ncdimdef(cdfid,recname.c_str(),recsize);
    int recvarid = ncvardef(cdfid,recname.c_str(),NC_REAL_REP,1,&recdimid);
    assert( recdimid != -1 && recvarid != -1 ); 
    
    // create the dimensions for each grid
    Box b;
    IndexType id;
    aString vartitle;
    DimCDF **xdim = new DimCDF *[ngrids];
    DimCDF **ydim = new DimCDF *[ngrids];
#if (BL_SPACEDIM == 3 )
    DimCDF **zdim = new DimCDF *[ngrids];
#endif
    for ( n = 0 ; n < ngrids ; n++ ) {
        b  = grow(grids[n],ngrow);
        id = b.ixType();
        const int *lo = b.loVect();
        const int *hi = b.hiVect();

        vartitle = aStringFormat( "x_%d", n );
        xdim[n]  = new DimCDF( cdfid, vartitle.c_str(), lo[0], hi[0],
                               id.test(0), 1.0 );
        
        vartitle = aStringFormat( "y_%d", n );
        ydim[n]  = new DimCDF( cdfid, vartitle.c_str(), lo[1], hi[1],
                               id.test(1), 1.0 );
#if (BL_SPACEDIM == 3 )
        vartitle = aStringFormat( "z_%d", n );
        zdim[n]  = new DimCDF( cdfid, vartitle.c_str(), lo[2], hi[2],
                               id.test(2), 1.0 );
#endif
    }
    
    // create the variables for each grid
    char *valid_range = "valid_range";
    char *dim_const   = "dim_const";
    char *boundaries  = "bounds";
    REAL vr[2];   vr[0] = -BADLIM;  vr[1] = BADLIM;
    int vdimid[MAXVDIMS];
    int varid,vdims,vsize;
    nc_type vtype;   vtype = NC_REAL_REP;
    for ( n = 0 ; n < ngrids ; n++ ) {
        b = grow(grids[n],ngrow);
        for ( i = 0 ; i < ncomps ; i++ ) {

            vdims = 0;
            vdimid[vdims++] = recdimid;
#if (BL_SPACEDIM == 3)
            vdimid[vdims++] = zdim[n]->dimid;
#endif
            vdimid[vdims++] = ydim[n]->dimid;
            vdimid[vdims++] = xdim[n]->dimid;
            
            // define the variable
            // its name is psi_component_ grid number
            vartitle = aStringFormat( "psi_%d_%d", i, n );
            assert( vdims == BL_SPACEDIM+1 );
            varid = ncvardef( cdfid, vartitle.c_str(), vtype, vdims, vdimid);
            assert( varid != -1 );
            
            // define the valid range
            check = ncattput( cdfid, varid, valid_range, NC_REAL_REP, 2, vr );
            assert( check != -1 );
            
            // and its bounds
            WriteBox( varid, "bounds", b );
        }
    }
	
    // create the global attributes
    WriteInt(    NC_GLOBAL, "tdomains",  ngrids );
    WriteInt(    NC_GLOBAL, "ngrids",    ngrids );
    WriteInt(    NC_GLOBAL, "ncomps",    ncomps );
    WriteInt(    NC_GLOBAL, "nio_per",   nioper );  
    WriteInt(    NC_GLOBAL, "ratio",     2      ); // a hack
    WriteInt(    NC_GLOBAL, "MultiFab",  1      );
    WriteInt(    NC_GLOBAL, "spacedim",  BL_SPACEDIM );
    
    // leave define mode
    ncendef(cdfid);
  
    // write out the dimension variables
    for ( i = 0 ; i < ngrids; i++ ) {
	xdim[i]->Write();  delete xdim[i];
    }
    delete xdim;

    for ( j = 0 ; j < ngrids; j++ ) {
	ydim[j]->Write();  delete ydim[j];
    }
    delete ydim;
#if ( BL_SPACEDIM == 3 )
    for ( k = 0 ; k < ngrids; k++ ) {
	zdim[k]->Write();  delete zdim[k];
    }
    delete zdim;
#endif

    // sync the data file
    ncsync( cdfid );
}




// barebones constructor which is data only
netcdfIO::netcdfIO( const char *title, const BOX &bounds,
                    int ngrow, int num_comp ) :
        title( title ), nioper(1), level(0), ngrids(1), ncomps(num_comp)
{
    int check,i;
    
    // create the netcdf file  
    aString filename = aStringFormat( "%s.cdf", title );
    cdfid = nccreate(filename.c_str(),NC_CLOBBER);
    assert( cdfid != -1 );
    
    // create the record dimension
    aString recname = aStringFormat( "time");
    long recsize = NC_UNLIMITED;
    int recdimid = ncdimdef(cdfid,recname.c_str(),recsize);
    int recvarid = ncvardef(cdfid,recname.c_str(),NC_REAL_REP,1,&recdimid);
    assert( recdimid != -1 && recvarid != -1 ); 
    
    // create the dimensions for the box
    BOX b         = grow(bounds,ngrow);
    IndexType id  = b.ixType();
    const int *lo = b.loVect();
    const int *hi = b.hiVect();
    DimCDF *xdim = new DimCDF( cdfid, "x_0", lo[0], hi[0],
                               id.test(0), 1.0 );
    
    DimCDF *ydim = new DimCDF( cdfid, "y_0", lo[1], hi[1],
                               id.test(1), 1.0 );
#if (BL_SPACEDIM == 3 )
    DimCDF *zdim = new DimCDF( cdfid, "z_0", lo[2], hi[2],
                               id.test(2), 1.0 );
#endif
    
    // create the variables for each grid
    aString vartitle;
    char *valid_range = "valid_range";
    char *dim_const   = "dim_const";
    char *boundaries  = "bounds";
    REAL vr[2];   vr[0] = -BADLIM;  vr[1] = BADLIM;
    int vdimid[MAXVDIMS];
    int varid,vdims,vsize;
    nc_type vtype;   vtype = NC_REAL_REP;
    for ( i = 0 ; i < ncomps ; i++ ) {

        vdims = 0;
        vdimid[vdims++] = recdimid;
#if (BL_SPACEDIM == 3)
        vdimid[vdims++] = zdim->dimid;
#endif
        vdimid[vdims++] = ydim->dimid;
        vdimid[vdims++] = xdim->dimid;
        
        // define the variable, its name is psi_component
        vartitle = aStringFormat( "psi_%d", i );
        assert( vdims == BL_SPACEDIM+1 );
        varid = ncvardef( cdfid, vartitle.c_str(), vtype, vdims, vdimid);
        assert( varid != -1 );
        
        // define the valid range
        check = ncattput( cdfid, varid, valid_range, NC_REAL_REP, 2, vr );
        assert( check != -1 );
        
        // and its bounds
        WriteBox( varid, "bounds", b );
    }
	
    // create the global attributes
    WriteInt(    NC_GLOBAL, "tdomains",  ngrids );
    WriteInt(    NC_GLOBAL, "ngrids",    ngrids );
    WriteInt(    NC_GLOBAL, "ncomps",    ncomps );
    WriteInt(    NC_GLOBAL, "nio_per",   nioper );  
    WriteInt(    NC_GLOBAL, "ratio",     2      );
    WriteInt(    NC_GLOBAL, "MultiFab",  0      );
    WriteInt(    NC_GLOBAL, "spacedim",  BL_SPACEDIM );
    
    // leave define mode
    ncendef(cdfid);

    // write out the dimension variables and garbage collection
    xdim->Write();  delete xdim;
    ydim->Write();  delete ydim;
#if ( BL_SPACEDIM == 3 )
    zdim->Write();  delete zdim;
#endif

    // sync the data file
    ncsync( cdfid );
}



// the destructor closes the file
netcdfIO::~netcdfIO()
{
    if ( cdfid != -1 )
        ncclose( cdfid );
}



// do IO in netCDF format
void netcdfIO::Write( const MultiFab &data, int it,
                      int dest_comp,
                      int src_comp, int num_comp )
{
    // process the timestepping
    int rec = (it)/nioper;
    if ( it%nioper != 0 ) {
        return;
    }

    // write out the data fields
    REAL time;
    int i,j,n,ind;
    int ndims,nvars,natts,xtendim;
    
    char dimname[LINE];
    long dimsize[MAXDSIZES];
    
    char varname[LINE];
    nc_type vartype;
    int check;
    int vardims;
    int vardimid[4];
    int varatts;
    long start[MAXFIELDS][4];
    long count[MAXFIELDS][4];
    
    // *** figure out what is in the netcdf data file
    ncinquire(cdfid,&ndims,&nvars,&natts,&xtendim);
    if ( nvars > MAXFIELDS ) {
        cout << "Raise MAXFIELDS in netcdfIO.C" << endl;
        exit(0);
    }
    if ( ndims > MAXDSIZES ) {
        cout << "Raise MAXDSIZES in netcdfIO.C" << endl;
        exit(0);
    }
    assert( ngrids >= 1 );
    if ( nvars < ndims+ (ngrids-1)*ncomps + dest_comp + num_comp ) {
        cout << "dest_comp too large in netcdfIO.C" << endl;
        exit(0);
    }
    if ( src_comp+num_comp > data.nComp() ) {
        cout << "src_comp+num_comp too large for FARRAYBOX in netcdfIO.C" << endl;
        exit(0);
    }
    assert( xtendim == 0 );
    for (i = 0; i < ndims; i++) {
        ncdiminq(cdfid,i,dimname,&dimsize[i] );
    }
    
    for (i = 0; i < nvars; i++) {
        ncvarinq(cdfid,i,varname,&vartype,&vardims,vardimid,&varatts);
        assert( vardims <= 4 );
        for (j = 0 ; j < vardims; j++ ) {
            start[i][j] = 0;
            count[i][j] = dimsize[ (vardimid[j]) ];
        }
    }

    // ***** set the record dimension *******
    // check to see that we have the right time record
    rec = ( nioper == 0 ? 0 : it/nioper );
    assert( rec <= dimsize[0] );
    start[0][0] = rec;
    count[0][0] = 1;
    time = it*0.0;
    check  = ncvarput( cdfid, 0, start[0], count[0], &time );
    assert( check != -1 );
    for ( i = ndims ; i < nvars ; i++ ) {
        start[i][0] = rec;
        count[i][0] = 1;
    }
    
    // ******* write out the volume info *********
    for ( n = 0 ; n < ngrids ; n++ ) {
        for ( i = 0 ; i < num_comp ; i++ ) {
            ind = ndims + n*ncomps + dest_comp + i;
            check = ncvarput( cdfid, ind, &start[ind][0], &count[ind][0],
                              data[n].dataPtr(src_comp+i) );
            assert( check != -1 );
        }
    }
  
    // ******** synchronize the data set *********
    ncsync( cdfid );
}




// do IO in netCDF format for a FAB
void netcdfIO::Write( const FARRAYBOX &data, int it,
                      int gridno, int dest_comp,
                      int src_comp, int num_comp )
{
    // process the timestepping
    int rec = (it)/nioper;
    if ( it%nioper != 0 ) {
        return;
    }

    // write out the data fields
    REAL time;
    int i,j,n,ind;
    int ndims,nvars,natts,xtendim;
    
    char dimname[LINE];
    long dimsize[MAXDSIZES];
    
    char varname[LINE];
    nc_type vartype;
    int check;
    int vardims;
    int vardimid[4];
    int varatts;
    long start[MAXFIELDS][4];
    long count[MAXFIELDS][4];
    
    // *** figure out what is in the netcdf dtat file
    ncinquire(cdfid,&ndims,&nvars,&natts,&xtendim);
    if ( nvars > MAXFIELDS ) {
        cout << "Raise MAXFIELDS in netcdfIO.C" << endl;
        exit(0);
    }
    if ( ndims > MAXDSIZES ) {
        cout << "Raise MAXDSIZES in netcdfIO.C" << endl;
        exit(0);
    }
    if ( nvars < ndims+ gridno*ncomps + dest_comp + num_comp ) {
        cout << "dest_comp too large in netcdfIO.C" << endl;
        exit(0);
    }
    if ( src_comp+num_comp > data.nComp() ) {
        cout << "src_comp+num_comp too large for FARRAYBOX in netcdfIO.C" << endl;
        exit(0);
    }
    assert( xtendim == 0 );
    for (i = 0; i < ndims; i++) {
        ncdiminq(cdfid,i,dimname,&dimsize[i] );
    }
    for (i = 0; i < nvars; i++) {
        ncvarinq(cdfid,i,varname,&vartype,&vardims,vardimid,&varatts);
        assert( vardims <= 4 );
        for (j = 0 ; j < vardims; j++ ) {
            start[i][j] = 0;
            count[i][j] = dimsize[ (vardimid[j]) ];
        }
    }

    // ***** set the record dimension *******
    // check to see that we have the right time record
    rec = ( nioper == 0 ? 0 : it/nioper );
    assert( rec <= dimsize[0] );
    start[0][0] = rec;
    count[0][0] = 1;
    time = it*0.0;  // a hack
    check  = ncvarput( cdfid, 0, start[0], count[0], &time );
    assert( check != -1 );
    for ( i = ndims ; i < nvars ; i++ ) {
        start[i][0] = rec;
        count[i][0] = 1;
    }
    
    // ******* write out the FAB data
    for ( i = 0 ; i < num_comp ; i++ ) {
        ind = ndims + gridno*ncomps + dest_comp + i;
        check = ncvarput( cdfid, ind, &start[ind][0], &count[ind][0],
                          data.dataPtr(src_comp+i) );
        assert( check != -1 );
    }
  
    // ******** synchronize the data set *********
    ncsync( cdfid );
}





void netcdfIO::WriteBox( int varid, char *name, BOX& b )
{
    // process the box into a set of integers
    int temp[3*BL_SPACEDIM + 1];
    int ntemp    = 3*BL_SPACEDIM + 1;
    int ind      = 0;
    IndexType ix = b.ixType();
    int i;
    for ( i = 0 ; i < BL_SPACEDIM ; i++ ) {
        temp[ind++] = b.smallEnd(i);
        temp[ind++] = b.bigEnd(i);
        temp[ind++] = (ix.ixType(i) ==  IndexType::CELL ? 1 : 0 );
    }
    temp[ind++] = level;
    assert( ind == ntemp );
    
    // write out the box as a set of integers
    int check;
    check = ncattput( cdfid, varid, name, NC_LONG, ntemp, temp );
    assert(check != -1);
}


void netcdfIO::WriteReal( int varid, char *name, REAL a )
{
    int check;
    check = ncattput( cdfid, varid, name, NC_REAL_REP, 1, &a );
    assert(check != -1);
}


void netcdfIO::WriteInt( int varid, char *name, int a )
{
    int check,ltemp;
    ltemp = a;
    check   = ncattput( cdfid, varid, name, NC_LONG, 1, &ltemp );
    assert(check != -1);
}

void netcdfIO::WriteString( int varid, char *name, const char *a )
{
    int check;
    int len = strlen(a);
    if ( len > 0 ) {
        check = ncattput( cdfid, varid, name, NC_CHAR, len, a );
    } else {
        check = ncattput( cdfid, varid, name, NC_CHAR, 2, "na" );
    }
    assert(check != -1);
}


void netcdfIO::WriteIArray( int varid, char *name, int *a, int len )
{
    // process the array
    int *temp = new int[len];
    for ( int i = 0 ; i < len ; i++ ) {
        temp[i] = (int) a[i];
    }

    // write out the array
    int check;
    check = ncattput( cdfid, varid, name, NC_LONG, len, temp );
    assert( check != -1);
    delete [] temp;
}


void netcdfIO::WriteRArray( int varid, char *name, REAL *a, int len )
{
    // write out the array
    int check;
    check = ncattput( cdfid, varid, name, NC_REAL_REP, len, a );
    assert(check != -1);
}


void netcdfIO::OpenAtts()
{
    // reenter define mode
    int check;
    check = ncredef(cdfid);
    assert( check != -1 );
}

void netcdfIO::CloseAtts()
{
    // leave define mode
    ncendef(cdfid);
}


// ==============================================================
// debug functions for dumping out multifab data
// ==============================================================


// MultiFab debug function
void mfab2cdf( const char *t, MultiFab *d )
{
    netcdfIO *temp = new netcdfIO( t, d->boxArray(), d->nComp(), d->nGrow() );
    temp->Write( *d, 0, 0, 0, d->nComp() );
    delete temp;
}



// FARRAYBOX debug function
void fab2cdf(  const char *t, FARRAYBOX *d )
{
    netcdfIO *temp = new netcdfIO( t, d->box(), 0, d->nComp() );
    temp->Write( *d, 0, 0, 0, 0, d->nComp() );
    delete temp;
}



// FARRAYBOX edge debug function
#if ( BL_SPACEDIM == 2 )
void edge2cdf( const char *base, FARRAYBOX *u, FARRAYBOX *v  )
#else
void edge2cdf( const char *base, FARRAYBOX *u, FARRAYBOX *v, FARRAYBOX *w  )
#endif
{
    aString title;
    netcdfIO *temp;
    
    title = aStringFormat( "%su", base );
    temp  = new netcdfIO( title.c_str(), u->box(), 0, u->nComp() );
    temp->Write( *u, 0, 0, 0, 0, u->nComp() );
    delete temp;

    title = aStringFormat( "%sv", base );
    temp  = new netcdfIO( title.c_str(), v->box(), 0, v->nComp() );
    temp->Write( *v, 0, 0, 0, 0, v->nComp() );
    delete temp;
#if ( BL_SPACEDIM == 3 )
    title = aStringFormat( "%sw", base );
    temp  = new netcdfIO( title.c_str(), w->box(), 0, w->nComp() );
    temp->Write( *w, 0, 0, 0, 0, w->nComp() );
    delete temp;
#endif
}

// Multifab edge debug function
#if ( BL_SPACEDIM == 2 )
void medge2cdf( const char *base, MultiFab *u, MultiFab *v  )
#else
void medge2cdf( const char *base, MultiFab *u, MultiFab *v, MultiFab *w  )
#endif
{
    aString title;
    netcdfIO *temp;
    
    title = aStringFormat( "%su", base );
    temp  = new netcdfIO( title.c_str(), u->boxArray(), u->nComp(), u->nGrow() );
    temp->Write( *u, 0, 0, 0, u->nComp() );
    delete temp;

    title = aStringFormat( "%sv", base );
    temp  = new netcdfIO( title.c_str(), v->boxArray(), v->nComp(), v->nGrow() );
    temp->Write( *v, 0, 0, 0, v->nComp() );
    delete temp;
#if ( BL_SPACEDIM == 3 )
    title = aStringFormat( "%sw", base );
    temp  = new netcdfIO( title.c_str(), w->boxArray(), w->nComp(), w->nGrow() );
    temp->Write( *w, 0, 0, 0, w->nComp() );
    delete temp;
#endif
}



// This function reads a component of a grid for a FARRAYBOX
// from a netCDF file
extern "C" void cdf2fab( const char *prefix, FARRAYBOX *fab, int grd, int comp )
{
    // declare variables
    int i,status, cdfid,varid,ndims,natts,dimids[10];
    nc_type vtype;
    long sz,start[10],count[10];
    aString fname,vname;
    fname = aStringFormat( "%s.cdf",    prefix );
    if (grd < 0 )
        vname = aStringFormat( "psi_%d", comp );
    else
        vname = aStringFormat( "psi_%d_%d", comp, grd );

    // process the data file parameters
    cdfid  = ncopen( fname.c_str(), NC_NOWRITE );  assert(cdfid != -1 );
    varid  = ncvarid( cdfid, vname.c_str());       assert(varid != -1 );
    status = ncvarinq( cdfid, varid, NULL, &vtype, &ndims, dimids, &natts);
    assert( status != -1 );
    
    sz = 1;
    for ( i = 0 ; i < ndims ; i++ ) {
        start[i] = 0;
        status = ncdiminq( cdfid, dimids[i], NULL, &count[i]);
        assert( status != -1 );
        sz *= count[i];
    }

    // process the target variable parameters
    BOX b      = fab->box();          assert( sz   == b.numPts() );
    int ncomps = fab->nComp();        assert( -1 < comp && comp < ncomps );
    REAL *data = fab->dataPtr(comp);  assert( data != NULL );
    
    // read the data and close the file
    status = ncvarget( cdfid, varid, start, count, data );
    assert( status != -1 );
    ncclose( cdfid );
}


// fill a multifab from a data file
extern "C" void cdf2mfab( const char *prefix, MultiFab *target )
{
    for ( int grd = 0 ; grd < target->length() ; grd++ ) {
        for ( int comp = 0 ; comp < target->nComp() ; comp++ ) {
            FARRAYBOX &fab = (*target)[grd];
            cdf2fab( prefix, &fab, grd, comp );
        }
    }
}



