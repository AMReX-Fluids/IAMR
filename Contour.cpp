#if    (BL_USE_WINDOWS)

#include <limits.h>

#include <Contour.H>
#include <Amr.H>
#include <Geometry.H>
#include <BoxArray.H>
#include <Array.H>
#include <TagBox.H>
#include <ParmParse.H>

#include <MAKESLICE_F.H>

int Contour::verbose = 0;

REAL Contour::tol = 1.e-8;

static char* Xwin_Colors[] = {"black","sea green","red",
                              "gold","orange","blue","cyan"};
static char* Vector_Colors[] = {"dark slate gray","dark green","dark orange",
                                "brown","orange","blue","cyan"};
// to allow more than 6 levels
static int num_xwin_colors = 7;
static int num_vector_colors = 7;

static int border_wid = 2;

ContFrame::ContFrame(const aString& file_nm, const aString& var_nm,
		     int dev, int freq, int ncont, int lev,
		     int draw_grids, 
#if (BL_SPACEDIM == 3)
                     double Xslice, double Yslice, double Zslice, 
#endif
                     int sz)
{
    int i;

    file_name  = file_nm;
    var_name   = var_nm;
    interval   = freq;
    num_cont   = ncont;
    level      = lev;
    show_grids = draw_grids;
    size       = sz;
    for (i=0; i<3; i++) win[i] = NULL;

#if (BL_SPACEDIM == 3)
    xslice     = Xslice;
    yslice     = Yslice;
    zslice     = Zslice;
#endif

    if (dev == 1) {
	device = xWinDevice;
    } else if (dev == 2) {
	device = psDevice;
    } else if (dev == 3) {
	device = xWinDevice | psDevice;
    } else {
	device = noDevice;
    };
}

ContFrame::~ContFrame()
{
int i;
#if (BL_SPACEDIM==2)
    delete win[0];
#else
    for (i=0;i<3;i++)
     delete win[i];
#endif
}

//   ############################################################
//   ##### Contour members
//   ############################################################

Contour::Contour(Amr& amrsys ){

#if (BL_SPACEDIM == 3)
     // we assume slice values start at zero
    double Xslice=-1.0;
    double Yslice=-1.0;
    double Zslice=-1.0;
#endif

    amrptr = &amrsys;
    num_X_windows = 0;

      // parse input file
    ParmParse pp("contour");
    verbose = pp.contains("verbose");
    int lev = -1;
    pp.query("level",lev);
    pp.query("tol",tol);
#if (BL_SPACEDIM==3)
    pp.query("xslice",Xslice);
    pp.query("yslice",Yslice);
    pp.query("zslice",Zslice);
#endif
    int n = pp.countname("plot");
    int k;
    for (k = 0; k < n; k++) {
	int freq = -1;
	pp.getkth("plot",k,freq,3);
	if (freq > 0) {
	    aString fname;
	    aString vname;
	    int dev, nc, sho, sz;
	    pp.getkth("plot",k,fname,0);
	    pp.getkth("plot",k,vname,1);
	    pp.getkth("plot",k,dev,2);
	    pp.getkth("plot",k,nc,4);
	    pp.getkth("plot",k,sho,5);
	    pp.getkth("plot",k,sz,6);
	    ContFrame *cf= new ContFrame(fname,vname,dev,freq,nc,
					 lev,sho,
#if (BL_SPACEDIM == 3)
					 Xslice,Yslice,Zslice,
#endif
					 sz);
	    frames.append(cf);
	    if (dev == 1) num_X_windows++;
	}
    }
}

Contour::~Contour()
{
    clear();
}

void Contour::clear()
{
    if (num_X_windows > 0) {
	cout << "Enter <cr> to remove Contour plots" << endl;
	while (cin.get() != '\n');
    }
    ListIterator<ContFrame*> li(frames);
    for (; li; ++li) {
	ContFrame *c = li();
	delete c;
	frames[li] = 0;
    }
    frames.clear();
    num_X_windows = 0;
}


ostream& operator << (ostream &os, Contour &cont)
{
    ListIterator<ContFrame*> li(cont.frames);
    for (; li; ++li) {
	ContFrame *c = li();
	os << "contplot  " << c->var_name << "  " << c->device
	   << "  " << c->interval
	   << "  " << c->num_cont << "  " << c->level
	   << "  " << c->show_grids << endl;
    };
    return os;
}

#define MAX_LEV 10
int Contour::draw(REAL time, int nstep, int force_draw)
{

    int slice;
    int boxinrange;
    char f_name[200];

#if (BL_SPACEDIM == 3)
    double slice_val,data_min,data_max;
#endif

      // loop through frames, creating and writing images
    int num_frames = 0;
    ListIterator<ContFrame*> li(frames);
    for (; li; ++li) {

#if (BL_SPACEDIM == 2)
    slice = 0;
#elif (BL_SPACEDIM == 3)
    for (slice = 0; slice <= 2; slice++) {
#endif

	ContFrame& r = *(li());

#if (BL_SPACEDIM == 3)
        if ( (slice == 0 && r.zslice != -1.0) ||
             (slice == 1 && r.yslice != -1.0) ||
             (slice == 2 && r.xslice != -1.0) ) {
#endif

	int intval = r.interval;
	char* dummystr = strcpy(f_name, (r.file_name).c_str() );

	if (force_draw || ((intval > 0) && (nstep%intval == 0))) {

		const REAL* dx_lev = amrptr->Geom(0).CellSize();
		double xborder = border_wid*dx_lev[0];
		double yborder = border_wid*dx_lev[1];
		double xmax = Geometry::ProbLength(0)+xborder;
		double ymax = Geometry::ProbLength(1)+yborder;
#if (BL_SPACEDIM==3)
                double zborder = border_wid*dx_lev[2];
                double zmax = Geometry::ProbLength(2)+zborder;
                if (slice==0) {
                 char* tempstr=strcat(f_name,":xy");
                 slice_val=(zmax-zborder)/2.0;
                 if (r.zslice>=0.0) slice_val=r.zslice;
                 data_min=0.0;
                 data_max=zmax-zborder;
                } else if (slice==1) {
                 char* tempstr=strcat(f_name,":xz");
                 slice_val=(ymax-yborder)/2.0;
                 if (r.yslice>=0.0) slice_val=r.yslice;
                 data_min=0.0;
                 data_max=ymax-yborder;
                 ymax=zmax;yborder=zborder;
                } else if (slice==2) {
                 char* tempstr=strcat(f_name,":yz");
                 slice_val=(xmax-xborder)/2.0;
                 if (r.xslice>=0.0) slice_val=r.xslice;
                 data_min=0.0;
                 data_max=xmax-xborder;
                 xmax=ymax;xborder=yborder;
                 ymax=zmax;yborder=zborder;
                }
                if (verbose)
                 cout << "3d slices, slice, slice_val " << slice << " " <<
                  slice_val << endl;
#endif
	      // if window not created yet, create it
	    if (r.win[slice] == NULL) {
		r.win[slice] = new GraphTool(-xborder,-yborder,xmax,ymax,
				             f_name,r.size,r.device);
	    };

	      // increment number of frames drawn
	    num_frames++;

	    r.win[slice]->newPage();

	      // derive data and compute min and max
	    int maxlev = r.level;
	    int amrlev = amrptr->finestLevel();
	    if (maxlev < 0) maxlev = amrlev;
	    maxlev = Min(maxlev,amrlev);
	    if (maxlev > MAX_LEV) {
		cerr << "Contour::draw: maxlev > MAX_LEV" << endl;
		abort();
	    }
	    int n_cont = r.num_cont;
	    double v_min =  1.e30;
	    double v_max = -1.e30;
	    PArray<FARRAYBOX> *soln[MAX_LEV];
	    int i, lev;
	    for (lev = 0; lev <= maxlev; lev++) {
		soln[lev] = amrptr->derive(r.var_name,time,lev);
		const PArray<FARRAYBOX> &grids = *soln[lev];
		int n_boxes = grids.length();
		for (i = 0; i < n_boxes; i++) {
		    REAL g_min = grids[i].min(0);
		    REAL g_max = grids[i].max(0);
		    v_min = Min(v_min,g_min);
		    v_max = Max(v_max,g_max);
		}
	    };

	    REAL v_range = v_max - v_min;
	    REAL v_off = v_min + 0.5*v_range/float(n_cont);
	    if (verbose) {
		cout << "CONTOUR " << r.var_name <<
		    ", min = " << v_min << " max = " << v_max << endl;
	    }
	    if (v_max-v_min < tol) {
		cout << "CONTOUR: constant field  = " << v_min
		     << " for " << r.var_name << endl;
		for (lev = 0; lev <= maxlev; lev++) {
		    delete soln[lev];
		    soln[lev] = NULL;
		}
		continue;
	    }

	      // draw box surrounding problem domain
	    REAL prob_x = Geometry::ProbLength(0);
	    REAL prob_y = Geometry::ProbLength(1);
#if (BL_SPACEDIM==3)
            REAL prob_z = Geometry::ProbLength(2);
            if (slice==1) {
             prob_y=prob_z;
            } else if (slice==2) {
             prob_x=prob_y;
             prob_y=prob_z;
            }
#endif
	    r.win[slice]->setfgColor(Xwin_Colors[0]);
	    r.win[slice]->drawBox(0.0,0.0,prob_x,prob_y);

	    for (lev = 0; lev <= maxlev; lev++) {
		  // set color for this level
	      // to allow more than 6 levels
		r.win[slice]->setfgColor(Xwin_Colors[lev%num_xwin_colors]);

		  // get BoxAssoc to index into array
		const PArray<FARRAYBOX> &grids = *soln[lev];
		const BoxArray &bs = amrptr->boxArray(lev);
		const REAL* hx = amrptr->Geom(lev).CellSize();

		  // draw box borders
		if (r.show_grids > 0) {
		    for (i = 0; i < bs.length(); i++) {
			const BOX&  bx = bs[i];
			int ilo = bx.smallEnd(0);
			int ihi = bx.bigEnd(0);
			int jlo = bx.smallEnd(1);
			int jhi = bx.bigEnd(1);
			double xlft = ilo*hx[0];
			double xrgt = (ihi+1)*hx[0];
			double ybot = jlo*hx[1];
			double ytop = (jhi+1)*hx[1];
                        boxinrange = 1;
#if (BL_SPACEDIM==3)
                        int klo = bx.smallEnd(2);
                        int khi = bx.bigEnd(2);
                        double zbot = klo*hx[2];
                        double ztop = (khi+1)*hx[2];
                        if (slice==0) {
                         if ((zbot>slice_val)||(ztop<=slice_val))
                          boxinrange=0;
                        } else if (slice==1) {
                         if ((ybot>slice_val)||(ytop<=slice_val))
                          boxinrange=0;
                         ybot=zbot;ytop=ztop;
                        } else if (slice==2) {
                         if ((xlft>slice_val)||(xrgt<=slice_val))
                          boxinrange=0;
                         xlft=ybot;
                         xrgt=ytop;
                         ybot=zbot;
                         ytop=ztop;
                        }
#endif
                        if (boxinrange != 0)
			  r.win[slice]->drawBox(xlft,ybot,xrgt,ytop);
		    }
		}

		if (lev == maxlev) {
		    for (i = 0; i < bs.length(); i++) {
			const BOX &bx = bs[i];
			int ilo = bx.smallEnd(0);
			int ihi = bx.bigEnd(0);
			int jlo = bx.smallEnd(1);
			int jhi = bx.bigEnd(1);
			double xlft = (ilo + 0.5)*hx[0];
			double xrgt = (ihi + 0.5)*hx[0];
			double ybot = (jlo + 0.5)*hx[1];
			double ytop = (jhi + 0.5)*hx[1];

                        boxinrange=1;
#if (BL_SPACEDIM==2)
                        const REAL* dat = grids[i].dataPtr();
                        const BOX bxreal(bx);
#elif (BL_SPACEDIM==3)
                        double xlft2= (ilo)*hx[0];
                        double xrgt2= (ihi+1.0)*hx[0];
                        double ybot2= (jlo)*hx[1];
                        double ytop2= (jhi+1.0)*hx[1];

                        int klo = bx.smallEnd(2);
                        int khi = bx.bigEnd(2);

                        double zbot = (klo+0.5)*hx[2];
                        double ztop = (khi+0.5)*hx[2];
                        double zbot2= (klo)*hx[2];
                        double ztop2= (khi+1.0)*hx[2];

                        BOX bxthin(bx);
                        bxthin.setSmall(2,0);
                        bxthin.setBig(2,0);
                        if (slice==0) {
                         if ((zbot2>slice_val)|| (ztop2<=slice_val))
                          boxinrange=0;
                        } else if (slice==1) {
                         if ((ybot2>slice_val)||(ytop2<=slice_val))
                          boxinrange=0;
                         ybot=zbot;ytop=ztop;
                         bxthin.setSmall(1,klo);
                         bxthin.setBig(1,khi);
                        } else if (slice==2) {
                         if ((xlft2>slice_val)||(xrgt2<=slice_val))
                          boxinrange=0;
                         xlft=ybot;
                         xrgt=ytop;
                         ybot=zbot;
                         ytop=ztop;
                         bxthin.setSmall(1,klo);
                         bxthin.setBig(1,khi);
                         bxthin.setSmall(0,jlo);
                         bxthin.setBig(0,jhi);
                        }
                        FARRAYBOX thindata;
                        thindata.resize(bxthin,1);

                        const int* lo1=bx.loVect();
                        const int* hi1=bx.hiVect();
                        const int* lo2=bxthin.loVect();
                        const int* hi2=bxthin.hiVect();

                        if (boxinrange!=0)
                         FORT_MAKESLICE(grids[i].dataPtr(),lo1,hi1,
                                        thindata.dataPtr(),lo2,hi2,
                                        &data_min,&data_max,
                                        &slice_val,&slice,hx);
                        const REAL* dat = thindata.dataPtr();
                        const BOX bxreal(bxthin);
#endif
			const int *mask_array = NULL;
                        if (boxinrange != 0)
			 for (int icont = 0; icont < n_cont; icont++) {
			    double frac = (float) icont / (float) n_cont;
			    double value = v_off + frac*(v_max - v_min);
			    r.win[slice]->contour(dat,value,0,mask_array,
					          bxreal,bxreal,
                                                  xlft,ybot,xrgt,ytop);
			 }
		    }

		} else {

		    IntVect lratio = amrptr->refRatio(lev);

		    for (i = 0; i < bs.length(); i++) {
			const BOX &bx = bs[i];

				// construct mask array.  must be size FAB.
			TAGBOX mask(bx);
			mask.setVal(0);
			const BoxArray &bsf = amrptr->boxArray(lev+1);
                        int j;
			for (j = 0; j < bsf.length(); j++) {
			    BOX b_crse(coarsen(bsf[j],lratio));
			    if (b_crse.intersects(bx)) {
                                b_crse &= bx;
				mask.setVal(1,b_crse,0);
			    }
			}

			int ilo = bx.smallEnd(0);
			int ihi = bx.bigEnd(0);
			int jlo = bx.smallEnd(1);
			int jhi = bx.bigEnd(1);
			double xlft = (ilo + 0.5)*hx[0];
			double xrgt = (ihi + 0.5)*hx[0];
			double ybot = (jlo + 0.5)*hx[1];
			double ytop = (jhi + 0.5)*hx[1];
                        double xlft2= (ilo)*hx[0];
                        double xrgt2= (ihi+1.0)*hx[0];
                        double ybot2= (jlo)*hx[1];
                        double ytop2= (jhi+1.0)*hx[1];
                        boxinrange = 1;

#if (BL_SPACEDIM==2)
                        const int *mask_array = mask.dataPtr();
                        const REAL *dat = grids[i].dataPtr();
                        BOX bxreal(bx);
#elif (BL_SPACEDIM == 3)
                        int klo = bx.smallEnd(2);
                        int khi = bx.bigEnd(2);
                        double zbot = (klo+0.5)*hx[2];
                        double ztop = (khi+0.5)*hx[2];
                        double zbot2= (klo)*hx[2];
                        double ztop2= (khi+1.0)*hx[2];
                        BOX bxthin(bx);
                        bxthin.setSmall(2,0);
                        bxthin.setBig(2,0);
                        if (slice==0) {
                         if ((zbot2>slice_val)||(ztop2<=slice_val))
                          boxinrange=0;
                        } else if (slice==1) {
                         if ((ybot2>slice_val)||(ytop2<=slice_val))
                          boxinrange=0;
                         ybot=zbot;ytop=ztop;
                         bxthin.setSmall(1,klo);
                         bxthin.setBig(1,khi);
                        } else if (slice==2) {
                         if ((xlft2>slice_val)||(xrgt2<=slice_val))
                          boxinrange=0;
                         xlft=ybot;
                         xrgt=ytop;
                         ybot=zbot;
                         ytop=ztop;
                         bxthin.setSmall(1,klo);
                         bxthin.setBig(1,khi);
                         bxthin.setSmall(0,jlo);
                         bxthin.setBig(0,jhi);
                        }
                        FARRAYBOX thindata;
                        TAGBOX thinmask(bxthin);
                        thindata.resize(bxthin,1);
                        const int* lo1=bx.loVect();
                        const int* hi1=bx.hiVect();
                        const int* lo2=bxthin.loVect();
                        const int* hi2=bxthin.hiVect();

                        if (boxinrange!=0) {
                         FORT_MAKESLICE(grids[i].dataPtr(),lo1,hi1,
                                        thindata.dataPtr(),lo2,hi2,
                                        &data_min,&data_max,
                                        &slice_val,&slice,hx);
                         FORT_MAKEMASKSLICE(mask.dataPtr(),lo1,hi1,
                                            thinmask.dataPtr(),lo2,hi2,
                                            &data_min,&data_max,
                                            &slice_val,&slice,hx);
                        }
                        const REAL *dat = thindata.dataPtr();
                        const int *mask_array = thinmask.dataPtr();
                        BOX bxreal(bxthin);
#endif

                        if (boxinrange != 0) {
                         int icont;
			 for (icont = 0; icont < n_cont; icont++) {
			    double frac = (float) icont / (float) n_cont;
			    double value = v_off + frac*(v_max - v_min);
			    r.win[slice]->contour(dat,value,1,mask_array,
					          bxreal,bxreal,
                                                  xlft,ybot,xrgt,ytop);
			 }
			}
		    }  // looping i, all the grids at this level

		}  // if finest level

	    } // loop over levels

	    for (lev = 0; lev <= maxlev; lev++) {
		delete soln[lev];
		soln[lev] = NULL;
	    }

            if (r.show_grids == 2) {
#if (BL_SPACEDIM == 2)
                draw_vector_field(r, time);
#elif (BL_SPACEDIM == 3)
                draw_vector_field(r, time, 
                                  slice, slice_val,
                                  data_min,data_max);
#endif
            }

      }   
#if (BL_SPACEDIM == 3)
      }   // if test
     }   // looping slice
#endif
    }

    return num_frames;

}

#if (BL_SPACEDIM == 2)
void Contour::draw_vector_field(const ContFrame& r, REAL time)
{
    int lev, i, j;
    int maxlev = r.level;
    int amrlev = amrptr->finestLevel();
    if (maxlev < 0) maxlev = amrlev;
    maxlev = Min(maxlev,amrlev);
    if (maxlev > MAX_LEV) {
	cerr << "Contour::draw: maxlev > MAX_LEV" << endl;
	abort();
    }

      // get velocity field
    PArray<FARRAYBOX> *u[MAX_LEV];
    PArray<FARRAYBOX> *v[MAX_LEV];
    for (lev = 0; lev <= maxlev; lev++) {
	u[lev] = amrptr->derive("x_velocity",time,lev);
	v[lev] = amrptr->derive("y_velocity",time,lev);
    };

      // zero out region covered by fine level
    for (lev = 0; lev < maxlev; lev++) {
	const BoxArray &bs = amrptr->boxArray(lev);
	PArray<FARRAYBOX> &U = *u[lev];
	PArray<FARRAYBOX> &V = *v[lev];
	const BoxArray &bsf = amrptr->boxArray(lev+1);
	IntVect lrat = amrptr->refRatio(lev);
	for (i = 0; i < bs.length(); i++) {
	    for (j = 0; j < bsf.length(); j++) {
		BOX bfine(coarsen(bsf[j],lrat));
		bfine &= bs[i];
		if (bfine.ok()) {
		    U[i].setVal(0.0,bfine,0);
		    V[i].setVal(0.0,bfine,0);
		}
	    }
	}
    }

      // compute maximum speed
    REAL smax = 0.0;
    for (lev = 0; lev <= maxlev; lev++) {
	const BoxArray &bs = amrptr->boxArray(lev);
	PArray<FARRAYBOX> &U = *u[lev];
	PArray<FARRAYBOX> &V = *v[lev];
	for (i = 0; i < bs.length(); i++) {
	      // compute max speed
            long t_long = bs[i].numPts();
            assert(t_long < INT_MAX);
	    int npts = int(t_long);
	    const REAL* udat = U[i].dataPtr();
	    const REAL* vdat = V[i].dataPtr();
            int k;
	    for (k = 0; k < npts; k++) {
		REAL s = sqrt(udat[k]*udat[k]+vdat[k]*vdat[k]);
		smax = Max(smax,s);
	    }
	}
    }

    if (smax < 1.0e-8) {
	cout << "CONTOUR: zero velocity field" << endl;
    } else {

	  // set color of vectors
	/*
	const REAL* dx0 = amrptr->Geom(0).CellSize();
	int  percent = 10;
	REAL prob_x = Geometry::ProbLength(0);
	REAL Amax = prob_x/double(percent);
	REAL eps = 1.0e-3;
	double alen = 0.25;
	int  nxprob = (int) (prob_x/dx0[0]);
	int  stride = nxprob/(2*percent);
	if (stride < 1) stride = 1;
	*/
        // --------- new sussman code
        const REAL* dx0 = amrptr->Geom(0).CellSize();
        int maxpoints= 20;  // partition longest side into 20 parts
        double maxlen=Geometry::ProbLength(0);
        for (i = 1; i < BL_SPACEDIM; i++)
         if (Geometry::ProbLength(i)>maxlen)
           maxlen=Geometry::ProbLength(i);

        double sight_h=maxlen/maxpoints;
        int stride[BL_SPACEDIM];
        for (i = 0; i < BL_SPACEDIM; i++) {
          stride[i]= (int) (sight_h/dx0[i]);
          if (stride[i]<1) stride[i]=1;
        }
        REAL Amax=1.25*sight_h;
        REAL eps = 1.0e-3;
        double alen = 0.25;
        // --------- end new sussman code

	  // draw vectors
	for (lev = 0; lev <= maxlev; lev++) {
	  //r.win[0]->setfgColor(Vector_Colors[lev]);
	  // to allow more than 6 levels
	  r.win[0]->setfgColor(Vector_Colors[lev%num_vector_colors]);
	    if (lev > 0) {
	      for(i = 0; i < BL_SPACEDIM; i++) {
	        stride[i] *= 2;
	      }
	    }
	    const REAL* hx = amrptr->Geom(lev).CellSize();
	    const BoxArray &bs = amrptr->boxArray(lev);
	    PArray<FARRAYBOX> &U = *u[lev];
	    PArray<FARRAYBOX> &V = *v[lev];
            int k;
	    for (k = 0; k < bs.length(); k++) {
		const BOX& bx = bs[k];
		int ilo = bx.smallEnd(0);
		int ihi = bx.bigEnd(0);
		int nx = ihi-ilo+1;
		int jlo = bx.smallEnd(1);
		int jhi = bx.bigEnd(1);
		double xlft = (ilo + 0.5)*hx[0];
		double ybot = (jlo + 0.5)*hx[1];
		const REAL* udat = U[k].dataPtr();
		const REAL* vdat = V[k].dataPtr();
		for (j = jlo; j <= jhi; j+=stride[1]) {
		    double y = ybot + (j-jlo)*hx[1];
		    for (i = ilo; i <= ihi; i+=stride[0]) {
			double x = xlft + (i-ilo)*hx[0];
			double u1 = udat[i-ilo + nx*(j-jlo)];
			double v1 = vdat[i-ilo + nx*(j-jlo)];
			double s = sqrt(u1*u1 + v1*v1);
			if (s < eps) continue;
			double a = Amax*(u1/smax);
			double b = Amax*(v1/smax);
			double nlen = sqrt(a*a + b*b);
			double x2 = x + a;
			double y2 = y + b;
			r.win[0]->movePen(x,y);
			r.win[0]->drawLine(x2,y2);
			double p1 = x2 - alen*a;
			double p2 = y2 - alen*b;
			p1 = p1 - (alen/2.0)*b;
			p2 = p2 + (alen/2.0)*a;
			r.win[0]->drawLine(p1,p2);
		    }
		}

	    }
	}

    } // end else

      // free memory
    for (lev = 0; lev <= maxlev; lev++) {
	delete u[lev];
	delete v[lev];
	u[lev] = v[lev] = NULL;
    }

}
#elif (BL_SPACEDIM == 3)
void Contour::draw_vector_field(const ContFrame& r, REAL time,
                                int slice,
                                double slice_val,
                                double data_min, double data_max)
{
    int boxinrange;
    int lev, i, j;
    int maxlev = r.level;
    int amrlev = amrptr->finestLevel();
    if (maxlev < 0) maxlev = amrlev;
    maxlev = Min(maxlev,amrlev);
    if (maxlev > MAX_LEV) {
	cerr << "Contour::draw: maxlev > MAX_LEV" << endl;
	abort();
    }

      // get velocity field
    PArray<FARRAYBOX> *u[MAX_LEV];
    PArray<FARRAYBOX> *v[MAX_LEV];
    PArray<FARRAYBOX> *w[MAX_LEV];
    for (lev = 0; lev <= maxlev; lev++) {
       u[lev] = amrptr->derive("x_velocity",time,lev);
       v[lev] = amrptr->derive("y_velocity",time,lev);
       w[lev] = amrptr->derive("z_velocity",time,lev);
    }

      // zero out region covered by fine level
    for (lev = 0; lev < maxlev; lev++) {
	const BoxArray &bs = amrptr->boxArray(lev);
	PArray<FARRAYBOX> &U = *u[lev];
	PArray<FARRAYBOX> &V = *v[lev];
	PArray<FARRAYBOX> &W = *w[lev];
	const BoxArray &bsf = amrptr->boxArray(lev+1);
	IntVect lrat = amrptr->refRatio(lev);
	for (i = 0; i < bs.length(); i++) {
	    for (j = 0; j < bsf.length(); j++) {
		BOX bfine(coarsen(bsf[j],lrat));
		bfine &= bs[i];
		if (bfine.ok()) {
		    U[i].setVal(0.0,bfine,0);
		    V[i].setVal(0.0,bfine,0);
		    W[i].setVal(0.0,bfine,0);
		}
	    }
	}
    }

      // compute maximum speed
    REAL smax = 0.0;
    for (lev = 0; lev <= maxlev; lev++) {
	const BoxArray &bs = amrptr->boxArray(lev);
	PArray<FARRAYBOX> &U = *u[lev];
	PArray<FARRAYBOX> &V = *v[lev];
	PArray<FARRAYBOX> &W = *w[lev];
	for (i = 0; i < bs.length(); i++) {
	      // compute max speed
            long t_long = bs[i].numPts();
            assert(t_long < INT_MAX);
	    int npts = t_long;
	    const REAL* udat = U[i].dataPtr();
	    const REAL* vdat = V[i].dataPtr();
	    const REAL* wdat = W[i].dataPtr();
	    for (int k = 0; k < npts; k++) {
		REAL s = udat[k]*udat[k]+vdat[k]*vdat[k];
		s+=wdat[k]*wdat[k];
		s=sqrt(s);
		smax = Max(smax,s);
	    }
	}
    }

    if (smax < 1.0e-8) {
	cout << "CONTOUR: zero velocity field" << endl;
    } else {

	  // set color of vectors
        const REAL* dx0 = amrptr->Geom(0).CellSize();
        int maxpoints= 20;  // partition longest side into 20 parts
        double maxlen=Geometry::ProbLength(0);
        for (i=1;i<BL_SPACEDIM;i++)
         if (Geometry::ProbLength(i)>maxlen)
           maxlen=Geometry::ProbLength(i);

        double sight_h=maxlen/maxpoints;
        int stride[BL_SPACEDIM];
        for (i=0;i<BL_SPACEDIM;i++) {
          stride[i]=sight_h/dx0[i];
          if (stride[i]<1) stride[i]=1;
        }
        REAL Amax=1.25*sight_h;
        REAL eps = 1.0e-3;
        double alen = 0.25;

	  // draw vectors
	for (lev = 0; lev <= maxlev; lev++) {
	    r.win[slice]->setfgColor(Vector_Colors[lev]);
	    if (lev > 0) {
	      for(i=0;i<BL_SPACEDIM;i++) {
	        stride[i] *= 2;
	      }
	    }
	    const REAL* hx = amrptr->Geom(lev).CellSize();
	    const BoxArray &bs = amrptr->boxArray(lev);
	    PArray<FARRAYBOX> &U = *u[lev];
	    PArray<FARRAYBOX> &V = *v[lev];
	    PArray<FARRAYBOX> &W = *w[lev];
	    for (int k = 0; k < bs.length(); k++) {
		const BOX& bx = bs[k];
                int ilo = bx.smallEnd(0);
		int ihi = bx.bigEnd(0);
                int jlo = bx.smallEnd(1);
		int jhi = bx.bigEnd(1);
		double xlft = (ilo + 0.5)*hx[0];
                double xrgt = (ihi + 0.5)*hx[0];
		double ybot = (jlo + 0.5)*hx[1];
                double ytop = (jhi + 0.5)*hx[2];
                double xlft2= (ilo)*hx[0];
                double xrgt2= (ihi+1.0)*hx[0];
                double ybot2= (jlo)*hx[1];
                double ytop2= (jhi+1.0)*hx[2];

		int nx = ihi-ilo+1;
		int xindex=0;
		int yindex=1;
                boxinrange=1;
		int klo = bx.smallEnd(2);
		int khi = bx.bigEnd(2);
		double zbot = (klo+0.5)*hx[2];
                double ztop = (khi+0.5)*hx[2];
		double zbot2= (klo)*hx[2];
                double ztop2= (khi+1.0)*hx[2];
		BOX bxthin(bx);
		bxthin.setSmall(2,0);
		bxthin.setBig(2,0);
		FARRAYBOX UHOLD,VHOLD;
                if (slice==0) {
                 if ((zbot2>slice_val)||(ztop2<=slice_val))
                  boxinrange=0;
                } else if (slice==1) {
                 if ((ybot2>slice_val)||(ytop2<=slice_val))
                  boxinrange=0;
		 bxthin.setSmall(1,klo);
		 bxthin.setBig(1,khi);
		 jlo=klo;
		 jhi=khi;
		 ybot=zbot;
                 xindex=0;
		 yindex=2;
                } else if (slice==2) {
                 if ((xlft2>slice_val)||(xrgt2<=slice_val))
                  boxinrange=0;
		 bxthin.setSmall(0,jlo);
		 bxthin.setBig(0,jhi);
		 bxthin.setSmall(1,klo);
		 bxthin.setBig(1,khi);
		 ilo=jlo;
		 ihi=jhi;
		 nx=ihi-ilo+1;
		 jlo=klo;
		 jhi=khi;
		 xlft=ybot;
		 ybot=zbot;
		 xindex=1;
		 yindex=2;
		}
                if (boxinrange!=0) {

		 UHOLD.resize(bxthin,1);
		 VHOLD.resize(bxthin,1);

                 const int* lo1=bx.loVect();
                 const int* hi1=bx.hiVect();
                 const int* lo2=bxthin.loVect();
                 const int* hi2=bxthin.hiVect();

                 if (slice==0) {
                  FORT_MAKESLICE(U[k].dataPtr(),lo1,hi1,UHOLD.dataPtr(),
                   lo2,hi2,&data_min,&data_max,&slice_val,&slice,hx);
                  FORT_MAKESLICE(V[k].dataPtr(),lo1,hi1,VHOLD.dataPtr(),
                   lo2,hi2,&data_min,&data_max,&slice_val,&slice,hx);
                 } else if (slice==1) {
 		  FORT_MAKESLICE(U[k].dataPtr(),lo1,hi1,UHOLD.dataPtr(),
                   lo2,hi2,&data_min,&data_max,&slice_val,&slice,hx);
                  FORT_MAKESLICE(W[k].dataPtr(),lo1,hi1,VHOLD.dataPtr(),
                   lo2,hi2,&data_min,&data_max,&slice_val,&slice,hx);
                 } else if (slice==2) {
		  FORT_MAKESLICE(V[k].dataPtr(),lo1,hi1,UHOLD.dataPtr(),
                   lo2,hi2,&data_min,&data_max,&slice_val,&slice,hx);
                  FORT_MAKESLICE(W[k].dataPtr(),lo1,hi1,VHOLD.dataPtr(),
                   lo2,hi2,&data_min,&data_max,&slice_val,&slice,hx);
		 }
                 const REAL* udat = UHOLD.dataPtr();
                 const REAL* vdat = VHOLD.dataPtr();
		 for (j = jlo; j <= jhi; j+=stride[yindex]) {
		    double y = ybot + (j-jlo)*hx[yindex];
		    for (i = ilo; i <= ihi; i+=stride[xindex]) {
			double x = xlft + (i-ilo)*hx[xindex];
			double u1 = udat[i-ilo + nx*(j-jlo)];
			double v1 = vdat[i-ilo + nx*(j-jlo)];
			double s = sqrt(u1*u1 + v1*v1);
			if (s < eps) continue;
			double a = Amax*(u1/smax);
			double b = Amax*(v1/smax);
			double nlen = sqrt(a*a + b*b);
			double x2 = x + a;
			double y2 = y + b;
			r.win[slice]->movePen(x,y);
			r.win[slice]->drawLine(x2,y2);
			double p1 = x2 - alen*a;
			double p2 = y2 - alen*b;
			p1 = p1 - (alen/2.0)*b;
			p2 = p2 + (alen/2.0)*a;
	                r.win[slice]->drawLine(p1,p2);
		    }  // looping i
		 }  // looping j
                }  // if boxinrange

	    }  // looping k -- the grids at one level
	}  // looping levels

    } // end else (ok to plot vectors if they are big enough)

      // free memory
    for (lev = 0; lev <= maxlev; lev++) {
	delete u[lev];
	delete v[lev];
	u[lev] = v[lev] = NULL;
        delete w[lev];
        w[lev]=NULL;
    }

}
#endif
#endif

