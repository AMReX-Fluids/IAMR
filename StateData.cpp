
#include <StateData.H>
#include <ParallelDescriptor.H>

#if defined(BL_USE_FLOAT) && !defined(BL_ARCH_CRAY)
const REAL INVALID_TIME = -1.0e30;
#else
const REAL INVALID_TIME = -1.0e100;
#endif

const int MFNEWDATA = 0;
const int MFOLDDATA = 1;
#ifdef __GNUG__
    StateData::StateData(const StateData&)
{
    BoxLib::Error("Don't call this function!");
}
    StateData& StateData::operator= (const StateData&)
{
    BoxLib::Error("Don't call this function!");
}
#endif

// -------------------------------------------------------------
StateData::StateData() 
    : bc()
{
   desc = 0;
   new_data = old_data = 0;
   new_time.start = INVALID_TIME;
   new_time.stop  = INVALID_TIME;
   old_time.start = INVALID_TIME;
   old_time.stop  = INVALID_TIME;
}

// -------------------------------------------------------------
StateData::StateData(const BOX& p_domain, const BoxArray& grds,
                     const StateDescriptor* d, REAL time, REAL dt)
{
    define(p_domain, grds, *d, time, dt);
}

// -------------------------------------------------------------
void
StateData::define(const BOX& p_domain, const BoxArray& grds,
		  const StateDescriptor &d, REAL time, REAL dt)
{
    domain = p_domain;
    desc = &d;
    grids.define(grds);
      // convert to proper type
    IndexType typ(desc->getType());
    TimeCenter t_typ(desc->timeType());
    if ( ! typ.cellCentered()) {
	domain.convert(typ);
	grids.convert(typ);
    }
    if (t_typ == Point) {
	new_time.start = new_time.stop = time;
	old_time.start = old_time.stop = time - dt;
    } else {
	new_time.start = time;
	new_time.stop  = time+dt;
	old_time.start = time-dt;
	old_time.stop  = time;
    }
    int ncomp = desc->nComp();
    new_data = new MultiFab(grids,ncomp,desc->nExtra(),Fab_allocate);
    old_data = 0;
    buildBC();
}

// -------------------------------------------------------------
void
StateData::restart(istream &is, const StateDescriptor &d)
{

    desc = &d;
   
    is >> domain;
    grids.define(is);

    is >> old_time.start;
    is >> old_time.stop;
    is >> new_time.start;
    is >> new_time.stop;

    int nsets;
    is >> nsets;

    new_data = new MultiFab(is);
    if (nsets == 2) {
	old_data = new MultiFab(is);
    } else {
	old_data = 0;
    }
    buildBC();
}

// -------------------------------------------------------------
void
StateData::buildBC()
{
    int ncomp = desc->nComp();
    bc.resize(ncomp);
    int i;
    for (i = 0; i < ncomp; i++) {
	bc[i].resize(grids.length());
        int j;
	for (j = 0; j < grids.length(); j++) {
	    BCRec bcr;
	    setBC(grids[j],domain,desc->getBC(i),bcr);
	    bc[i].set(j,bcr);
	}
    }
}

// -------------------------------------------------------------
StateData::~StateData()
{
   desc = 0;
   if (new_data != NULL)
     delete new_data;
   if (old_data != NULL)
     delete old_data;
}

// -------------------------------------------------------------
void
StateData::setTimeLevel(REAL time, REAL dt_old, REAL dt_new)
{
    if (desc->timeType() == Point) {
	new_time.start = new_time.stop = time;
	old_time.start = old_time.stop = time - dt_old;
    } else {
	new_time.start = time;
	new_time.stop  = time+dt_new;
	old_time.start = time-dt_old;
	old_time.stop  = time;
    }
}

// -------------------------------------------------------------
void
StateData::swapTimeLevels(REAL dt)
{
    old_time = new_time;
    if (desc->timeType() == Point) {
	new_time.start += dt;
	new_time.stop  += dt;
    } else {
	new_time.start = new_time.stop;
	new_time.stop += dt;
    }

    Swap(old_data, new_data);
}

// -------------------------------------------------------------
void
StateData::allocOldData()
{
    if (old_data == 0) {
	old_data = new MultiFab(grids,desc->nComp(),desc->nExtra(),
				Fab_allocate);
    }
}

// -------------------------------------------------------------
void
StateData::removeOldData()
{
    delete old_data;
    old_data = 0;
}

// -------------------------------------------------------------
void
StateData::reset()
{
    new_time = old_time;
    old_time.start = old_time.stop = INVALID_TIME;

    Swap(old_data, new_data);
}

// -------------------------------------------------------------
void
StateData::FillBoundary(const REAL* dx, const REALBOX& prob_domain,
		    int src_comp, int num_comp, int do_new)
{
    REAL cur_time;
    if (desc->timeType() == Point) {
	cur_time = new_time.start;
	if (!do_new) cur_time = old_time.start;
    } else {
	cur_time = 0.5*(new_time.start + new_time.stop);
	if (!do_new) cur_time = 0.5*(old_time.start + old_time.stop);
    }
   
    //int g;
    //for (g = 0; g < grids.length(); g++) {
    for(MultiFabIterator new_datamfi(*new_data); new_datamfi.isValid();
	++new_datamfi)
    {
	DependentMultiFabIterator old_datamfi(new_datamfi, *old_data);
	//FARRAYBOX* dest = &((*new_data)[g]);
	//if (!do_new) dest = &((*old_data)[g]);
	FARRAYBOX* dest = &(new_datamfi());
	if (!do_new) dest = &(old_datamfi());
	const BOX& bx = dest->box();
	if (!domain.contains(bx)) {
	    const int* dlo = bx.loVect();
	    const int* dhi = bx.hiVect();
	    const int* plo = domain.loVect();
	    const int* phi = domain.hiVect();
	    REAL xlo[BL_SPACEDIM];
	    BCRec bcr;
	    const REAL* problo = prob_domain.lo();
            int i;
	    for (i = 0; i < BL_SPACEDIM; i++) {
		xlo[i] = problo[i] + dx[i]*(dlo[i]-plo[i]);
	    }
	    for (i = 0; i < num_comp; i++) {
		int sc = src_comp+i;
		REAL* dat = dest->dataPtr(sc);
		setBC(bx,domain,desc->getBC(sc),bcr);
                desc->bndryFill(sc)(dat,ARLIM(dlo),ARLIM(dhi),
                                    plo,phi,dx,xlo,
				    &cur_time,bcr.vect());
	    }	                          
	}
    }

}

// -------------------------------------------------------------
void
StateData::FillBoundary(FARRAYBOX& dest, REAL time, const REAL* dx,
                        const REALBOX& prob_domain,
                        int dest_comp, int src_comp, int num_comp)
{
    BOX dbox(dest.box());
    assert( dbox.ixType() == desc->getType() );
   
    if (domain.contains(dbox)) return;

    dbox &= domain;
    const BOX& bx = dest.box();
    const int* dlo = dest.loVect();
    const int* dhi = dest.hiVect();
    const int* plo = domain.loVect();
    const int* phi = domain.hiVect();
    REAL xlo[BL_SPACEDIM];
    BCRec bcr;
    const REAL* problo = prob_domain.lo();
    int i;
    for (i = 0; i < BL_SPACEDIM; i++) {
	xlo[i] = problo[i] + dx[i]*(dlo[i]-plo[i]);
    }
    for (i = 0; i < num_comp; i++) {
	int dc = dest_comp+i;
	int sc = src_comp+i;
	REAL* dat = dest.dataPtr(dc);
	setBC(bx,domain,desc->getBC(sc),bcr);
	desc->bndryFill(sc)(dat,ARLIM(dlo),ARLIM(dhi),
                            plo,phi,dx,xlo,&time,bcr.vect());
    }	                          
}

// -------------------------------------------------------------
void
StateData::linInterp(FARRAYBOX& dest, const BOX& subbox, REAL time,
		     int src_comp, int dest_comp, int num_comp,
		     bool extrap)
{
    REAL teps = (new_time.start - old_time.start)/1000.0;
    if (desc->timeType() == Point) {
	if (old_data == 0) {
	    teps = new_time.start/10000.0;
	    REAL dt = time - new_time.start;
	    if (dt < 0.0) dt = -dt;
	    if (extrap || dt < teps);
	    new_data->copy(dest,subbox,src_comp,dest_comp,num_comp);
	    return;
	}
                       
      ::linInterp(dest,subbox,*old_data,*new_data,
                  old_time.start,new_time.start,time,
	          src_comp,dest_comp,num_comp,
                  extrap);
   } else {
      if (time > new_time.start-teps && time < new_time.stop+teps) {
         new_data->copy(dest,subbox,src_comp,dest_comp,num_comp);
      } else if (old_data != 0 && time > old_time.start-teps &&
                 time < old_time.stop+teps) {
         old_data->copy(dest,subbox,src_comp,dest_comp,num_comp);
      } else {
          BoxLib::Error("StateData::linInterp: cannot interp");
      }
   }
}



// -------------------------------------------------------------
void StateData::RegisterData(MultiFabCopyDescriptor &multiFabCopyDesc,
                             Array<MultiFabId> &mfid)
{
  mfid.resize(2);
  mfid[MFNEWDATA] = multiFabCopyDesc.RegisterFabArray(new_data);
  mfid[MFOLDDATA] = multiFabCopyDesc.RegisterFabArray(old_data);
}


// -------------------------------------------------------------
void
StateData::linInterpAddBox(MultiFabCopyDescriptor &multiFabCopyDesc,
                     Array<MultiFabId> &mfid,
                     BoxList &unfillableBoxes,
                     Array<FillBoxId>  &returnedFillBoxIds,
                     const BOX &subbox,
                     REAL time, int src_comp, int dest_comp, int num_comp,
                     bool extrap)
{
    REAL teps = (new_time.start - old_time.start)/1000.0;
    if (desc->timeType() == Point) {
        if (old_data == 0) {
            teps = new_time.start/10000.0;
            REAL dt = time - new_time.start;
            if (dt < 0.0) dt = -dt;
            returnedFillBoxIds.resize(1);
            returnedFillBoxIds[0] = multiFabCopyDesc.AddBox(mfid[MFNEWDATA],subbox,
                                        unfillableBoxes,
                                        src_comp, dest_comp, num_comp);
            return;
        }

      ::linInterpAddBox(multiFabCopyDesc, unfillableBoxes, returnedFillBoxIds,
                  subbox, mfid[MFOLDDATA], mfid[MFNEWDATA],
                  old_time.start, new_time.start, time,
                  src_comp, dest_comp, num_comp,
                  extrap);
   } else {
      if (time > new_time.start-teps && time < new_time.stop+teps) {
         returnedFillBoxIds.resize(1);
         returnedFillBoxIds[0] = multiFabCopyDesc.AddBox(mfid[MFNEWDATA], subbox,
                                      unfillableBoxes,
                                      src_comp, dest_comp, num_comp);
      } else if (old_data != 0 && time > old_time.start-teps &&
                 time < old_time.stop+teps)
      {
         returnedFillBoxIds.resize(1);
         returnedFillBoxIds[0] = multiFabCopyDesc.AddBox(mfid[MFOLDDATA], subbox,
                                     unfillableBoxes,
                                     src_comp, dest_comp, num_comp);
      } else {
         BoxLib::Error("StateData::linInterp: cannot interp");
      }
   }
}


// -------------------------------------------------------------
void
StateData::linInterpFillFab(MultiFabCopyDescriptor &multiFabCopyDesc,
                     const Array<MultiFabId> &mfid,
                     const Array<FillBoxId> &fillBoxIds,
                     FARRAYBOX &dest,
                     REAL time, int src_comp, int dest_comp, int num_comp,
                     bool extrap)
{
    REAL teps = (new_time.start - old_time.start)/1000.0;
    if (desc->timeType() == Point) {
        if (old_data == 0) {
            teps = new_time.start/10000.0;
            REAL dt = time - new_time.start;
            if (dt < 0.0) dt = -dt;
            multiFabCopyDesc.FillFab(mfid[MFNEWDATA], fillBoxIds[0], dest);

            return;
        }

      ::linInterpFillFab(multiFabCopyDesc, fillBoxIds,
                         mfid[MFOLDDATA], mfid[MFNEWDATA],
                         dest,
                         old_time.start, new_time.start, time,
                         src_comp, dest_comp, num_comp,
                         extrap);
   } else {
      if (time > new_time.start-teps && time < new_time.stop+teps) {
         multiFabCopyDesc.FillFab(mfid[MFNEWDATA], fillBoxIds[0], dest);
      } else if (old_data != 0 && time > old_time.start-teps &&
                 time < old_time.stop+teps) {
         multiFabCopyDesc.FillFab(mfid[MFOLDDATA], fillBoxIds[0], dest);
      } else {
         BoxLib::Error("StateData::linInterp: cannot interp");
      }
   }
}


// -------------------------------------------------------------
void
StateData::checkPoint(ostream &os, bool dump_old) 
{
    if(dump_old == true && old_data == 0) {
	dump_old = false;
    }
    if(ParallelDescriptor::IOProcessor()) {
      os << domain << '\n';
      grids.writeOn(os);
      os << old_time.start << '\n';
      os << old_time.stop  << '\n';
      os << new_time.start << '\n';
      os << new_time.stop  << '\n';
      int nsets = 1;
      if (dump_old) nsets = 2;
      os << nsets << '\n';
    }
    assert(new_data);
    new_data->writeOn(os);
    if (dump_old) old_data->writeOn(os);
}

void
StateData::printTimeInterval(ostream &os)
{
    os << '[' << old_time.start << ' ' << old_time.stop << "] ["
       << new_time.start << ' ' << new_time.stop << ']' << '\n';
}
