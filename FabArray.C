//BL_COPYRIGHT_NOTICE

//
// $Id: FabArray.C,v 1.2 1997-07-12 00:09:47 vince Exp $
//

#include <Assert.H>
#include <Utility.H>
#include <BoxDomain.H> 


template <class T, class FAB>
FabArray<T, FAB>::FabArray ()
    : n_grow(0),
      n_comp(0),
      fabparray(0,PArrayManage)
{}

template <class T, class FAB>
FabArray<T, FAB>::FabArray (const BoxArray& bxs,
                            int             ngrow,
                            FabAlloc        alloc)
    : fabparray(0,PArrayManage)
{
    define(bxs,ngrow,alloc);
}

template <class T, class FAB>
FabArray<T, FAB>::FabArray (const BoxArray& bxs,
                            int             nvar,
                            int             ngrow,
               FabAlloc alloc)
    : fabparray(0, PArrayManage)
{
    define(bxs,nvar,ngrow,alloc);
}


template <class T, class FAB>
FabArray<T, FAB>::~FabArray ()
{
}

template <class T, class FAB>
bool
FabArray<T, FAB>::ok () const
{
    assert(boxarray.ready());
    bool isok = true;

    //if ( ! fabparray.defined(0))
        //isok = false;
    //else
    {
        //for (i = 0, nvar = fabparray[0].nComp(); i < length() && isok; ++i)
	// need to define nvar for check
	for(ConstFabArrayIterator<T, FAB> fai(*this); fai.isValid() && isok; ++fai)
        {
            if(fabparray.defined(fai.index())) {
                const FAB &f = fai();
                if(f.box() != ::grow(box(fai.index()), n_grow)) {
                    isok = false;
		}
		// cant make this check because nvar cannot be defined as above
                //if (nvar != f.nComp()) {
                    //isok = false;
		//}
            } else {
                isok = false;
	    }
        }
    }

    ParallelDescriptor::ReduceBoolAnd(isok);
    return isok;
}

template <class T, class FAB>
void
FabArray<T, FAB>::define (const BoxArray& bxs,
                          int             ngrow,
                          FabAlloc        alloc)
{
    assert( ! boxarray.ready());
    n_grow = ngrow;
    n_comp = 0;
    boxarray.define(bxs);
    distributionMap.define(ParallelDescriptor::NProcs(), boxarray);
    int nbox = boxarray.length();
    fabparray.resize(nbox);
    if (alloc == Fab_allocate) {
	n_comp = 1;
	AllocFabs(n_comp);
    }
}


template <class T, class FAB>
void
FabArray<T, FAB>::define (const BoxArray& bxs,
                          int             nvar,
                          int             ngrow,
                          FabAlloc        alloc)
{
    assert(!boxarray.ready());
    n_grow = ngrow;
    n_comp = nvar;
    boxarray.define(bxs);
    distributionMap.define(ParallelDescriptor::NProcs(), boxarray);
    int nbox = bxs.length();
    fabparray.resize(nbox);
    if (alloc == Fab_allocate) {
	AllocFabs(nvar);
    }
}


template <class T, class FAB>
void
FabArray<T, FAB>::AllocFabs(int nvar) {
  for(FabArrayIterator<T, FAB> fai(*this); fai.isValid(); ++fai) {
    Box tmp(::grow(fai.validbox(), n_grow));
    FAB *f = new FAB(tmp, nvar);
    if(f == 0) {
      BoxLib::OutOfMemory(__FILE__, __LINE__);
    }
    fabparray.set(fai.index(), f);
  }
}


template <class T, class FAB>
void
FabArray<T, FAB>::setFab (int  boxno,
                          FAB* elem)
{
    //
    // Must check it is of the proper size.
    //
    if (n_comp == 0)
        n_comp = elem->nComp();
    assert(n_comp == elem->nComp());
    assert(boxarray.ready());
    assert(elem->box() == ::grow(boxarray[boxno],n_grow));
    assert(!fabparray.defined(boxno));
    if(distributionMap.ProcessorMap()[boxno] == ParallelDescriptor::MyProc()) {
      fabparray.set(boxno,elem);
    } else {
      ParallelDescriptor::Abort("Error in FabArray<T, FAB>::setFab:  trying to setFab remotely");
    }
}

template <class T, class FAB>
void
FabArray<T, FAB>::setBndry (T val)
{
    setBndry(val, 0, n_comp);
}

template <class T, class FAB>
void
FabArray<T, FAB>::setBndry(T val, int strt_comp, int ncomp)
{
    if (n_grow > 0) {
	for(FabArrayIterator<T, FAB> fai(*this); fai.isValid(); ++fai) {
            fai().setComplement(val, fai.validbox(), strt_comp, ncomp);
	}
    }
}


//
// Self copy for given component only.
//

template <class T, class FAB>
void
FabArray<T, FAB>::copy (const FabArray<T, FAB>& farray) {
  if(farray.boxarray == boxarray) {
    for(FabArrayIterator<T, FAB> fai(*this); fai.isValid(); ++fai) {
      fai().copy(farray[fai.index()]);
    }
  } else {

    // the index into these arrays is the same index of the FabArray
    // these are intersections relative to this FabArray
    Array<BoxArray> boxIntersections(length());
    Array<BoxArray> boxIntersectionsLocal(length());
    Array<BoxArray> boxIntersectionsRemote(length());
    Array<BoxArray> boxNonIntersections(length());
    int myproc = ParallelDescriptor::MyProc();
    FabComTag fabComTag;
    int tagSize = sizeof(FabComTag);
    ParallelDescriptor::SetMessageHeaderSize(tagSize);

    int i, ii, iblocal;
    for(i = 0; i < length(); ++i) {  // loop over fabparray

      boxIntersections[i]    = intersect(farray.boxarray, boxarray[i]);
      boxNonIntersections[i] = complementIn(boxarray[i], farray.boxarray);

      //boxIntersectionsLocal[i].reserve(boxIntersections[i].length());
      if(distributionMap[i] == myproc) {  // fab[i] is local
	iblocal = 0;
	for(ii = 0; ii < farray.boxarray.length(); ++ii) {
	  if(farray.distributionMap[ii] == myproc) {    // both fabs are local
	    Box intersectBox(farray.boxarray[ii] & boxarray[i]);
	    if(intersectBox.ok()) {
              boxIntersectionsLocal[i].resize(iblocal + 1);
              boxIntersectionsLocal[i].set(iblocal, intersectBox);
	      ++iblocal;
	      fabparray[i].copy(farray[ii]);  // do the local copy
	    }
	  }
	}
      } else {                            // fab[i] is remote
	iblocal = 0;
	for(ii = 0; ii < farray.boxarray.length(); ++ii) {
	  if(farray.distributionMap[ii] == myproc) {    // farray[ii] is local
	    Box intersectBox(farray.boxarray[ii] & boxarray[i]);
	    if(intersectBox.ok()) {
              boxIntersectionsRemote[i].resize(iblocal + 1);
              boxIntersectionsRemote[i].set(iblocal, intersectBox);

	      // send the data
	      // send the intersection of farray[ii] to fab[i]s processor

              FAB tempFab(boxIntersectionsRemote[i][iblocal], farray.nComp());
              int tempDestComp = 0;
              int tempSrcComp  = 0;
              Box srcBox(boxIntersectionsRemote[i][iblocal]);
              Box destBox(boxIntersectionsRemote[i][iblocal]);
              tempFab.copy(farray[ii], srcBox, tempSrcComp,
                         destBox, tempDestComp, farray.nComp());

              fabComTag.fromProc = myproc;
              fabComTag.toProc   = distributionMap[i];
              fabComTag.fabIndex = i;
              fabComTag.destComp = tempDestComp;
              fabComTag.nComp    = farray.nComp();
              fabComTag.box      = destBox;

              ParallelDescriptor::SendData(fabComTag.toProc, fabComTag,
					   tempFab.dataPtr(),
                                           fabComTag.box.numPts() *
					   fabComTag.nComp * sizeof(T));

	      ++iblocal;
	    }  // end if(intersectBox.ok())

	  }  // end if(farray[ii] is local)
	}  // end for(ii...)
      }  // end if(fab[i] is local)
    }  // end for(i...)

    ParallelDescriptor::Synchronize();

    // now receive data if any was sent
    int dataWaitingSize;
    while(ParallelDescriptor::GetMessageHeader(dataWaitingSize, &fabComTag))
    {  // data was sent to this processor

      if(myproc != fabComTag.toProc) {
	cerr << "Error:  _in FabArray::copy(FabArray):  myproc!=fabComTag.toProc : "
	     << myproc << " != " << fabComTag.toProc << endl;
	ParallelDescriptor::Abort("Bad fabComTag.toProc");
      }
      int shouldReceiveBytes = fabComTag.box.numPts() * fabComTag.nComp * sizeof(T);
      if(dataWaitingSize != shouldReceiveBytes) {
	cerr << "Error:  _in FabArray::copy(FabArray):  "
	     << "dataWaitingSize != shouldReceiveBytes:  = "
	     << dataWaitingSize << " != " << shouldReceiveBytes << endl;
	ParallelDescriptor::Abort("Bad receive nbytes");
      }
      if( ! fabComTag.box.ok()) {
	cerr << "Error:  _in FabArray::copy(FabArray):  bad fabComTag.box" << endl;
	ParallelDescriptor::Abort("Bad receive nbytes");
      }

      FAB tempFab(fabComTag.box, fabComTag.nComp);
      ParallelDescriptor::ReceiveData(tempFab.dataPtr(),
	       fabComTag.box.numPts() * fabComTag.nComp * sizeof(T));
      int srcComp = 0;
      fabparray[fabComTag.fabIndex].copy(tempFab, fabComTag.box, srcComp,
			 fabComTag.box, fabComTag.destComp, fabComTag.nComp);
    }

    ParallelDescriptor::Synchronize();  // do we need this here?


    if(myproc == 0) {
      for(i = 0; i < length(); i++) {  // loop over fabparray
        cout << "boxIntersections[" << i << "] = "
             << boxIntersections[i] << endl;
      }
      for(i = 0; i < length(); i++) {  // loop over fabparray
        cout << "boxIntersectionsLocal[" << i << "] = "
             << boxIntersectionsLocal[i] << endl;
      }
      for(i = 0; i < length(); i++) {  // loop over fabparray
        cout << "boxIntersectionsRemote[" << i << "] = "
             << boxIntersectionsRemote[i] << endl;
      }
      for(i = 0; i < length(); i++) {  // loop over fabparray
        cout << "boxNonIntersections[" << i << "] = "
             << boxNonIntersections[i] << endl;
      }
    }

  }  // end if(same boxarray...)
}  // end FabArray<T, FAB>::copy (const FabArray<T, FAB>& farray)


template <class T, class FAB>
void
FabArray<T, FAB>::copy (const FabArray<T, FAB>& src,
                        int                     src_comp,
                        int                     dest_comp,
                        int                     num_comp,
                        int                     nghost)
{
    assert(src.boxarray == boxarray);
    assert(nghost <= n_grow);
    assert(nghost <= src.n_grow);
    int myproc = ParallelDescriptor::MyProc();
    for(FabArrayIterator<T, FAB> fai(*this); fai.isValid(); ++fai) {
        Box b(grow(fai.validbox(), nghost));
        fai().copy(src[fai.index()], b, src_comp, b, dest_comp, num_comp);
    }
}




//
// Copies to FABs, note that destination is first arg.
//

template <class T, class FAB>
void
FabArray<T, FAB>::copy (FAB& dest) const
{
    int srccomp  = 0;
    int destcomp = 0;
    copy(dest, dest.box(), srccomp, destcomp, dest.nComp());
}

template <class T, class FAB>
void
FabArray<T, FAB>::copy (FAB&       dest,
                        const Box& subbox) const
{
    int srccomp  = 0;
    int destcomp = 0;
    copy(dest, subbox, srccomp, destcomp, dest.nComp());
}


template <class T, class FAB>
void
FabArray<T, FAB>::copy (FAB& dest,
                        int  src_comp,
                        int  dest_comp,
                        int  num_comp) const
{
    copy(dest, dest.box(), src_comp, dest_comp, num_comp);
}


template <class T, class FAB>
void
FabArray<T, FAB>::copy (FAB&       dest,
                        const Box& subbox,
                        int        src_comp,
                        int        dest_comp,
                        int        num_comp) const
{
    assert(dest_comp + num_comp <= dest.nComp());

    int i, overlapLength = 0, overlapIndex = 0;
    int myproc = ParallelDescriptor::MyProc();
    for(i = 0; i < length(); i++) {
      if(subbox.intersects(boxarray[i])) {
	overlapLength++;
      }
    }

    BoxArray overlapBoxArray;  // make a boxarray of overlaps
    Array<int> overlapProcessors;
    Array<int> overlapFAIndex;
    PArray<FAB> overlapData(PArrayManage);

    if(overlapLength > 0) {
      overlapBoxArray.resize(overlapLength);  // make a boxarray of overlaps
      overlapProcessors.resize(overlapLength);
      overlapFAIndex.resize(overlapLength);
      for(i = 0; i < length(); i++) {
	Box overlapBox = subbox & boxarray[i];
	if(overlapBox.ok()) {
	  overlapBoxArray.set(overlapIndex, overlapBox);
	  overlapProcessors[overlapIndex] = distributionMap[i];
	  overlapFAIndex[overlapIndex] = i;
	  overlapIndex++;
	}
      }

      // cant use a FabArray because they are distributed
      overlapData.resize(overlapBoxArray.length());
      for(int tfi = 0; tfi < overlapData.length(); tfi++) {
	int tfigrow = 0;
	Box tfibox(::grow(overlapBoxArray[tfi], tfigrow));
	FAB *f = new FAB(tfibox, num_comp);
	if(f == 0) {
	  BoxLib::OutOfMemory(__FILE__, __LINE__);
	}
	overlapData.set(tfi, f);
      }

      for(i = 0; i < overlapBoxArray.length(); i++) {  // register for sharing
	ParallelDescriptor::ShareVar(overlapData[i].dataPtr(),
	    overlapBoxArray[i].numPts() * overlapData[i].nComp() * sizeof(T));
      }

      for(i = 0; i < overlapBoxArray.length(); i++) {
        if(myproc == overlapProcessors[i]) {  // local data
	  int overlapDestComp = 0;
	  overlapData[i].copy(fabparray[overlapFAIndex[i]],
			      src_comp, overlapDestComp, num_comp);
	}
      }
    }  // end if(overlapLength > 0)

    ParallelDescriptor::Synchronize();  // for ShareVar

    if(overlapLength > 0) {
      for(i = 0; i < overlapBoxArray.length(); i++) {
	bool overlapDataIsLocal = (overlapProcessors[i] == myproc);
	for(int np = 0; np < ParallelDescriptor::NProcs(); np++) {
          if(overlapDataIsLocal && myproc != np) {  // broadcast data
	    int doffset = 0;
	    ParallelDescriptor::WriteData(np, overlapData[i].dataPtr(),
	      overlapData[i].dataPtr(), doffset,
	      overlapData[i].box().numPts() * overlapData[i].nComp() * sizeof(T));
          }
	}
      }
    }  // end if(overlapLength > 0)

    ParallelDescriptor::Synchronize();  // make sure data has arrived

    if(overlapLength > 0) {
      for(i = 0; i < overlapBoxArray.length(); i++) {
	int overlapSrcComp = 0;
        dest.copy(overlapData[i], overlapBoxArray[i], overlapSrcComp,
				  overlapBoxArray[i], dest_comp, num_comp);
      }

      for(i = overlapBoxArray.length() - 1; i >= 0; i--) {         // unshare in
	ParallelDescriptor::UnshareVar(overlapData[i].dataPtr());  // reverse order
      }
    }  // end if(overlapLength > 0)

    ParallelDescriptor::Synchronize();
}  // end FabArray<T, FAB>::copy(...)




template <class T, class FAB>
void
FabArray<T, FAB>::setVal (T val)
{
    for(FabArrayIterator<T, FAB> fai(*this); fai.isValid(); ++fai) {
        fai().setVal(val);
    }
}

template <class T, class FAB>
void
FabArray<T, FAB>::setVal (T   val,
                          int comp,
                          int num_comp,
                          int nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);
    int k;
    for(FabArrayIterator<T, FAB> fai(*this); fai.isValid(); ++fai) {
        Box b(::grow(fai.validbox(), nghost));
        for(k = 0; k < num_comp; k++) {
            fai().setVal(val, b, comp + k);
	}
    }
}

template <class T, class FAB>
void
FabArray<T, FAB>::setVal (T          val,
                          const Box& region,
                          int        comp,
                          int        num_comp,
                          int        nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);
    int k;
    for(FabArrayIterator<T, FAB> fai(*this); fai.isValid(); ++fai) {
        Box b(::grow(fai.validbox(), nghost));
        b &= region;
        if(b.ok()) {
            for(k = 0; k < num_comp; k++) {
                fai().setVal(val, b, comp + k);
	    }
        }
    }
}



template <class T, class FAB>
void
FabArray<T, FAB>::FillBoundary () {
  FillBoundary(0, n_comp);
}


template <class T, class FAB>
void
FabArray<T, FAB>::FillBoundary (int start_comp, int num_comp) {

/*  original code
    for (int i = 0; i < length(); ++i) {
        FAB& fab = fs[i];
        for (int j=0; j < ba.nborMax(i); j++) {
            int nbor = ba.nborIndex(i,j);
            // For self copy, don't copy into self.
            if (nbor == i)
                continue;
            Box destbox(ba[nbor]);
            destbox &= fab.box();
            fab.copy(fs[nbor],destbox,scomp,destbox,scomp,num_comp);
        }
    }
*/

  FabArrayCopyDescriptor<T, FAB> facd(true);
  FabArrayId faid = facd.RegisterFabArray(this);
  List<FillBoxId> fillBoxIdList;

  for(FabArrayIterator<T, FAB> fai(*this); fai.isValid(); ++fai) {
    for(int j = 0; j < length(); j++) {
      if(j == fai.index()) {  // Dont copy into self.
        continue;
      }
      Box destbox(boxarray[j]);
      destbox &= fai().box();
      if(destbox.ok()) {
        //fai().copy(fs[j],destbox,scomp,destbox,scomp,num_comp);
	BoxList unfilledBoxes;  // unused here
	FillBoxId fbid;
	fbid = facd.AddBox(faid, destbox, unfilledBoxes,
			   start_comp, start_comp, num_comp);
	fillBoxIdList.append(fbid);
      }
    }
  }

  facd.CollectData();

  ListIterator<FillBoxId> fbidli(fillBoxIdList);

  for(FabArrayIterator<T, FAB> fai(*this); fai.isValid(); ++fai) {
    for(int j = 0; j < length(); j++) {
      if(j == fai.index()) {  // Dont copy into self.
        continue;
      }
      Box destbox(boxarray[j]);
      destbox &= fai().box();
      if(destbox.ok()) {
	assert(fbidli);
	FillBoxId fbid = fbidli();
	++fbidli;

	assert(destbox == fbid.box());
	FArrayBox overlapFab(fbid.box(), num_comp);
	facd.FillFab(faid, fbid, overlapFab);

        fai().copy(overlapFab, destbox);
      }
    }
  }
}




// ========================================================================

template <class T, class FAB>
FabArrayIterator<T, FAB>::FabArrayIterator(FabArray<T, FAB> &fabarray,
					   const int nghost)
                         : fabArray(fabarray),
			   nGhost(nghost),
			   currentIndex(0)
{
  // increment the currentIndex to start at the first valid index
  // for this ParallelDescriptor::MyProc

  int myproc = ParallelDescriptor::MyProc();
  while(currentIndex < fabArray.length() &&
        fabArray.DistributionMap().ProcessorMap()[currentIndex] != myproc) {
    ++currentIndex;
  }
}



template <class T, class FAB>
FabArrayIterator<T, FAB>::~FabArrayIterator() {
}


template <class T, class FAB>
FAB &FabArrayIterator<T, FAB>::operator()() {
  int myproc = ParallelDescriptor::MyProc();
  assert(fabArray.DistributionMap().ProcessorMap()[currentIndex] == myproc);
  return fabArray[currentIndex];
}


template <class T, class FAB>
const FAB &FabArrayIterator<T, FAB>::operator()() const {
  int myproc = ParallelDescriptor::MyProc();
  assert(fabArray.DistributionMap().ProcessorMap()[currentIndex] == myproc);
  return fabArray[currentIndex];
}


template <class T, class FAB>
FabArrayIterator<T, FAB> &FabArrayIterator<T, FAB>::operator++() {
  int myproc = ParallelDescriptor::MyProc();
  // go to the next index on this processor
  do {
    ++currentIndex;
  } while(currentIndex < fabArray.length() &&
	  fabArray.DistributionMap().ProcessorMap()[currentIndex] != myproc);

  return *this;
}



template <class T, class FAB>
int FabArrayIterator<T, FAB>::index() const {
  return currentIndex;
}


template <class T, class FAB>
const Box &FabArrayIterator<T, FAB>::validbox() const {
  return fabArray.box(currentIndex);
}


template <class T, class FAB>
const Box &FabArrayIterator<T, FAB>::fabbox() const {
  return fabArray[currentIndex].box();
}


template <class T, class FAB>
bool FabArrayIterator<T, FAB>::isValid() {
  bool valid = (currentIndex >= 0 && currentIndex < fabArray.length());
  if( ! valid) {
    //ParallelDescriptor::Synchronize("FAI::isValid() sync");
    ParallelDescriptor::Synchronize();
  }
  return valid;
}




// --------------------------
template <class T, class FAB>
DependentFabArrayIterator<T, FAB>::DependentFabArrayIterator(
	                    FabArrayIterator<T, FAB> &controllerfaiter,
	                    FabArray<T, FAB> &dependentfabarray)
                         : fabArray(dependentfabarray),
			   controller(&controllerfaiter)
{
}


template <class T, class FAB>
DependentFabArrayIterator<T, FAB>::DependentFabArrayIterator(
	                    FabArrayIterator<T, FAB> &controllerfaiter,
	                    const FabArray<T, FAB> &dependentfabarray)
                         : fabArray((FabArray<T, FAB> &) dependentfabarray),
			   controller(&controllerfaiter)
{
}


template <class T, class FAB>
DependentFabArrayIterator<T, FAB>::~DependentFabArrayIterator() {
}


template <class T, class FAB>
int DependentFabArrayIterator<T, FAB>::index() const {
  return controller->index();
}


template <class T, class FAB>
const Box &DependentFabArrayIterator<T, FAB>::validbox() const {
  return fabArray.box(controller->index());
}


template <class T, class FAB>
const Box &DependentFabArrayIterator<T, FAB>::fabbox() const {
  return fabArray[controller->index()].box();
}


template <class T, class FAB>
FAB &DependentFabArrayIterator<T, FAB>::operator()() {
  return fabArray[controller->index()];
}


template <class T, class FAB>
const FAB &DependentFabArrayIterator<T, FAB>::operator()() const {
  return fabArray[controller->index()];
}


// --------------------------
// --------------------------
template <class T, class FAB>
ConstFabArrayIterator<T, FAB>::ConstFabArrayIterator(
				   const FabArray<T, FAB> &fabarray,
				   const int nghost)
                         : fabArray(fabarray),
			   nGhost(nghost),
			   currentIndex(0)
{
  // increment the currentIndex to start at the first valid index
  // for this ParallelDescriptor::MyProc

  int myproc = ParallelDescriptor::MyProc();
  while(currentIndex < fabArray.length() &&
        fabArray.DistributionMap()[currentIndex] != myproc) {
    ++currentIndex;
  }
}


template <class T, class FAB>
ConstFabArrayIterator<T, FAB>::~ConstFabArrayIterator() {
}


template <class T, class FAB>
const FAB &ConstFabArrayIterator<T, FAB>::operator()() const {
  return fabArray[currentIndex];
}


template <class T, class FAB>
ConstFabArrayIterator<T, FAB> &ConstFabArrayIterator<T, FAB>::operator++() {
  int myproc = ParallelDescriptor::MyProc();
  // go to the next index on this processor
  do {
    ++currentIndex;
  } while(currentIndex < fabArray.length() &&
	  fabArray.DistributionMap()[currentIndex] != myproc);

  return *this;
}

template <class T, class FAB>
int ConstFabArrayIterator<T, FAB>::index() const {
  return (currentIndex);
}


template <class T, class FAB>
const Box &ConstFabArrayIterator<T, FAB>::validbox() const {
  return fabArray.box(currentIndex);
}


template <class T, class FAB>
const Box &ConstFabArrayIterator<T, FAB>::fabbox() const {
  return fabArray[currentIndex].box();
}


template <class T, class FAB>
bool ConstFabArrayIterator<T, FAB>::isValid() const {
  bool valid = (currentIndex >= 0 && currentIndex < fabArray.length());
  if( ! valid) {
    //ParallelDescriptor::Synchronize("CFAI::isValid sync");
    ParallelDescriptor::Synchronize();
  }
  return valid;
}


// --------------------------
template <class T, class FAB>
ConstDependentFabArrayIterator<T, FAB>::ConstDependentFabArrayIterator(
	                    ConstFabArrayIterator<T, FAB> &controllerfaiter,
	                    const FabArray<T, FAB> &dependentfabarray)
                         : fabArray(dependentfabarray),
			   controller(&controllerfaiter)
{
}


template <class T, class FAB>
ConstDependentFabArrayIterator<T, FAB>::~ConstDependentFabArrayIterator() {
}


template <class T, class FAB>
int ConstDependentFabArrayIterator<T, FAB>::index() const {
  return (controller->index());
}

template <class T, class FAB>
const Box &ConstDependentFabArrayIterator<T, FAB>::validbox() const {
  return fabArray.box(controller->index());
}


template <class T, class FAB>
const Box &ConstDependentFabArrayIterator<T, FAB>::fabbox() const {
  return fabArray[controller->index()].box();
}


template <class T, class FAB>
const FAB &ConstDependentFabArrayIterator<T, FAB>::operator()() const {
  return fabArray[controller->index()];
}
// ========================================================================

template <class T, class FAB>
FabArrayCopyDescriptor<T, FAB>::FabArrayCopyDescriptor(bool cacheremotedata)
                         : nextFabArrayId(0),
                           nextFillBoxId(0),
			   dataAvailable(true),
			   cacheRemoteData(cacheremotedata),
			   completelyFilled(false),
			   totalRemoteBoxes(0),
			   //nSubBoxes(0),
			   totalRemoteBytes(0)
{
  if(cacheRemoteData == false) {
    cerr << "Error in FabArrayCopyDescriptor(bool):  "
	 << "uncached remote data not implemented." << endl;
    ParallelDescriptor::Abort("Exiting.");
  }

  //nRemoteBoxes.resize(fabArray.length(), 0);
  //nRemoteBytes.resize(fabArray.length(), 0);

  FabComTag fabComTag;
  int tagSize = sizeof(FabComTag);
  ParallelDescriptor::SetMessageHeaderSize(tagSize);

  fabArrays.reserve(16);
  fabComTagList.clear();

}  // end FabArrayCopyDescriptor(...)



// ----------------------------------------------------------------------------
template <class T, class FAB>
FabArrayId
FabArrayCopyDescriptor<T, FAB>::RegisterFabArray(FabArray<T, FAB> *fabarray) {
  assert(nextFabArrayId == fabArrays.length());
  fabArrays.resize(nextFabArrayId + 1);
  fabArrays[nextFabArrayId] = fabarray;
  return FabArrayId(nextFabArrayId++);
}



// ----------------------------------------------------------------------------
template <class T, class FAB>
FillBoxId
FabArrayCopyDescriptor<T, FAB>::AddBox(const FabArrayId &fabarrayid,
			               const Box &destFabBox,
				       BoxList &returnedUnfilledBoxes)
{
  return AddBox(fabarrayid, destFabBox, returnedUnfilledBoxes,
		0, 0, fabArrays[fabarrayid.Id()]->nComp());
}


// ----------------------------------------------------------------------------
template <class T, class FAB>
FillBoxId
FabArrayCopyDescriptor<T, FAB>::AddBox(const FabArrayId &fabarrayid,
			               const Box &destFabBox,
				       BoxList &returnedUnfilledBoxes,
				       int srccomp, int destcomp, int numcomp)
{
  int myproc = ParallelDescriptor::MyProc();
  FabComTag fabComTag;

  //nRemoteBoxes.resize(fabArray.length(), 0);
  //nRemoteBytes.resize(fabArray.length(), 0);

  BoxDomain unfilledBoxDomain(destFabBox.ixType());
  BoxList filledBoxes(destFabBox.ixType());
  unfilledBoxDomain.add(destFabBox);
  FabArray<T, FAB> *fabArray = fabArrays[fabarrayid.Id()];
  for(int i = 0; i < fabArray->length(); ++i) {
    const Box intersectBox = destFabBox & fabArray->box(i);
    if(intersectBox.isValid()) {
      filledBoxes.add(intersectBox);
      int remoteProc = fabArray->DistributionMap().ProcessorMap()[i];
      FabCopyDescriptor *fcd = new FabCopyDescriptor;
      fcd->fabArrayId = fabarrayid.Id();
      fcd->fillBoxId  = nextFillBoxId;
      fcd->subBox = intersectBox;
      fcd->myProc = myproc;
      fcd->copyFromProc = remoteProc;
      fcd->copyFromIndex = i;
      fcd->srcComp  = srccomp;
      fcd->destComp = destcomp;
      fcd->nComp    = numcomp;
      if(myproc == remoteProc) {  // data is local
        fcd->fillType = FillLocally;
	fcd->localFabSource = &((*fabArray)[i]);
      } else {                    // data is remote
	dataAvailable = false;
        fcd->fillType = FillRemotely;
	if(cacheRemoteData) {
	  fcd->localFabSource = new FAB(intersectBox, numcomp);
	  //fcd->localFabSource->setVal(1.e28);
	  fcd->cacheDataAllocated = true;
	  // send request to send data
	  fabComTag.fabArrayId = fabarrayid.Id();
	  fabComTag.fillBoxId  = nextFillBoxId;
	  fabComTag.fabIndex   = i;
	  fabComTag.procThatNeedsData = myproc;
	  fabComTag.procThatHasData   = remoteProc;
	  fabComTag.subBox     = intersectBox;
	  fabComTag.srcComp    = srccomp;
	  fabComTag.destComp   = destcomp;
	  fabComTag.nComp      = numcomp;

	  // dont SendData yet
	  fabComTagList.append(fabComTag);
	}
        //nRemoteBoxes[i]++;
        //nRemoteBytes[i]  += intersectBox.numPts() * fabArray->nComp() * sizeof(T);
        totalRemoteBoxes++;
        totalRemoteBytes += intersectBox.numPts() * fabArray->nComp() * sizeof(T);
      }

      fabCopyDescList.append(fcd);
      unfilledBoxDomain.rmBox(intersectBox);
      //nSubBoxes++;
    }  // end if(intersectBox...)
  }
  returnedUnfilledBoxes.clear();
  returnedUnfilledBoxes = unfilledBoxDomain.boxList();
  //if(unfilledBoxes.length() == 0) {
    //completelyFilled = true;
  //}

  return FillBoxId(nextFillBoxId++, destFabBox, filledBoxes);
}  // end AddBox



// ----------------------------------------------------------------------------
template <class T, class FAB>
FillBoxId
FabArrayCopyDescriptor<T, FAB>::AddBox(const FabArrayId &fabarrayid,
			               const Box &destFabBox,
				       BoxList &returnedUnfilledBoxes,
				       int fabarrayindex,
				       int srccomp, int destcomp, int numcomp)
{

  int myproc = ParallelDescriptor::MyProc();
  FabComTag fabComTag;

  BoxDomain unfilledBoxDomain(destFabBox.ixType());
  BoxList filledBoxes(destFabBox.ixType());
  unfilledBoxDomain.add(destFabBox);
  FabArray<T, FAB> *fabArray = fabArrays[fabarrayid.Id()];
  assert(fabarrayindex >= 0 && fabarrayindex < fabArray->length());
  //for(int i = 0; i < fabArray->length(); ++i) {
  int i = fabarrayindex;
    const Box &intersectBox = destFabBox & fabArray->box(i);
    if(intersectBox.isValid()) {
      filledBoxes.add(intersectBox);
      int remoteProc = fabArray->DistributionMap().ProcessorMap()[i];
      FabCopyDescriptor *fcd = new FabCopyDescriptor;
      fcd->fabArrayId = fabarrayid.Id();
      fcd->fillBoxId  = nextFillBoxId;
      fcd->subBox = intersectBox;
      fcd->myProc = myproc;
      fcd->copyFromProc = remoteProc;
      fcd->copyFromIndex = i;
      fcd->srcComp  = srccomp;
      fcd->destComp = destcomp;
      fcd->nComp    = numcomp;
      if(myproc == remoteProc) {  // data is local
        fcd->fillType = FillLocally;
	fcd->localFabSource = &((*fabArray)[i]);
      } else {                    // data is remote
	dataAvailable = false;
        fcd->fillType = FillRemotely;
	if(cacheRemoteData) {
	  fcd->localFabSource = new FAB(intersectBox, numcomp);
	  //fcd->localFabSource->setVal(1.e28);
	  fcd->cacheDataAllocated = true;
	  // send request to send data
	  fabComTag.fabArrayId = fabarrayid.Id();
	  fabComTag.fillBoxId  = nextFillBoxId;
	  fabComTag.fabIndex   = i;
	  fabComTag.procThatNeedsData = myproc;
	  fabComTag.procThatHasData   = remoteProc;
	  fabComTag.subBox     = intersectBox;
	  fabComTag.srcComp    = srccomp;
	  fabComTag.destComp   = destcomp;
	  fabComTag.nComp      = numcomp;

	  // dont SendData yet
	  fabComTagList.append(fabComTag);
	}
        //nRemoteBoxes[i]++;
        //nRemoteBytes[i]  += intersectBox.numPts() * fabArray->nComp() * sizeof(T);
        totalRemoteBoxes++;
        totalRemoteBytes += intersectBox.numPts() * fabArray->nComp() * sizeof(T);
      }

      fabCopyDescList.append(fcd);
      unfilledBoxDomain.rmBox(intersectBox);
      //nSubBoxes++;
    }  // end if(intersectBox...)
  //}
  returnedUnfilledBoxes.clear();
  returnedUnfilledBoxes = unfilledBoxDomain.boxList();
  //if(unfilledBoxes.length() == 0) {
    //completelyFilled = true;
  //}

  return FillBoxId(nextFillBoxId++, destFabBox, filledBoxes);
}  // end AddBox



// ----------------------------------------------------------------------------
template <class T, class FAB>
FabArrayCopyDescriptor<T, FAB>::~FabArrayCopyDescriptor() {
  for(ListIterator<FabCopyDescriptor *> fli(fabCopyDescList); fli; ++fli) {
    delete fli();
  }
}


// ----------------------------------------------------------------------------
template <class T, class FAB>
void FabArrayCopyDescriptor<T, FAB>::CollectData() {

  int dataWaitingSize;
  int myproc = ParallelDescriptor::MyProc();
  FabComTag fabComTag;

  int tagSize = sizeof(FabComTag);
  ParallelDescriptor::SetMessageHeaderSize(tagSize);
  int *nullptr = NULL;

  // go through the fabComTagList and send all fab requests
  if(ParallelDescriptor::NProcs() == 1) {
    assert(fabComTagList.length() == 0);
  }
  for(ListIterator<FabComTag> fctli(fabComTagList); fctli; ++fctli) {
    ParallelDescriptor::SendData(fctli().procThatHasData, fctli(), NULL, 0);
  }
  fabComTagList.clear();

  ParallelDescriptor::Synchronize();

  // check for consistency and correct number of received and expected messages
  if( ! dataAvailable) {

    while(ParallelDescriptor::GetMessageHeader(dataWaitingSize, &fabComTag))
    {  // data was sent to this processor

      // checks
      if(dataWaitingSize != 0) {
        cerr << "Error in FabArrayCopyDescriptor<T, FAB>::CollectData:  ";
        cerr << "data payload size for data send request is nonzero." << endl;
        cerr << "data payload size = " << dataWaitingSize << endl;
        ParallelDescriptor::Abort("CollectData:  bad data send request size");
      }

      if(fabComTag.procThatHasData != myproc) {
	cerr << "Error 1 in CollectData..." << endl;
	ParallelDescriptor::Abort("Error 1 in CollectData...");
      }
      if( ! fabComTag.subBox.ok()) {
	cerr << "Error 2 in CollectData..." << endl;
	ParallelDescriptor::Abort("Error 2 in CollectData...");
      }
      if( ! fabArrays[fabComTag.fabArrayId]->
	      box(fabComTag.fabIndex).contains(fabComTag.subBox))
      {
	cerr << "Error 3 in CollectData..." << endl;
	ParallelDescriptor::Abort("Error 3 in CollectData...");
      }
      // other checks here

      ParallelDescriptor::ReceiveData(nullptr, 0);  // to advance msg header

      // now send the requested data from the fabs on this processor

      FAB tempFab(fabComTag.subBox, fabComTag.nComp);
      //tempFab.setVal(1.e29);
      tempFab.copy((*fabArrays[fabComTag.fabArrayId])[fabComTag.fabIndex],
		   fabComTag.subBox, fabComTag.srcComp,
		   fabComTag.subBox, 0, fabComTag.nComp);
				  // ^ copy to zero component of temp fab
      ParallelDescriptor::SendData(fabComTag.procThatNeedsData, fabComTag,
				   tempFab.dataPtr(),
				   tempFab.box().numPts() *
				   tempFab.nComp() * sizeof(T));

    }
  }

  ParallelDescriptor::Synchronize();

  // now collect all remote data into local fab data caches
  while(ParallelDescriptor::GetMessageHeader(dataWaitingSize, &fabComTag))
  {  // data was sent to this processor

    // check for consistency and correct number of received and expected messages
    if(myproc != fabComTag.procThatNeedsData) {
      ParallelDescriptor::Abort("FillFab:  bad procThatNeedsData");
    }

    // find the box in the list and move the data
    bool matchFound = false;
    for(ListIterator<FabCopyDescriptor *> fli(fabCopyDescList);
	fli && ! matchFound;
	++fli)
    {
      if(fli()->subBox == fabComTag.subBox &&
	 fabComTag.fabArrayId == fli()->fabArrayId &&
	 fabComTag.fillBoxId == fli()->fillBoxId)
      {
	// copy the data
        ParallelDescriptor::ReceiveData(fli()->localFabSource->dataPtr(),
	         fabComTag.subBox.numPts() * fabComTag.nComp * sizeof(T));
	matchFound = true;
      }
    }
    if(matchFound == false) {
      ParallelDescriptor::Abort("FillFab:  match not found");
    }
  }

  dataAvailable = true;

  ParallelDescriptor::Synchronize();

}  // end CollectData()


// ----------------------------------------------------------------------------
template <class T, class FAB>
void FabArrayCopyDescriptor<T, FAB>::FillFab(const FabArrayId &fabarrayid,
		                             const FillBoxId  &fillboxid,
		                             FAB &destFab)
{
  // at this point, all remote data should be in local caches

  if( ! dataAvailable) {
    cerr << "Error in FabArrayCopyDescriptor<T, FAB>::FillFab:  ";
    cerr << "data not available" << endl;
    ParallelDescriptor::Abort("FillFab:  data not available");
  }

  // optimize this vvvvvvvvvvvvvv
  for(ListIterator<FabCopyDescriptor *> fli(fabCopyDescList); fli; ++fli) {
    if(fli()->fabArrayId == fabarrayid.Id() && fli()->fillBoxId == fillboxid.Id())
    {
      int localSrcComp;
      if(fli()->fillType == FillLocally) {
        localSrcComp = fli()->srcComp;
      } else {  // FillRemotely
        localSrcComp = 0;  // copy from zero component of local temp fab
      }
      destFab.copy(*(fli()->localFabSource), fli()->subBox, localSrcComp,
		   fli()->subBox, fli()->destComp, fli()->nComp);
    }
  }
}  // end FillFab


// ----------------------------------------------------------------------------
template <class T, class FAB>
void FabArrayCopyDescriptor<T, FAB>::PrintStats() {
  int myproc = ParallelDescriptor::MyProc();

  cout << "----- " << myproc << ":  "
       << "Parallel stats for FabArrayCopyDescriptor:" << endl;
  cout << "----- " << myproc << ":  "
       << "    totalRemoteBoxes = " << totalRemoteBoxes << endl;
  cout << "----- " << myproc << ":  "
       << "    totalRemoteBytes = " << totalRemoteBytes << endl;
  for(int fa = 0; fa < fabArrays.length(); ++fa) {
    cout << "fabArrays[" << fa << "]->boxArray() = " << fabArrays[fa]->boxArray()
	 << endl;
  }


  /*
  int i;
  cout << "           nSubBoxes = " << nSubBoxes << endl;
  for(i = 0; i < nRemoteBoxes.length(); ++i) {
    cout << "    nRemoteBoxes[" << i << "] = " << nRemoteBoxes[i] << endl;
  }
  for(i = 0; i < nRemoteBytes.length(); ++i) {
    cout << "    nRemoteBytes[" << i << "] = " << nRemoteBytes[i] << endl;
  }
  */
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ========================================================================
