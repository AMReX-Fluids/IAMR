// --------------------------------------------------------------------
// FabSet.C
// --------------------------------------------------------------------
#include <FabSet.H>
#include <Looping.H>

// --------------------------------------------------------------------
FabSet::FabSet() : MultiFab() {
}


// --------------------------------------------------------------------
FabSet::~FabSet() {
}


// --------------------------------------------------------------------
FabSet::FabSet(int _len) : MultiFab() {
  fabparray.resize(_len);
}


// --------------------------------------------------------------------
FabSet::FabSet(istream &is) : MultiFab() {
    readFrom(is);
}


// --------------------------------------------------------------------
void FabSet::setFab(int boxno, FArrayBox *fab) {
      if(n_comp == 0) {
        n_comp = fab->nComp();
      }
      assert(n_comp == fab->nComp());
      assert(boxarray.ready());
      assert( ! fabparray.defined(boxno));
      if(distributionMap.ProcessorMap()[boxno] == ParallelDescriptor::MyProc()) {
        fabparray.set(boxno, fab);
      } else {
        ParallelDescriptor::Abort("Error: FabArray<T, FAB>::setFab: nonlocal set");
      }
      fabboxarray.convert(fab->box().ixType());
      fabboxarray.set(boxno, fab->box());
}


// --------------------------------------------------------------------
const FabSet &FabSet::copyTo(FArrayBox &dest) const {
    this->copy(dest);
    return *this;
}


// --------------------------------------------------------------------
const FabSet &FabSet::copyTo(FArrayBox& dest, int src_comp,
			     int dest_comp, int num_comp) const
{
    this->copy(dest, src_comp, dest_comp, num_comp);
    return *this;
}


// --------------------------------------------------------------------
const FabSet &FabSet::copyTo(FArrayBox &dest, const Box &subbox,
			     int src_comp, int dest_comp,
			     int num_comp) const
{
    this->copy(dest, subbox, src_comp, dest_comp, num_comp);
    return *this;
}

// --------------------------------------------------------------------
FabSet &FabSet::copyFrom(const FArrayBox &src) {
    for(FabSetIterator fsi(*this); fsi.isValid(); ++fsi) {
	fsi().copy(src);
    }
    return *this;
}

// --------------------------------------------------------------------
FabSet &FabSet::copyFrom(const FArrayBox& src, int src_comp,
			 int dest_comp, int num_comp)
{
    for(FabSetIterator fsi(*this); fsi.isValid(); ++fsi) {
	fsi().copy(src,src_comp,dest_comp,num_comp);
    }
    return *this;
}


// --------------------------------------------------------------------
FabSet &FabSet::copyFrom(const FArrayBox &src, const Box &subbox,
		         int src_comp, int dest_comp, int num_comp)
{
    const Box &sbox = src.box();
    assert( sbox.contains(subbox) );
    for(FabSetIterator fsi(*this); fsi.isValid(); ++fsi) {
	FArrayBox &fab = fsi();
	Box dbox = fab.box();
	dbox &= subbox;
	if(dbox.ok()) {
	    fsi().copy(src,dbox,src_comp,dbox,dest_comp,num_comp);
	}
    }
    return *this;
}


// --------------------------------------------------------------------
// the following are different from MultiFab only in the return value

// --------------------------------------------------------------------
FabSet &FabSet::plus(Real v, int comp, int num_comp) {
    this->plus(v, comp, num_comp);
    return *this;
}


// --------------------------------------------------------------------
FabSet &FabSet::plus(Real v, const Box &subreg, int comp, int num_comp)
{
    this->plus(v, subreg, comp, num_comp);
    return *this;
}


// --------------------------------------------------------------------
FabSet &FabSet::mult(Real v, int comp, int num_comp) {
    this->mult(v, comp, num_comp);
    return *this;
}


// --------------------------------------------------------------------
FabSet &FabSet::mult(Real v, const Box &subreg, int comp, int num_comp)
{
    this->mult(v, subreg, comp, num_comp);
    return *this;
}

// --------------------------------------------------------------------
FabSet &FabSet::copyFrom(const FabSet &src) {
    this->copy(src);
    return *this;
}

// --------------------------------------------------------------------
FabSet &FabSet::copyFrom(const FabSet &src, int src_comp, int dest_comp,
		         int num_comp)
{
    this->copy(src, src_comp, dest_comp, num_comp);
    return *this;
}

// --------------------------------------------------------------------
FabSet &FabSet::copyFrom(const FabSet &src, const Box &subreg,
		         int src_comp, int dest_comp, int num_comp)
{
cerr << "FabSet::copyFrom(FabSet, Box, ...) not yet implemented" << endl;
ParallelDescriptor::Abort("Exiting.");
    return *this;
}


// --------------------------------------------------------------------
FabSet &FabSet::copyFrom(const MultiFab &src, int nghost, int src_comp,
		         int dest_comp, int num_comp)
{
cerr << "FabSet::copyTo(MultiFab, nghost, ...) not yet implemented" << endl;
ParallelDescriptor::Abort("Exiting.");
    return *this;
}

// --------------------------------------------------------------------
const FabSet &FabSet::copyTo(MultiFab &dest, int nghost, int src_comp,
	                     int dest_comp, int num_comp) const
{
cerr << "FabSet::copyFrom(MultiFab, nghost, ...) not yet implemented" << endl;
ParallelDescriptor::Abort("Exiting.");
    return *this;
}

// --------------------------------------------------------------------
FabSet &FabSet::plusFrom(const MultiFab &src, int nghost, int src_comp,
		         int dest_comp, int num_comp)
{
cerr << "FabSet::copyFrom(MultiFab, nghost, ...) not yet implemented" << endl;
ParallelDescriptor::Abort("Exiting.");
/*
    int slen = src.length();
    int dlen = length();
    assert (nghost <= src.nGrow());
    const BoxArray& sba = src.boxArray();

 // turn this loop inside out for parallel (mfiterate over *this)
 // to fill dest locally

    for(ConstMultiFabIterator cmfi(src); cmfi.isValid(); ++cmfi) {
	const FArrayBox& sfab = cmfi();
	Box sbox = grow(cmfi.validbox(), nghost);
	for (int d = 0; d < dlen; d++) {
	    FArrayBox& dfab = (*this)[d];
	    Box ovlp = dfab.box();
	    ovlp &= sbox;
	    if (ovlp.ok()) {
		dfab.plus(sfab,ovlp,src_comp,dest_comp,num_comp);
	    }
	}
    }
*/
    return *this;
}

// --------------------------------------------------------------------
// Linear combination this := a*this + b*src
// Note: corresponding fabsets must be commensurate
FabSet &FabSet::linComb(Real a, Real b, const FabSet &src, int src_comp,
		        int dest_comp, int num_comp)
{
    assert(length() == src.length());
    for(FabSetIterator fsi(*this); fsi.isValid(); ++fsi) {
      DependentFabSetIterator dfsi(fsi, src);
	FArrayBox &dfab = fsi();
	const FArrayBox &sfab = dfsi();
	const Box &dbox = dfab.box();
	const Box &sbox = sfab.box();
	assert( dbox == sbox );
	  // WARNING: same fab used as src and dest here
	dfab.linComb(dfab,dbox,dest_comp,sfab,sbox,src_comp,
		     a,b,dbox,dest_comp,num_comp);
    }
    return *this;
}

// --------------------------------------------------------------------
FabSet &FabSet::linComb(Real a, const MultiFab &mfa, int a_comp,
		        Real b, const MultiFab &mfb, int b_comp,
		        int dest_comp, int num_comp, int n_ghost)
{
  if(ParallelDescriptor::NProcs() > 1) {
    cerr << "FabSet::linComb(2) not yet implemented" << endl;
    ParallelDescriptor::Abort("Exiting.");
  }

    const BoxArray &bxa = mfa.boxArray();
    const BoxArray &bxb = mfb.boxArray();
    assert( bxa == bxb);

 // turn this loop inside out for parallel (mfiterate over *this)
 // to fill dest locally

    assert( n_ghost <= mfa.nGrow());
    assert( n_ghost <= mfb.nGrow());

    int nfab = bxa.length();
    int nreg = length();

    for (int grd = 0; grd < nfab; grd++) {
        const BOX& grd_box = grow(bxa[grd],n_ghost);
        const FARRAYBOX& a_fab = mfa[grd];
        const FARRAYBOX& b_fab = mfb[grd];
        for (int reg = 0; reg < nreg; reg++) {
            FARRAYBOX &reg_fab = (*this)[reg];
            BOX ovlp(reg_fab.box());
            ovlp &= grd_box;
            if (ovlp.ok()) {
                reg_fab.linComb(a_fab,ovlp,a_comp,b_fab,ovlp,b_comp,
                                a,b,ovlp,dest_comp,num_comp);
            }
        }
    }

    return *this;
}

// --------------------------------------------------------------------
ostream &
operator << (ostream &os, const FabSet &mf)
{
cerr << "FabSet::operator<<() not yet implemented" << endl;
ParallelDescriptor::Abort("Exiting.");
/*
    int nfab = mf.length();
    os << "(FabSet "
       << nfab << '\n';
    for (int i = 0; i < nfab; i++) {
	os << mf[i] << '\n';
    }
    os << ")" << flush;
*/
    return os;
}

// --------------------------------------------------------------------
ostream &
FabSet::writeOn(ostream &os) const
{
cerr << "FabSet::writeOn() not yet implemented" << endl;
ParallelDescriptor::Abort("Exiting.");
/*
    assert( ready() );
    int nfab = length();
    os << nfab << '\n';
    for (int i = 0; i < nfab; i++) {
	get(i).writeOn(os);
    }
*/
    return os;
}

// --------------------------------------------------------------------
istream &
FabSet::readFrom(istream &is)
{
cerr << "FabSet::readFrom() not yet implemented" << endl;
ParallelDescriptor::Abort("Exiting.");
/*
    if (ready()) {
	clear();
    }
    int ngrd;
    is >> ngrd;
    while (is.get() != '\n');
    resize(ngrd);
    for (int i = 0; i < ngrd; i++) {
	FArrayBox* tmp = new FArrayBox;
	tmp->readFrom(is);
	set(i,tmp);
    }
*/
    return is;
}



// --------------------------------------------------------------------
FabSetIterator::FabSetIterator(FabSet &fabset)
               : MultiFabIterator(fabset)
{
}

// --------------------------------------------------------------------
FabSetIterator::~FabSetIterator() {
}


// --------------------------------------------------------------------
DependentFabSetIterator::DependentFabSetIterator(FabSetIterator &controllerfsiter,
                                                 FabSet &fabset)
                        : DependentMultiFabIterator(controllerfsiter, fabset)
{
}

// --------------------------------------------------------------------
DependentFabSetIterator::DependentFabSetIterator(FabSetIterator &controllerfsiter,
                                                 const FabSet &fabset)
                        : DependentMultiFabIterator(controllerfsiter, fabset)
{
}

// --------------------------------------------------------------------
DependentFabSetIterator::DependentFabSetIterator(MultiFabIterator &controllerfsiter,
                                                 FabSet &fabset)
                        : DependentMultiFabIterator(controllerfsiter, fabset)
{
}

// --------------------------------------------------------------------
DependentFabSetIterator::DependentFabSetIterator(MultiFabIterator &controllerfsiter,
                                                 const FabSet &fabset)
                        : DependentMultiFabIterator(controllerfsiter, fabset)
{
}

// --------------------------------------------------------------------
DependentFabSetIterator::~DependentFabSetIterator() {
}


// --------------------------------------------------------------------
ConstFabSetIterator::ConstFabSetIterator(const FabSet &fabset)
                    : ConstMultiFabIterator(fabset)
{
}

// --------------------------------------------------------------------
ConstFabSetIterator::~ConstFabSetIterator() {
}


// --------------------------------------------------------------------
ConstDependentFabSetIterator::ConstDependentFabSetIterator(
                                 ConstFabSetIterator &controllerfsiter,
                                 const FabSet &fabset)
                : ConstDependentMultiFabIterator(controllerfsiter, fabset)
{
}

// --------------------------------------------------------------------
ConstDependentFabSetIterator::ConstDependentFabSetIterator(
                                 ConstMultiFabIterator &controllerfsiter,
                                 const FabSet &fabset)
                : ConstDependentMultiFabIterator(controllerfsiter, fabset)
{
}

// --------------------------------------------------------------------
ConstDependentFabSetIterator::~ConstDependentFabSetIterator() {
}



// --------------------------------------------------------------------
FabSetCopyDescriptor::FabSetCopyDescriptor(bool cacheremotedata)
                     : MultiFabCopyDescriptor(cacheremotedata)
{
}

// --------------------------------------------------------------------
FabSetCopyDescriptor::~FabSetCopyDescriptor() {
}

// --------------------------------------------------------------------
// --------------------------------------------------------------------
