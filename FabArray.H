//BL_COPYRIGHT_NOTICE

#ifndef _FABARRAY_H_
#define _FABARRAY_H_

//
// $Id: FabArray.H,v 1.1 1997-07-08 23:08:06 vince Exp $
//

#include <Assert.H>
#include <aString.H>
#include <PArray.H>
#include <Array.H>

#include <Boolean.H>
#include <Box.H>
#include <BoxLib.H>
#include <BoxArray.H>
#include <FArrayBox.H>
#include <DistributionMapping.H>
#include <ParallelDescriptor.H>

/*@ManDoc: An enumumeration that controls whether or not the memory for a FAB
           will actually be allocated on construction of a FabArray.
           Possible values are: Fab\_noallocate and Fab\_allocate.
*/

enum FabAlloc
{
    Fab_noallocate = 0,
    Fab_allocate
};

//
//@Man:
//@Memo: A Collection of Fortran Array-like Objects
/*@Doc:

  The FabArray<T,FAB> class implements a collection (stored as an array) of
  Fortran array-like objects.  The parameterized type FAB is intended to be
  any class derived from BaseFab<T>.  For example, FAB may be a BaseFab of
  integers, so we could write:

    FabArray< int,BaseFab<int> > int\_fabs;

  Then int\_fabs is a FabArray that can hold a collection of BaseFab<int>
  objects.

  FabArray is not just a general container class for Fortran arrays.  It is
  intended to hold "grid" data for use in finite difference calculations in
  which the data is defined on a union of (usually disjoint) rectangular
  regions embedded in a uniform index space.  This region, called the valid
  region, is represented by a BoxArray.  For the purposes of this discussion,
  the Kth Box in the BoxArray represents the interior region of the Kth grid.

  Since the intent is to be used with finite difference calculations a
  FabArray also includes the notion of a boundary region for each grid.  The
  boundary region is specified by the ngrow parameter which tells the FabArray
  to allocate each FAB to be ngrow cells larger in all directions than the
  underlying Box.  The larger region covered by the union of all the FABs is
  called the region of definition.  The underlying notion is that the valid
  region contains the grid interior data and the region of definition includes
  the interior region plus the boundary areas.

  Operations are available to copy data from the valid regions into these
  boundary areas where the two overlap.  The number of components, that is,
  the number of values that can be stored in each cell of a FAB, is either
  given as an argument to the constructor or is inherent in the definition of
  the underlying FAB.  Each FAB in the FabArray will have the same number of
  components.

  In summary, a FabArray is an array of FABs.  The Kth element contains a FAB
  that holds the data for the Kth grid, a Box that defines the valid region
  of the Kth grid.

  A typical use for a FabArray would be to hold the solution vector or
  right-hand-side when solving a linear system of equations on a union of
  rectangular grids.  The copy operations would be used to copy data from the
  valid regions of neighboring grids into the boundary regions after each
  relaxation step of the iterative method.  If a multigrid method is used, a
  FabArray could be used to hold the data at each level in the multigrid
  hierarchy.

  This class is a concrete class not a polymorphic one.

  This class does NOT provide a copy constructor or assignment operator.
*/


template <class T, class FAB>
class FabArray
{

public:
    //
    //@ManDoc: Constructs an empty FabArray<T,FAB>.
    //
    FabArray ();

    /*@ManDoc: Construct a FabArray<T,FAB> with a valid region defined by bxs
               and a region of definition defined by the grow factor ngrow.
               If mem\_mode is defined to be Fab\_allocate then FABs are
               allocated for each Box in the BoxArray.  The size of the Kth
               FAB is given by bxs[K] grown by ngrow.  If mem\_mode is defined
               to be Fab\_noallocate, then no FABs are allocated at this time,
               but can be defined later.  The number of components in each
               FAB is not specified and is expected to be implicit in the
               definition of the FAB class.  That is, the FAB constructor will
               take only a Box argument.  Call this constructor number two.
    */
    FabArray (const BoxArray& bxs,
              int             ngrow,
              FabAlloc        mem_mode = Fab_allocate);

    /*@ManDoc: This constructor is similar to constructor number two except
               that it also specifies the number of components for each FAB.
               The FAB constructor used here will take both a Box and nvar
               (the number of components) as arguments.
    */
public:
    FabArray (const BoxArray& bxs,
              int             nvar,
              int             ngrow,
              FabAlloc        mem_mode = Fab_allocate);


public:
    //
    //@ManDoc: The (virtual) destructor -- deletes all FABs in the array.
    //
    virtual ~FabArray ();

    /*@ManDoc: Define this FabArray identically to that performed by
               the constructor having an analogous function signature.
               This is only valid if this FabArray was defined using
               the default constructor.
    */
    void define (const BoxArray& bxs,
                 int             ngrow,
                 FabAlloc        mem_mode);

    /*@ManDoc: Define this FabArray identically to that performed by
               the constructor having an analogous function signature.
               This is only valid if this FabArray was defined using
               the default constructor.
    */
    void define (const BoxArray& bxs,
                 int             nvar,
                 int             ngrow,
                 FabAlloc        mem_mode);

    /*@ManDoc: Returns true if the FabArray is well-defined.  That is,
               if FABs are allocated for each Box in the BoxArray and the
               sizes of the FABs and the number of components are consistent
               with the definition of the FabArray.
    */
    bool ok () const;

    //
    //@ManDoc: Returns the grow factor that defines the region of definition.
    //
    inline int nGrow () const;

    /*@ManDoc: Returns a constant reference to the BoxArray that defines the
               valid region associated with this FabArray.
    */
    inline const BoxArray &boxArray () const;

    /*@ManDoc: Returns a constant reference to the Kth Box in the BoxArray.
               That is, the valid region of the Kth grid.
    */
    virtual inline const Box& box (int K) const;


    /*@ManDoc: Returns a constant reference to the Kth Fabs Box in the FabArray.
               That is, the region the Kth fab is actually defined on.
    */
    virtual inline const Box fabbox (int K) const;


    //
    //@ManDoc: Returns the length of this FabArray; i.e. the number of FABs.
    //
    inline int length () const;

    /*@ManDoc: Returns a constant reference to the FAB associated with the
               Kth element.
    */
//protected:
    inline const FAB& operator[] (int K) const;
    //
    //@ManDoc: Returns a reference to the FAB associated with the Kth element.
    //
    inline FAB& operator[] (int K);

public:
    //@ManDoc: Returns number of variables associated with each point (nvar).
    //
    inline int nComp () const;
    //
    //@ManDoc: Explicitly set the Kth FAB in the FabArray to point to elem.
    //
//private:
    void setFab (int  K,
		 FAB* elem);

public:

    /*@ManDoc: Removes the Kth FAB from the FabArray and returns it to the
               caller.  It leaves the Kth element undefined.  It is the
               callers responsibility to delete the pointer when no
               longer needed.
    */
    inline FAB* remove (int K);

    /*@ManDoc: Delete the Kth element in the FabArray only if the array is
	       managed.  It leaves the Kth element in the FabArray undefined
	       for both managed and unmanaged FabArrays.
    */
    inline void clear (int K);
    inline void clear ();

    /*@ManDoc: Delete all elements in the FabArray.  Use this only
	       on unmanaged FabArrays.  It leaves all elements
	       in the FabArray undefined.
    */
    inline void clearUnmanaged ();


    /*@ManDoc: Set all components in the valid region of each FAB in the
               FabArray to val.
    */
    void setVal (T val);

    /*@ManDoc: Set the value of num\_comp components in the valid region of
               each FAB in the FabArray, starting at component comp to val.
               Also set the value of nghost boundary cells.
    */
    void setVal (T   val,
                 int comp,
                 int num_comp,
                 int nghost = 0);

    /*@ManDoc: Set the value of num\_comp components in the valid region of
               each FAB in the FabArray, starting at component comp, as well
               as nghost boundary cells, to val, provided they also intersect
               with the Box region.
    */
    void setVal (T          val,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);

    /*@ManDoc: Set all components in the valid region of each FAB in the
               FabArray to val, including nghost boundary cells.
    */
    inline void setVal (T   val,
                        int nghost);

    /*@ManDoc: Set all components in the valid region of each FAB in the
               FabArray to val, including nghost boundary cells, that also
               intersect the Box region.
    */
    inline void setVal (T          val,
                        const Box& region,
                        int        nghost);
    //
    //@ManDoc: Set all values in the boundary region to val.
    //
    void setBndry (T val);

    /*@ManDoc: Set ncomp values in the boundary region, starting at
               start\_comp to val.
    */
    void setBndry (T   val,
                   int strt_comp,
                   int ncomp);

    /*@ManDoc: Copy on intersection within FabArray.  Data is copied from
               valid regions to intersecting regions of definition.  The
               purpose is to fill in the boundary regions of each FAB in
               the FabArray.
    */
    inline void FillBoundary();

    /*@ManDoc: Same as FillBoundary(), but only copies
               num\_comp components starting at start\_comp.
    */
    inline void FillBoundary(int start_comp, int num_comp);


    /*@ManDoc: This function copies data from fa to this FabArray.  Each FAB
               in fa is intersected with all FABs in this FabArray and a copy
               is performed on the region of intersection.  The intersection
               is restricted to the valid region of each FAB.
    */
    void copy (const FabArray<T,FAB>& fa);

    /*@ManDoc: This function copies data from src to this FabArray.  Each FAB
               in src is intersected with all FABs in this FabArray and a copy
               is performed on the region of intersection.  The intersection
               is restricted to the num\_comp components starting at src\_comp
               in the FabArray src, with the destination components in this
               FabArray starting at dest\_comp.  nghost boundary zones are
               included in the intersection as well.  This assumes that
               the source and destination FabArray have identical valid
               regions.

    */
    void copy (const FabArray<T,FAB>& src,
               int                    src_comp,
               int                    dest_comp,
               int                    num_comp,
               int                    nghost=0);

    /*@ManDoc: Copies the values contained in the intersection of the
               valid region of this FabArray with the FAB dest into dest.
    */
    void copy (FAB& dest) const;

    /*@ManDoc: Copies the values contained in the intersection of the
               valid region of this FabArray with the FAB dest and the Box
               subbox into that subregion of dest.
    */
    void copy (FAB&       dest,
               const Box& subbox) const;

    /*@ManDoc: Copies the values contained in the intersection of the
               num\_comp component valid region of this FabArray, starting at
               component src\_comp, with the FAB dest into dest, starting at
               component dest\_comp in dest.
    */
    void copy (FAB& dest,
               int  src_comp,
               int  dest_comp,
               int  num_comp) const;

    /*@ManDoc: Copies the values contained in the intersection of the
               num\_comp component valid region of this FabArray, starting at
               component src\_comp, with the FAB dest and the Box subbox, into
               dest, starting at component dest\_comp in dest.
    */
    void copy (FAB&       dest,
               const Box& subbox,
               int        src_comp,
               int        dest_comp,
               int        num_comp) const;


    const DistributionMapping &DistributionMap() const { return distributionMap; }

    aString FAName() { return faName; }
    void    SetFAName(const aString &faname) { faName = faname; }

protected:
    int         n_grow;
    int         n_comp;
    BoxArray    boxarray;
    PArray<FAB> fabparray;
    DistributionMapping distributionMap;
    aString faName;

private:
    //
    // These are disallowed.
    //
    FabArray (const FabArray<T,FAB>&);
    FabArray<T,FAB>& operator= (const FabArray<T,FAB>&);

    // This is used locally in all define functions
    void AllocFabs(int nvar);
};


// ========================================================================

template <class T, class FAB> class ConstFabArrayIterator;

// ------------------------------------------------------------------------
// FabArrayIterator
// ------------------------------------------------------------------------
template <class T, class FAB>
class FabArrayIterator {

  public:
    FabArrayIterator(FabArray<T, FAB> &fabarray, const int nghost = 0);
    virtual ~FabArrayIterator();

    inline const FAB &operator()() const;
    inline FAB &operator()();
    virtual inline const Box &validbox() const;
    virtual inline const Box &fabbox() const;
    virtual inline FabArrayIterator<T, FAB> &operator++();
    virtual inline bool isValid();
    virtual inline int index() const;

    inline int NGhost() const;

  protected:
    FabArray<T, FAB> &fabArray;
    int nGhost;
    int currentIndex;

  private:
    // These are disallowed.
    FabArrayIterator();
    FabArrayIterator(const FabArrayIterator<T,FAB> &);
    FabArrayIterator<T,FAB> &operator=(const FabArrayIterator<T,FAB> &);
};


// ------------------------------------------------------------------------
// DependentFabArrayIterator
// ------------------------------------------------------------------------
template <class T, class FAB>
class DependentFabArrayIterator {

  public:
    DependentFabArrayIterator(FabArrayIterator<T, FAB> &controllerfaiter,
		     FabArray<T, FAB> &dependentfabarray);
    DependentFabArrayIterator(FabArrayIterator<T, FAB> &controllerfaiter,
		     const FabArray<T, FAB> &dependentfabarray);
    virtual ~DependentFabArrayIterator();

    virtual inline const FAB &operator()() const;
    virtual inline FAB &operator()();
    virtual inline const Box &validbox() const;
    virtual inline const Box &fabbox() const;
    virtual inline int index() const;

  protected:
    FabArray<T, FAB> &fabArray;
    FabArrayIterator<T, FAB> *controller;

  private:
    // These are disallowed.
    DependentFabArrayIterator();
    DependentFabArrayIterator(const DependentFabArrayIterator<T,FAB> &);
    DependentFabArrayIterator<T,FAB>
		      &operator=(const DependentFabArrayIterator<T,FAB> &);
};


// ------------------------------------------------------------------------
// ConstFabArrayIterator
// ------------------------------------------------------------------------
template <class T, class FAB>
class ConstFabArrayIterator {

  public:
    ConstFabArrayIterator(const FabArray<T, FAB> &fabarray, const int nghost = 0);
    virtual ~ConstFabArrayIterator();

    virtual inline const FAB &operator()() const;
    virtual inline const Box &validbox() const;
    virtual inline const Box &fabbox() const;
    virtual inline ConstFabArrayIterator<T, FAB> &operator++();
    virtual inline int index() const;
    virtual inline bool isValid() const;
    inline int NGhost() const;

  protected:
    const FabArray<T, FAB> &fabArray;
    int nGhost;
    int currentIndex;

  private:
    // These are disallowed.
    ConstFabArrayIterator();
    ConstFabArrayIterator(const ConstFabArrayIterator<T,FAB> &);
    ConstFabArrayIterator<T,FAB> &operator=(const ConstFabArrayIterator<T,FAB> &);
};


// ------------------------------------------------------------------------
// ConstDependentFabArrayIterator
// ------------------------------------------------------------------------
template <class T, class FAB>
class ConstDependentFabArrayIterator {

  public:
    ConstDependentFabArrayIterator(
		     ConstFabArrayIterator<T, FAB> &controllerfaiter,
		     const FabArray<T, FAB> &dependentfabarray);
    virtual ~ConstDependentFabArrayIterator();

    virtual inline const FAB &operator()() const;
    virtual inline const Box &validbox() const;
    virtual inline const Box &fabbox() const;
    virtual inline int index() const;

  protected:

    const FabArray<T, FAB> &fabArray;
    ConstFabArrayIterator<T, FAB> *controller;

  private:
    // These are disallowed.
    ConstDependentFabArrayIterator();
    ConstDependentFabArrayIterator(const ConstDependentFabArrayIterator<T,FAB> &);
    ConstDependentFabArrayIterator<T,FAB>
		      &operator=(const ConstDependentFabArrayIterator<T,FAB> &);
};

// ========================================================================

class FabArrayId {
  public:
    FabArrayId() : fabArrayId(-1) { }
    FabArrayId(int newid) : fabArrayId(newid) { }
    ~FabArrayId() { }
    int Id() const    { return fabArrayId; }
  protected:
    void SetId(int id) { fabArrayId = id;   }
  private:
    int fabArrayId;
};


class FillBoxId {
  public:
    FillBoxId() : fillBoxId(-1) { }
    FillBoxId(int newid, const Box &fillbox, const BoxList &filledboxes) :
	     fillBoxId(newid), fillBox(fillbox), filledBoxes(filledboxes) { }
    ~FillBoxId() { }
    int Id() const                       { return fillBoxId;   }
    const Box &box() const               { return fillBox;     }
    const BoxArray &FilledBoxes() const  { return filledBoxes; }

  protected:
    void SetId(int id) { fillBoxId = id;   }
  private:
    int fillBoxId;
    Box fillBox;
    BoxArray filledBoxes;
};


// this class orchestrates filling a destination fab of size destFabBox
// from fabarray on the local processor (myProc)

template <class T, class FAB>
class FabArrayCopyDescriptor {

  public:

    FabArrayCopyDescriptor(bool cacheremotedata = true);
    virtual ~FabArrayCopyDescriptor();

    FabArrayId RegisterFabArray(FabArray<T, FAB> *fabarray);
    FillBoxId AddBox(const FabArrayId &fabarrayid, const Box &destFabBox,
		     BoxList &unfilledBoxes);
    FillBoxId AddBox(const FabArrayId &fabarrayid, const Box &destFabBox,
		     BoxList &unfilledBoxes,
		     int srccomp, int destcomp, int numcomp);
    // add a box but only from FabArray[fabarrayindex]
    FillBoxId AddBox(const FabArrayId &fabarrayid, const Box &destFabBox,
		     BoxList &unfilledBoxes,
		     int fabarrayindex,
		     int srccomp, int destcomp, int numcomp);
    //Array<FillBoxId> AddFabArrayBoxes(const FabArrayId &fabarrayid,
		     //int srccomp, int destcomp, int numcomp);
    void CollectData();
    void FillFab(const FabArrayId &fabarrayid, const FillBoxId &fillboxid,
		 FAB &destFab);

    void PrintStats();
    //const BoxList &UnfilledBoxes()  { return unfilledBoxes;    }
    //const Box &box() const          { return destBox;          }
    bool DataAvailable() const      { return dataAvailable;    }
    //bool CompletelyFilled()         { return completelyFilled; }

    enum FillType { FillLocally, FillRemotely, Unfillable };

  protected:

    struct FabCopyDescriptor {
        FabCopyDescriptor() : myProc(-1), copyFromProc(-1), copyFromIndex(-1),
			      fabArrayId(-1), fillBoxId(-1),
			      fillType(Unfillable),
			      cacheDataAllocated(false)
        {
        }

        ~FabCopyDescriptor() {
	  if(cacheDataAllocated) {
	    delete localFabSource;
	  }
        }

        int fabArrayId;
        int fillBoxId;
        int copyFromIndex;
        int myProc;
        int copyFromProc;
        int srcComp;
        int destComp;
        int nComp;
        Box subBox;
        FillType fillType;
        bool cacheDataAllocated;
        FAB *localFabSource;  // either a local fab from this fabarray or a
			      // temp local fab with data from a remote copy
    };

    Array< FabArray<T, FAB> * > fabArrays;
    Array<BoxList>              unfilledBoxes;     // for each fabArray
    Array<bool>                 completelyFilled;  // for each fabArray
    List<FabCopyDescriptor *>   fabCopyDescList;
    List<FabComTag>             fabComTagList;
    int nextFabArrayId;
    int nextFillBoxId;
    bool dataAvailable;
    bool cacheRemoteData;

    // stats
    Array<int>  nRemoteBoxes;  // per fabArray[index]
    Array<long> nRemoteBytes;  // per fabArray[index]
    int  totalRemoteBoxes;
    long totalRemoteBytes;     // total bytes to communicate to fill destFab
    //int nSubBoxes;             // how many pieces is destFabBox broken into

  private:
    // These are disallowed.
    FabArrayCopyDescriptor();
    FabArrayCopyDescriptor(const FabArrayCopyDescriptor<T,FAB> &);
    FabArrayCopyDescriptor<T,FAB> &operator=(const FabArrayCopyDescriptor<T,FAB> &);
};

// ========================================================================



#ifdef BL_USE_LPARX_CONSTRUCTS
#undef BL_USE_LPARX_CONSTRUCTS
/*@ManDoc:
  The macro BL\_USE\_LPARX\_CONSTRUCTS enables the use of the "forallthis",
  "forall", "pforallthis" and "pforall" macros in some earlier versions
  of BoxLib.  The BL\_USE\_LPARX\_CONSTRUCTS macro is not turned on by default.
  If you use any of the four above macros in your code, you should probably
  replace them with the equivalent macro expansion, as the
  BL\_USE\_LPARX\_CONSTRUCTS macro will probably go away in the future.
*/
#define BL_USE_LPARX_CONSTRUCTS 1
#endif /*BL_USE_LPARX_CONSTRUCTS*/

#ifdef BL_USE_LPARX_CONSTRUCTS
#define forallthis(i) {                             \
    for (int i = 0; i < length(); i++) {

#define forall(i, fabarray) {                       \
    for (int i = 0; i < (fabarray).length(); i++) {
#define endforall }}

#define pforallthis(i) {                            \
    for (int i = 0; i < length(); i++) {

#define pforall(i, fabarray) {                      \
    for (int i = 0; i < (fabarray).length(); i++) {
#define pendforallnosync }}
#define pendforall pendforallnosync
#endif



//
// Inlines.
//
template <class T, class FAB>
inline
int
FabArray<T,FAB>::nGrow () const
{
    return n_grow;
}

template <class T, class FAB>
inline
const BoxArray&
FabArray<T,FAB>::boxArray () const
{
    return boxarray;
}


template <class T, class FAB>
inline
const Box&
FabArray<T,FAB>::box (int K) const
{
    return boxarray[K];
}

template <class T, class FAB>
inline
const Box
FabArray<T,FAB>::fabbox (int K) const
{
    // dont use fabparray[K] because it may not be valid in parallel
    return (::grow(boxarray[K], n_grow));
}

template <class T, class FAB>
inline
int
FabArray<T,FAB>::length () const
{
    return boxarray.length();
}

template <class T, class FAB>
inline
const FAB&
FabArray<T,FAB>::operator[] (int K) const
{
    return fabparray[K];
}

template <class T, class FAB>
inline
FAB&
FabArray<T,FAB>::operator[] (int K)
{
    return fabparray[K];
}

template <class T, class FAB>
inline
int
FabArray<T,FAB>::nComp () const
{
    return n_comp;
}

template <class T, class FAB>
inline
FAB*
FabArray<T,FAB>::remove (int K)
{
    return fabparray.remove(K);
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::clear (int K)
{
    fabparray.clear(K);
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::clear ()
{
    fabparray.clear();
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::clearUnmanaged ()
{
  for(FabArrayIterator<T, FAB> fai(*this); fai.isValid(); ++fai) {
    FAB *fab = fabparray.remove(fai.index());
    delete fab;
  }
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::setVal (T   val,
                         int nghost)
{
    setVal(val,0,n_comp,nghost);
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::setVal (T          val,
                         const Box& region,
                         int        nghost)
{
    setVal(val,region,0,n_comp,nghost);
}

#ifdef BL_AUTO_INSTANTIATE
#include <FabArray.C>
#endif

#endif
