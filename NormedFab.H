//BL_COPYRIGHT_NOTICE

#ifndef NORMEDFAB_H
#define NORMEDFAB_H

//
// $Id: NormedFab.H,v 1.1 1997-07-08 23:08:07 vince Exp $
//
#include <iostream.h>

#include <BaseFab.H>
#include <REAL.H>

//
//@Man:
//@Memo: A BaseFab with Norms
/*@Doc:

  This class provides various L-p norms to a BaseFab.

  This class does NOT provide a copy constructor or assignment operator.
 */

template <class T>
class NormedFab
    : public virtual BaseFab<T>
{
public:
    //
    //@ManDoc: The default constructor.
    //
    NormedFab ();
    //
    //@ManDoc: Construct a NormedFab from a Box.
    //
    NormedFab (const Box& b,
               int        n = 1);
    //
    //@ManDoc: Construct a NormedFab from a BaseFab.
    //
    NormedFab (BaseFab<T>& bf,
               Box         subb,
               int         ns,
               int         nc);
    //
    //@ManDoc: The (virtual) destructor.
    //
    virtual ~NormedFab();

    /*@ManDoc: Compute the Lp-norm of this FAB using components
               (comp : comp+numcomp-1).
               p < 0  -> ERROR
               p = 0  -> infinity norm (max norm)
               p = 1  -> sum of ABS(FAB)
    */
    virtual Real norm (int p,
                       int comp=0,
                       int numcomp=1) const;
    //
    //@ManDoc: Same as above except only on given subbox.
    //
    virtual Real norm (const Box& subbox,
                       int        p,
                       int        comp=0,
                       int        numcomp=1) const;
    //
    //@ManDoc:  Compute absolute value for all components of this FAB.
    //
    void abs ();
    //
    //@ManDoc: Same as above except only for components (comp: comp+numcomp-1)
    //
    void abs (int comp,
              int numcomp=1);
    //
    //@ManDoc: Calculate abs() on subbox for given component range.
    //
    void abs (const Box& subbox,
              int        comp = 0,
              int        numcomp=1);
private:
    //
    // These are not allowed.
    //
    NormedFab (const NormedFab<T>&);
    NormedFab<T>& operator= (const NormedFab<T>&);
};

#ifdef BL_AUTO_INSTANTIATE
#include <NormedFab.C>
#endif

#endif /*NORMEDFAB_H*/
