//BL_COPYRIGHT_NOTICE

//
// $Id: NormedFab.C,v 1.1 1997-07-08 23:08:07 vince Exp $
//

#include <math.h>
#include <BoxLib.H>
#include <Looping.H>

template <class T>
NormedFab<T>::NormedFab ()
    : BaseFab<T>()
{}

template <class T>
NormedFab<T>::NormedFab (const Box& b,
                         int        n)
    : BaseFab<T>(b,n)
{}

template <class T>
NormedFab<T>::NormedFab (BaseFab<T>& bf,
                         Box         subb,
                         int         ns,
                         int         nc)
    : BaseFab<T>(bf,subb,ns,nc)
{}

//
// This isn't inlined as it's virtual.
//

template <class T>
NormedFab<T>::~NormedFab()
{}

template <class T>
void
NormedFab<T>::abs ()
{
    ForAllThis(Real)
    {
        thisR = Abs(thisR);
    } EndFor
}

template <class T>
void
NormedFab<T>::abs (int comp,
                   int numcomp)
{
    ForAllThisNN(Real,comp,numcomp)
    {
        thisR = Abs(thisR);
    } EndFor
}

template <class T>
void
NormedFab<T>::abs (const Box& subbox,
                   int        comp,
                   int        numcomp)
{
    ForAllThisBNN(Real,subbox,comp,numcomp)
    {
        thisR = Abs(thisR);
    } EndFor
}

// ----------------------  NORM MEMBERS ------------------------

//
// This isn't inlined as it's virtual.
//

template <class T>
Real
NormedFab<T>::norm (int p,
                    int comp,
                    int numcomp) const
{
    return norm(domain,p,comp,numcomp);
}


template <class T>
Real
NormedFab<T>::norm (const Box& subbox,
                    int        p,
                    int        comp,
                    int        numcomp) const
{
    assert(comp >= 0 && comp+numcomp <= nComp());
    assert(p >= 0);

    Real* tmp  = 0;
    int tmplen = 0;
    Real nrm   = 0;
    if (p == 0)
    {
        ForAllThisCPencil(T,subbox,comp,numcomp)
        {
            const T* row = &thisR;
            if (tmp == 0)
            {
                tmp = new Real[thisLen];
                if (tmp == 0)
                    BoxLib::OutOfMemory(__FILE__, __LINE__);
                tmplen = thisLen;
                for (int i = 0; i < thisLen; i++) {
                    tmp[i] = Abs(Real(row[i]));
                    //tmp[i] = fabs(row[i]);
		}
            }
            else
            {
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = Max(tmp[i],Real(Abs(row[i])));
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++)
            nrm = Max(nrm, tmp[i]);
    }
    else if (p == 1)
    {
        ForAllThisCPencil(T,subbox,comp,numcomp)
        {
            const T* row = &thisR;
            if (tmp == 0)
            {
                tmp = new Real[thisLen];
                if (tmp == 0)
                    BoxLib::OutOfMemory(__FILE__, __LINE__);
                tmplen = thisLen;
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = Abs(Real(row[i]));
            }
            else
            {
                for (int i = 0; i < thisLen; i++)
                    tmp[i] += Abs(Real(row[i]));
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++)
            nrm += tmp[i];
    }
    else
      BoxLib::Error("NormedFab::norm(): only p == 0 or p == 1 are supported");

    delete [] tmp;

    return nrm;
}


#ifdef BL_ARCH_CRAY
#define RESTRICT restrict
// template specialization for Real
Real
NormedFab<Real>::norm (const Box& subbox,
		       int        p,
		       int        comp,
		       int        numcomp) const
{
    assert(comp >= 0 && comp+numcomp <= nComp());
    assert(p >= 0);

    Real* RESTRICT tmp  = 0;
    int tmplen = 0;
    Real nrm   = 0;
    if (p == 0)
    {
        ForAllThisCPencil(Real,subbox,comp,numcomp)
        {
            const Real* RESTRICT row = &thisR;
            if (tmp == 0)
            {
                tmp = new Real[thisLen];
                if (tmp == 0)
                    BoxLib::OutOfMemory(__FILE__, __LINE__);
                tmplen = thisLen;
#pragma _CRI ivdep
                for (int i = 0; i < thisLen; i++) {
                    tmp[i] = fabs(row[i]);
		}
            }
            else
            {
#pragma _CRI ivdep
	      for (int i = 0; i < thisLen; i++){
		Real a = fabs(row[i]);
		tmp[i] = tmp[i]>a ? tmp[i] : a ;
	      }
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++) {
	    Real a = tmp[i];
	    nrm = nrm>a ? nrm : a ;
	}
    }
    else if (p == 1)
    {
        ForAllThisCPencil(Real,subbox,comp,numcomp)
        {
            const Real* row = &thisR;
            if (tmp == 0)
            {
                tmp = new Real[thisLen];
                if (tmp == 0)
                    BoxLib::OutOfMemory(__FILE__, __LINE__);
                tmplen = thisLen;
#pragma _CRI ivdep
                for (int i = 0; i < thisLen; i++){
		    tmp[i] = fabs(row[i]);
		}
            }
            else
            {
#pragma _CRI ivdep
	      for (int i = 0; i < thisLen; i++) {
		    tmp[i] += fabs( row[i] );
	      }
            }
        } EndForPencil
        nrm = tmp[0];
#pragma _CRI ivdep
        for (int i = 1; i < tmplen; i++)
            nrm += tmp[i];
    }
    else
      BoxLib::Error("NormedFab<Real>::norm(): only p == 0 or p == 1 are supported");

    delete [] tmp;

    return nrm;
}
#endif
