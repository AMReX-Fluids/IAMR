//BL_COPYRIGHT_NOTICE

//
// $Id: ArithFab.C,v 1.1 1997-07-08 23:08:04 vince Exp $
//

#include <iostream.h>
#include <iomanip.h>

#include <Misc.H>
#include <ArithFab.H>
#include <Looping.H>

template <class T>
ArithFab<T>::ArithFab ()
    : BaseFab<T>()
{}

template <class T>
ArithFab<T>::ArithFab (const Box &b,
                       int        n)
    : BaseFab<T>(b,n)
{}

template <class T>
ArithFab<T>::ArithFab (BaseFab<T> &bf,
                       Box         subb,
                       int         ns,
                       int         nc)
    : BaseFab<T>(bf,subb,ns,nc)
{}

//
// This is not inlined as it's virtual.
//

template <class T>
ArithFab<T>::~ArithFab ()
{}

template <class T>
void
ArithFab<T>::patternFill (int mark)
{
    ForAllThis(T)
    {
        thisR = D_TERM(iR*100, +jR*10, + kR) + 1000*nR + 10000*mark;
    } EndFor
}

// ----------------------  COPYREV MEMBERS ------------------------

template <class T>
void
ArithFab<T>::copyRev (const Box&         destbox,
                      const ArithFab<T>& src,
                      const Box&         srcbox,
                      int                reversal_index,
                      T*                 multiplier)
{
    ArithFab<T> &dest = *this;
    int sc = 0;
    int dc = 0;
    int nv = nComp();
    ForAllRevXBNYCBNNN(T,dest,destbox,dc,src,srcbox,sc,nv,reversal_index)
    {
        destR = multiplier[_n]*srcR;
    } EndFor
}

template <class T>
void
copyFABRev (const ArithFab<T>& src,
            const Box&         srcbox,
            ArithFab<T>&       dest,
            const Box&         destbox,
            int                reversal_index,
            int                comp,
            T                  multiplier)
{
    int sc = comp;
    int dc = comp;
    int nv = 1;

    ForAllRevXBNYCBNNN(T,dest,destbox,dc,src,srcbox,sc,nv,reversal_index)
    {
        destR = multiplier * srcR;
    } EndFor
}

//
// ---------------------  Sum
//
template <class T>
T
ArithFab<T>::sum (int comp,
                  int numcomp) const
{
    T *_sum_row = 0;
    int _sum_len = 0;
    ForAllThisCPencil(T,domain,comp,numcomp)
    {
        const T* _row = &thisR;
        if (_sum_row == 0)
        {
            _sum_row = new T[thisLen];
            if (_sum_row == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
            _sum_len = thisLen;
            for (int i = 0; i < thisLen; i++)
                _sum_row[i] = _row[i];
        }
        else
        {
            for (int i = 0; i < thisLen; i++)
                _sum_row[i] += _row[i];
        }
    } EndForPencil;

    T _sum = _sum_row[0];
    for (int i = 1; i < _sum_len; i++)
        _sum += _sum_row[i];
    delete [] _sum_row;

    return _sum;
}

template <class T>
T
ArithFab<T>::sum (const Box& subbox,
                  int        comp,
                  int        numcomp) const
{
    T *_sum_row = 0;
    int _sum_len = 0;
    ForAllThisCPencil(T,subbox,comp,numcomp)
    {
        const T* _row = &thisR;
        if (_sum_row == 0)
        {
            _sum_row = new T[thisLen];
            if (_sum_row == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
            _sum_len = thisLen;
            for (int i = 0; i < thisLen; i++)
                _sum_row[i] = _row[i];
        }
        else
        {
            for (int i = 0; i < thisLen; i++)
            {
                _sum_row[i] += _row[i];
            }
        }
    } EndForPencil;

    T _sum = _sum_row[0];
    for (int i = 1; i < _sum_len; i++)
        _sum += _sum_row[i];
    delete [] _sum_row;

    return _sum;
}

//
// ---------------------  Negate
//

template <class T>
ArithFab<T>&
ArithFab<T>::negate ()
{
    ForAllThis(T)
    {
        thisR = - thisR;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::negate (int comp,
                     int numcomp)
{
    ForAllThisNN(T,comp,numcomp)
    {
        thisR = - thisR;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::negate (const Box& b,
                     int        comp,
                     int        numcomp)
{
    ForAllThisBNN(T,b,comp,numcomp)
    {
        thisR = - thisR;
    } EndFor
    return *this;
}

//
// ---------------------  Invert
//

template <class T>
ArithFab<T>&
ArithFab<T>::invert (T r)
{
    ForAllThis(T)
    {
        thisR = r/thisR;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::invert (T   r,
                     int comp,
                     int numcomp)
{
    ForAllThisNN(T,comp,numcomp)
    {
        thisR = r/thisR;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::invert (T          r,
                     const Box& b,
                     int        comp,
                     int        numcomp)
{
    ForAllThisBNN(T,b,comp,numcomp)
    {
        thisR = r/thisR;
    } EndFor
    return *this;
}

//
// ---------------------  Scalar Addition
//

template <class T>
ArithFab<T>&
ArithFab<T>::operator += (T r)
{
    ForAllThis(T)
    {
        thisR += r;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::plus (T   r,
                   int comp,
                   int numcomp)
{
    ForAllThisNN(T,comp,numcomp)
    {
        thisR += r;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::plus (T          r,
                   const Box& b,
                   int        comp,
                   int        numcomp)
{
    ForAllThisBNN(T,b,comp,numcomp)
    {
        thisR += r;
    } EndFor
    return *this;
}

//
// ---------------------  Fab Addition
//

template <class T>
ArithFab<T>&
ArithFab<T>::operator += (const ArithFab<T>& x)
{
    ForAllThisXC(T,x)
    {
        thisR += xR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::plus (const ArithFab<T>& src,
                   int                srccomp,
                   int                destcomp,
                   int                numcomp)
{
    ForAllThisBNNXC(T,domain,destcomp,numcomp,src,srccomp)
    {
        thisR += srcR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::plus (const ArithFab<T>& src,
                   const Box&         subbox,
                   int                srccomp,
                   int                destcomp,
                   int                numcomp)
{
    ForAllThisBNNXC(T,subbox,destcomp,numcomp,src,srccomp)
    {
        thisR += srcR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::plus (const ArithFab<T>& src,
                   const Box&         srcbox,
                   const Box&         destbox,
                   int                srccomp,
                   int                destcomp,
                   int                numcomp)
{
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
        thisR += srcR;
    } EndForTX
    return *this;
}

//
// ---------------------  Fab Subtraction
//

template <class T>
ArithFab<T>&
ArithFab<T>::operator -= (const ArithFab<T>& x)
{
    ForAllThisXC(T,x)
    {
        thisR -= xR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::minus (const ArithFab<T>& src,
                    int                srccomp,
                    int                destcomp,
                    int                numcomp)
{
    ForAllThisBNNXC(T,domain,destcomp,numcomp,src,srccomp)
    {
        thisR -= srcR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::minus (const ArithFab<T>& src,
                    const Box&         subbox,
                    int                srccomp,
                    int                destcomp,
                    int                numcomp)
{
    ForAllThisBNNXC(T,subbox,destcomp,numcomp,src,srccomp)
    {
        thisR -= srcR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::minus (const ArithFab<T>& src,
                    const Box&         srcbox,
                    const Box&         destbox,
                    int                srccomp,
                    int                destcomp,
                    int                numcomp)
{
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
        thisR -= srcR;
    } EndForTX
    return *this;
}

//
// ---------------------  Scalar Multiplication
//

template <class T>
ArithFab<T>&
ArithFab<T>::operator *= (T r)
{
    ForAllThis(T)
    {
        thisR *= r;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::mult (T   r,
                   int comp,
                   int numcomp)
{
    ForAllThisNN(T,comp,numcomp)
    {
        thisR *= r;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::mult (T          r,
                   const Box& b,
                   int        comp,
                   int        numcomp)
{
    ForAllThisBNN(T,b,comp,numcomp)
    {
        thisR *= r;
    } EndFor
    return *this;
}

//
// ---------------------  Fab Multiplication
//

template <class T>
ArithFab<T>&
ArithFab<T>::operator *= (const ArithFab<T> &x)
{
    ForAllThisXC(T,x)
    {
        thisR *= xR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::mult (const ArithFab<T>& src,
                   int                srccomp,
                   int                destcomp,
                   int                numcomp)
{
    ForAllThisBNNXC(T,domain,destcomp,numcomp,src,srccomp)
    {
        thisR *= srcR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::mult (const ArithFab<T>& src,
                   const Box&         subbox,
                   int                srccomp,
                   int                destcomp,
                   int                numcomp)
{
    ForAllThisBNNXC(T,subbox,destcomp,numcomp,src,srccomp)
    {
        thisR *= srcR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::mult (const ArithFab<T>& src,
                   const Box&         srcbox,
                   const Box&         destbox,
                   int                srccomp,
                   int                destcomp,
                   int                numcomp)
{
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
        thisR *= srcR;
    } EndForTX
    return *this;
}

//
// ---------------------  Scalar Division
//

template <class T>
ArithFab<T>&
ArithFab<T>::operator /= (T r)
{
    ForAllThis(T)
    {
        thisR /= r;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::divide (T   r,
                     int comp,
                     int numcomp)
{
    ForAllThisNN(T,comp,numcomp)
    {
        thisR /= r;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::divide (T          r,
                     const Box& b,
                     int        comp,
                     int        numcomp)
{
    ForAllThisBNN(T,b,comp,numcomp)
    {
        thisR /= r;
    } EndFor
    return *this;
}

//
// ---------------------  Fab Division
//

template <class T>
ArithFab<T>&
ArithFab<T>::operator /= (const ArithFab<T> &x)
{
    ForAllThisXC(T,x)
    {
        thisR /= xR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::divide (const ArithFab<T>& src,
                     int                srccomp,
                     int                destcomp,
                     int                numcomp)
{
    ForAllThisBNNXC(T,domain,destcomp,numcomp,src,srccomp)
    {
        thisR /= srcR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::divide (const ArithFab<T>& src,
                     const Box&         subbox,
                     int                srccomp,
                     int                destcomp,
                     int                numcomp)
{
    ForAllThisBNNXC(T,subbox,destcomp,numcomp,src,srccomp)
    {
        thisR /= srcR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::divide (const ArithFab<T>& src,
                     const Box&         srcbox,
                     const Box&         destbox,
                     int                srccomp,
                     int                destcomp,
                     int                numcomp)
{
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
        thisR /= srcR;
    } EndForTX
    return *this;
}

//
// Linear Interpolation / Extrapolation
// Result is (t2-t)/(t2-t1)*f1 + (t-t1)/(t2-t1)*f2
// Data is taken from b1 region of f1, b2 region of f2
// and stored in b region of this FAB.
// Boxes b, b1 and b2 must be the same size.
// Data is taken from component comp1 of f1, comp2 of f2,
// and stored in component comp of this FAB.
// This fab is returned as a reference for chaining.
//

template <class T>
ArithFab<T>&
ArithFab<T>::linInterp (const ArithFab<T>& f1,
                        const Box&         b1,
                        int                comp1,
                        const ArithFab<T>& f2,
                        const Box&         b2,
                        int                comp2,
                        Real               t1,
                        Real               t2,
                        Real               t,
                        const Box&         b,
                        int                comp,
                        int                numcomp)

{
    Real alpha = (t2-t)/(t2-t1);
    Real beta = (t-t1)/(t2-t1);
    ForAllThisBNNXCBNYCBN(T,b,comp,numcomp,f1,b1,comp1,f2,b2,comp2)
    {
        thisR = (T) (alpha*Real(f1R) + beta*Real(f2R));
    } EndForTX
    return *this;
}

//
// Linear combination, Result is alpha*f1 + beta*f2
// Data is taken from b1 region of f1, b2 region of f2
// and stored in b region of this FAB.
// Boxes b, b1 and b2 must be the same size.
// Data is taken from component comp1 of f1, comp2 of f2,
// and stored in component comp of this FAB.
// This fab is returned as a reference for chaining.
//

template <class T>
ArithFab<T>&
ArithFab<T>::linComb (const ArithFab<T>& f1,
                      const Box&         b1,
                      int                comp1,
                      const ArithFab<T>& f2,
                      const Box&         b2,
                      int                comp2,
                      Real               alpha,
                      Real               beta,
                      const Box&         b,
                      int                comp,
                      int                numcomp)

{
    ForAllThisBNNXCBNYCBN(T,b,comp,numcomp,f1,b1,comp1,f2,b2,comp2)
    {
        thisR = (T) (alpha*Real(f1R) + beta*Real(f2R));
    } EndForTX
    return *this;
}




