//BL_COPYRIGHT_NOTICE

//
// $Id: Array.C,v 1.1 1997-07-08 23:08:04 vince Exp $
//

#include <Assert.H>
#include <Misc.H>

template <class T>
Array<T>::Array (long     len,
                 const T& initialValue)
{
    assert(len >= 0);
    nelem = len;
    vp    = new T[len];
    if (vp == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
    truesize = nelem;
    for(long i = 0; i < nelem; ++i)
        vp[i] = initialValue;
}

template <class T>
Array<T>::Array (const T* vec,
                 long     len)
{
    assert(len >= 0);
    nelem = len;
    vp    = new T[len];
    if (vp == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
    truesize = nelem;
    for(long i = 0; i < nelem; ++i)
        vp[i] = vec[i];
}

template <class T>
Array<T>::Array (const Array<T>& a)
{
    nelem = a.nelem;
    vp    = new T[nelem];
    if (vp == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
    truesize = nelem;
    for (long i = 0; i < nelem; i++)
        vp[i] = a.vp[i];
}

template <class T>
Array<T>&
Array<T>::operator= (const Array<T>& sa)
{
    if (this != &sa)
    {
        clear();
        vp       = new T[sa.nelem];
        nelem    = sa.nelem;
        truesize = nelem;
        if (vp == 0)
            BoxLib::OutOfMemory(__FILE__, __LINE__);
        for(long i = 0; i < nelem; i++)
            vp[i] = sa.vp[i];
    }
    return *this;
}

template <class T>
void
Array<T>::resize (long newlen)
{
    if (newlen == nelem)
        return;
    if (newlen <= truesize)
    {
        nelem = newlen;
        return;
    }
    T* newvp = new T[newlen];
    long len = Min(newlen,nelem);
    if (newvp == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
    for (long i = 0; i < len; i++)
        newvp[i] = vp[i];
    delete [] vp;
    vp = newvp;
    nelem = newlen;
    truesize = newlen;
}

template <class T>
void Array<T>::resize (long     newlen,
                       const T& initialValue)
{
    if (newlen == nelem)
        return;
    if (newlen <= truesize)
    {
        for(long i = nelem; i < newlen; ++i)
            vp[i] = initialValue;
        nelem = newlen;
        return;
    }
    T* newvp = new T[newlen];
    if (newvp == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
    long len = Min(newlen,nelem);
    long i;
    for (i = 0; i < len; i++)
        newvp[i] = vp[i];
    for(i = len; i < newlen; ++i)
        newvp[i] = initialValue;
    delete [] vp;
    vp = newvp;
    nelem = newlen;
    truesize = newlen;
}

template <class T>
void
Array<T>::reserve (long _truesize)
{
    if (_truesize > truesize)
    {
        T* newvp = new T[_truesize];
        if (newvp == 0)
            BoxLib::OutOfMemory(__FILE__, __LINE__);
        for (long i = 0; i < nelem; i++)
            newvp[i] = vp[i];
        delete [] vp;
        vp = newvp;
        truesize = _truesize;
    }
}

template <class T>
void
Array<T>::shrinkWrap ()
{
    if (nelem != truesize)
    {
        T* newvp = new T[nelem];
        if (newvp == 0)
            BoxLib::OutOfMemory(__FILE__, __LINE__);
        for (long i = 0; i < nelem; i++)
            newvp[i] = vp[i];
        delete [] vp;
        vp = newvp;
        truesize = nelem;
    }
}

template <class T>
bool
Array<T>::operator== (const Array<T>& rhs) const
{
    bool rc = true;
    if (length() != rhs.length())
        rc = false;
    else
    {
        for (long i = 0; i < length() && rc; ++i)
            if (!((*this)[i] == rhs[i]))
                rc = false;
    }
    return rc;
}
