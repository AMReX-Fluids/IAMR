//BL_COPYRIGHT_NOTICE

//
// $Id: AliasedDPtr.C,v 1.1 1997-07-08 23:08:02 vince Exp $
//

#include <AliasedDPtr.H>

template <class T>
AliasedDPtr<T>::AliasedDPtr (size_t _size)
{
    dptr        = new SimpleDPtr<T>(_size);
    currentsize = _size;
    offset      = 0;

    if (dptr.isNull())
        BoxLib::OutOfMemory(__FILE__, __LINE__);
}

template <class T>
AliasedDPtr<T>::AliasedDPtr (AliasedDPtr<T>* _asdptr,
                             size_t          _size,
                             int             _offset)
{
    dptr        = _asdptr->dptr;
    offset      = _offset;
    currentsize = _size;
}

template <class T>
AliasedDPtr<T>::~AliasedDPtr ()
{}

template <class T>
size_t
AliasedDPtr<T>::size () const
{
    return currentsize;
}

template <class T>
T&
AliasedDPtr<T>::operator[] (long n) const
{
    return (*dptr)[n + offset];
}

template <class T>
void
AliasedDPtr<T>::define (size_t _size)
{
    dptr        = new SimpleDPtr<T>(_size);
    currentsize = _size;
    offset      = 0;

    if (dptr.isNull())
        BoxLib::OutOfMemory(__FILE__, __LINE__);
}

template <class T>
void
AliasedDPtr<T>::resize (size_t _size)
{
    dptr        = new SimpleDPtr<T>(_size);
    currentsize = _size;
    offset      = 0;

    if (dptr.isNull())
        BoxLib::OutOfMemory(__FILE__, __LINE__);
}

template <class T>
void
AliasedDPtr<T>::clear ()
{
    dptr = 0;
}
