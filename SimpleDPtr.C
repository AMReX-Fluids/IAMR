//BL_COPYRIGHT_NOTICE

//
// $Id: SimpleDPtr.C,v 1.1 1997-07-08 23:08:07 vince Exp $
//

#include <BArena.H>

template <class T>
SimpleDPtr<T>::SimpleDPtr (size_t _size,
                           Arena* _arena)
    : DPtrRep<T>(_arena),
      dp(0)
{
    if (arena == 0)
    {
        static BArena _builtin_BArena;
        arena = &_builtin_BArena;
    }
    define(_size);
}

template <class T>
SimpleDPtr<T>::~SimpleDPtr ()
{
    clear();
}

template <class T>
void
SimpleDPtr<T>::clear ()
{
    if (arena)
        arena->free(dp);
    dp = 0;
    currentsize = 0;
}

template <class T>
void
SimpleDPtr<T>::resize (size_t _size)
{
    if (_size > currentsize)
    {
        T *g = (T*) arena->alloc(_size * sizeof(T));
        for (size_t i = 0; i < currentsize; ++i)
            g[i] = dp[i];
        delete [] dp;
        dp = g;
        currentsize = _size;
    }
}

template <class T>
void
SimpleDPtr<T>::define (size_t _size)
{
    assert(dp == 0);
    dp = (T*) arena->alloc(_size*sizeof(T));
    currentsize = _size;
    assert(dp != 0);
}

template <class T>
T&
SimpleDPtr<T>::operator[] (long n) const
{
    assert(n >= 0 && n < currentsize);
    return dp[n];
}

template <class T>
size_t
SimpleDPtr<T>::length ()
{
    return currentsize;
}

template <class T>
size_t
SimpleDPtr<T>::size () const
{
    return currentsize;
}
