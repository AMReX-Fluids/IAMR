//BL_COPYRIGHT_NOTICE

//
// $Id: Tuple.C,v 1.1 1997-07-08 23:08:07 vince Exp $
//

template <class T, size_t N>
Tuple<T,N>::Tuple (const T* v)
{
    assert(v != 0);
    for (size_t i = 0; i < N; ++i)
        vect[i] = v[i];
}

template <class T, size_t N>
Tuple<T,N>::Tuple (const Tuple<T,N>& rhs)
{
    for (size_t i = 0; i < N; ++i)
        vect[i] = rhs.vect[i];
}

template <class T, size_t N>
Tuple<T,N>&
Tuple<T,N>::operator= (const Tuple<T,N>& rhs)
{
    for (size_t i = 0; i < N; ++i)
        vect[i] = rhs.vect[i];
    return *this;
}

