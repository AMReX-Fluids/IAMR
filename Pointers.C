//BL_COPYRIGHT_NOTICE

//
// $Id: Pointers.C,v 1.1 1997-07-08 23:08:07 vince Exp $
//

#if 0
//
// This stuff is commented out for now.  If someone really uses it,
// we'll uncomment it as well as document it.
//

template<class T>
VcClassPtr<T>::VcClassPtr (const VcClassPtr<T>& _a)
    :  ptr(_a.isNull() ?  0 : _a.ptr->clone())
{}

template<class T>
VcClassPtr<T>&
VcClassPtr<T>::operator= (T* _ptr)
{
    delete ptr;
    ptr = _ptr;
    return *this;
}

template<class T>
VcClassPtr<T>&
VcClassPtr<T>::operator= (const VcClassPtr<T>& _r)
{
    if (ptr != _r.ptr)
    {
        delete ptr;
        ptr = _r.isNull() ? 0 : _r.ptr->clone();
    }
    return *this;
}

template<class T>
T*
VcClassPtr<T>::release()
{
    T* old = ptr;
    ptr = 0;
    return old;
}
#endif
