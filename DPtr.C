//BL_COPYRIGHT_NOTICE

//
// $Id: DPtr.C,v 1.1 1997-07-08 23:08:06 vince Exp $
//

template <class T>
DPtrRep<T>::~DPtrRep ()
{}

template <class T>
T*
DPtrRep<T>::dataPtr ()
{
    return &(*this)[0];
}

template <class T>
void
DPtrRep<T>::Lock ()
{}

template <class T>
void
DPtrRep<T>::UnLock ()
{}

template <class T>
void
DPtrRep<T>::Demand ()
{}

template <class T>
void
DPtrRep<T>::Queue ()
{}
