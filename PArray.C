//BL_COPYRIGHT_NOTICE

//
// $Id: PArray.C,v 1.1 1997-07-08 23:08:07 vince Exp $
//

#include <iostream.h>

template <class T>
PArray<T>::PArray (int    len,
                   PArrayPolicy _managed)
    : vp(len),
      managed(_managed)
{
    for (int i = 0; i < length(); i++)
        vp[i] = 0;
}

template <class T>
void
PArray<T>::clear ()
{
    if (managed)
    {
        for (int i = 0; i < length(); ++i) {
	  if(vp[i] == 0) {
	    //cerr << "Error in PArray<T>::clear vp[" << i << "] = 0 (fix parallel)"
		 //<< endl;
	  } else {
            delete ((T*)(vp[i]));
	  }
	}
    }
    for (int i = 0; i < length(); ++i)
        vp[i] = 0;
}


template <class T>
PArray<T>::~PArray ()
{
    clear();
}

template <class T>
void
PArray<T>::resize (int newlen)
{
    void **ovp = vp.ready() ? &vp[0] : 0;
    int onelem = vp.length();
    vp.resize(newlen);
    for (int i = onelem; i < newlen; ++i)
        vp[i] = 0;
    if (managed)
    {
        for (int i = newlen; i < onelem; i++)
        {
	  if(ovp[i] == 0) {
	    //cerr << "Error in PArray<T>::clear ovp[" << i << "] = 0 (fix parallel)"
		 //<< endl;
	  } else {
            delete ((T*)(ovp[i]));
            ovp[i]  = 0;
	  }
        }
    }
}

template <class T>
void
PArray<T>::clear (int n)
{
    if (managed) {
	  if(vp[n] == 0) {
	    //cerr << "Error in PArray<T>::clear vp[" << n << "] = 0 (fix parallel)"
		 //<< endl;
	  } else {
            delete ((T*)(vp[n]));
	  }
    }
    vp[n] = 0;
}
