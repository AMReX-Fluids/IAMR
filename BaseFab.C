//BL_COPYRIGHT_NOTICE

//
// $Id: BaseFab.C,v 1.1 1997-07-08 23:08:05 vince Exp $
//

#include <stdlib.h>
#include <Assert.H>
#include <Box.H>
#include <BoxList.H>
#include <Looping.H>
#include <AliasedDPtr.H>

#ifdef BL_USE_POINTLIB
#include <PointLooping.H>
#endif

template <class T>
BaseFab<T>::BaseFab ()
    : domain(Box()),
      pdomain(Box()),
      nvar(0),
      noffset(0),
      pnvar(0),
      numpts(0)
{ }

template <class T>
BaseFab<T>::BaseFab (const Box& bx,
                     int        n)
    : domain(bx),
      pdomain(bx)
{
    pnvar   = nvar = n;
    pnumpts = numpts = bx.numPts();
    noffset = 0;
    define();
}

#ifdef BL_USE_POINTLIB
#ifndef BL_CRAY_BUG_DEFARG
template <class T>
BaseFab<T>::BaseFab (const PointBaseFab<PointDomain,T>& pbf,
                     T                                  val)
{
    pdomain = domain = pbf.minimalBox();
    pnvar   = nvar = pbf.nComp();
    noffset = 0;
    pnumpts = numpts = domain.numPts();
    resize(domain, nvar);
    setVal(val);
    PointForAllCX(PointDomain,T,pbf)
    {
        this->operator()(ivR,nR) = pbfR;
    } EndPointForAll
}
#endif /*BL_CRAY_BUG_DEFARG*/
#endif /*BL_USE_POINTLIB*/

template <class T>
void
BaseFab<T>::define ()
{
    assert(nvar > 0);
    assert(numpts > 0);
    if (!dptr.ok())
    {
        AliasedDPtr<T>* ptr = new AliasedDPtr<T>(nvar*numpts);
        if (ptr == 0)
            BoxLib::OutOfMemory(__FILE__, __LINE__);
        dptr.setDPtr(ptr);
    }
}

//
// Alias.
//
template <class T>
BaseFab<T>::BaseFab (BaseFab<T>& fr,
                     const Box&  subb,
                     int         ns,
                     int         nc)
    : domain(subb),
      pdomain(fr.box())
{
    pnvar   = fr.nComp();
    pnumpts = fr.pnumpts;
    nvar    = nc;
    noffset = ns;
    numpts  = subb.numPts();

    AliasedDPtr<T>* p = new AliasedDPtr<T>((AliasedDPtr<T>*)fr.dptr.getDPtr(),
                                           nvar*numpts);
    if (p == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);

    dptr.setDPtr(p);
}

//
// This isn't inlined as it's virtual.
//

template <class T>
BaseFab<T>::~BaseFab ()
{}

template <class T>
void
BaseFab<T>::resize (const Box& b,
                    int        n)
{
    pnvar    = nvar   = n;
    pdomain  = domain = b;
    pnumpts  = numpts = domain.numPts();
    int size = nvar*numpts;
    if (!dptr.ok())
        define();
    else
        dptr->resize(size);
}

template <class T>
void
BaseFab<T>::clear ()
{
    if (dptr.ok())
        dptr->clear();
    domain = Box();
    nvar   = 0;
    numpts = 0;
}

// performCopy has been rewritten here so we can insert pragma's which
// will enhance vectorization on the Cray's.  The downside is that the
// code is greatly expanded and rather incompreshensible.
template <class T>
void
BaseFab<T>::performCopy (const BaseFab<T>& src,
                         const Box&        srcbox,
                         int               srccomp,
                         const Box&        destbox,
                         int               destcomp,
                         int               numcomp)
{
    assert(src.box().contains(srcbox));
    assert(box().contains(destbox));
    assert(destbox.sameSize(srcbox));
    assert(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    assert(destcomp >= 0 && destcomp+numcomp <= nComp());

#if (BL_SPACEDIM == 1)
{                                                                       
    assert((destcomp) >= 0 && (destcomp) + (numcomp) <= nComp());
    assert((srccomp) >= 0 && (srccomp) + (numcomp) <= (src).nComp());                
    Box _subbox_ = box();                                            
    _subbox_ &= destbox;                                                  
    assert(srcbox.sameSize(_subbox_));                                 
    if(_subbox_.ok()) {                                           
        const int *_th_plo = pLoVect();                          
        const int *_th_plen = pLength();                        
        const int *_x_plo = (src).pLoVect();                     
        const int *_x_plen = (src).pLength();                   
        const int *_subbox_lo = _subbox_.loVect();           
        const int *_subbox_len = _subbox_.length();         
        const int *_bx_lo = (srcbox).loVect();                 
        const int *_bx_len = (srcbox).length();               
        T* _th_p = dataPtr(destcomp);                          
        const T* _x_p  = (src).dataPtr(srccomp);              
        for(int _n = 0; _n < (numcomp); ++_n) {             
            int nR = _n + destcomp; nR += 0;                
            int nxR = _n + srccomp; nxR += 0;    
            T *_th_pp = _th_p                       
                + ((_subbox_lo[0] - _th_plo[0])    
                   + _n * _th_plen[0]);           
            const T *_x_pp = _x_p                
                + ((_bx_lo[0] - _x_plo[0])      
                   + _n * _x_plen[0]);         
#ifdef BL_ARCH_CRAY
#     pragma _CRI ivdep
#endif
            for(int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp) {      
                int iR = _i + _subbox_lo[0]; iR += 0;                  
                int ixR = _i + _bx_lo[0]; ixR += 0;           
                T &thisR = * _th_pp; const T & srcR = _x_pp[_i];
#elif (BL_SPACEDIM == 2)
{                                                                       
    assert((destcomp) >= 0 && (destcomp) + (numcomp) <= nComp());                       
    assert((srccomp) >= 0 && (srccomp) + (numcomp) <= (src).nComp());                
    Box _subbox_ = box();                                            
    _subbox_ &= destbox;                                                  
    assert(srcbox.sameSize(_subbox_));                                 
    if(_subbox_.ok()) {                                           
        const int *_th_plo = pLoVect();                          
        const int *_th_plen = pLength();                        
        const int *_x_plo = (src).pLoVect();                     
        const int *_x_plen = (src).pLength();                   
        const int *_subbox_lo = _subbox_.loVect();           
        const int *_subbox_len = _subbox_.length();         
        const int *_bx_lo = (srcbox).loVect();                 
        const int *_bx_len = (srcbox).length();               
        T* _th_p = dataPtr(destcomp);                          
        const T* _x_p  = (src).dataPtr(srccomp);              
        for(int _n = 0; _n < (numcomp); ++_n) {             
            int nR = _n + destcomp; nR += 0;                
            int nxR = _n + srccomp; nxR += 0;    
            for(int _j = 0; _j < _subbox_len[1]; ++_j) {                
                const int jR = _j + _subbox_lo[1];                     
                const int jxR = _j + _bx_lo[1];                   
                T *_th_pp = _th_p                                    
                    + ((_subbox_lo[0] - _th_plo[0])                 
                       + _th_plen[0]*(                             
                           (jR - _th_plo[1])                      
                           + _n * _th_plen[1]));                 
                const T *_x_pp = _x_p                           
                    + ((_bx_lo[0] - _x_plo[0])                 
                       + _x_plen[0]*(                         
                           (jxR - _x_plo[1])             
                           + _n * _x_plen[1]));             
#ifdef BL_ARCH_CRAY
#     pragma _CRI ivdep
#endif
                for(int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp) {  
                    int iR = _i + _subbox_lo[0]; iR += 0;              
                    int ixR = _i + _bx_lo[0]; ixR += 0; 
                    T &thisR = * _th_pp; const T & srcR = _x_pp[_i];
#elif (BL_SPACEDIM == 3)
{ 
    assert((destcomp) >= 0 && (destcomp) + (numcomp) <= nComp());                        
    assert((srccomp) >= 0 && (srccomp) + (numcomp) <= (src).nComp());                 
    Box _subbox_(box());                                              
    _subbox_ &= destbox;                                                   
    assert((srcbox).sameSize(_subbox_));                                
    if(_subbox_.ok()) {                                            
        const int *_th_plo = pLoVect();                           
        const int *_th_plen = pLength();                         
        const int *_x_plo = (src).pLoVect();                      
        const int *_x_plen = (src).pLength();                    
        const int *_subbox_lo = _subbox_.loVect();            
        const int *_subbox_len = _subbox_.length();          
        const int *_bx_lo = (srcbox).loVect();                  
        const int *_bx_len = (srcbox).length();                
        T* _th_p = dataPtr(destcomp);                           
        const T* _x_p  = (src).dataPtr(srccomp);               
        for(int _n = 0; _n < (numcomp); ++_n) {              
            int nR = _n + destcomp; nR += 0;                 
            int nxR = _n + srccomp; nxR += 0;     
            for(int _k = 0; _k < _subbox_len[2]; ++_k) {              
                const int kR = _k + _subbox_lo[2];                   
                const int kxR = _k + _bx_lo[2];                 
                for(int _j = 0; _j < _subbox_len[1]; ++_j) {       
                    const int jR = _j + _subbox_lo[1];            
                    const int jxR = _j + _bx_lo[1];         
                    T *_th_pp = _th_p                          
                        + ((_subbox_lo[0] - _th_plo[0])       
                           + _th_plen[0]*(                   
                               (jR - _th_plo[1])            
                               + _th_plen[1]*(             
                                   (kR - _th_plo[2])      
                                   + _n * _th_plen[2])));
                    const T *_x_pp = _x_p               
                        + ((_bx_lo[0] - _x_plo[0])     
                           + _x_plen[0]*(             
                               (jxR - _x_plo[1]) 
                               + _x_plen[1]*(       
                                   (kxR - _x_plo[2])               
                                   + _n * _x_plen[2])));              
#ifdef BL_ARCH_CRAY
#     pragma _CRI ivdep
#endif
                    for(int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp) {
                        int iR = _i + _subbox_lo[0]; iR += 0;          
                        int ixR = _i + _bx_lo[0]; ixR += 0;   
                        T &thisR = * _th_pp; const T & srcR = _x_pp[_i];
#endif
    {
        thisR = srcR;
    }
#if (BL_SPACEDIM == 1)
     }}}}
#elif (BL_SPACEDIM == 2)
     }}}}}
#elif (BL_SPACEDIM == 3)
     }}}}}}
#endif
}

#ifdef BL_USE_POINTLIB
#ifndef BL_CRAY_BUG_DEFARG
template <class T>
void
BaseFab<T>::performCopy (const PointBaseFab<PointDomain,T>& src,
                         const Box&                         srcbox,
                         int                                srccomp,
                         const Box&                         destbox,
                         int                                destcomp,
                         int                                numcomp)
{
    assert(box().contains(destbox));
    assert(destbox.sameSize(srcbox));
    assert(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    assert(destcomp >= 0 && destcomp+numcomp <= nComp());

    destcomp -= srccomp; // Adjust
    Box bx(domain);
    bx &= srcbox;
    bx &= destbox;
    if (bx.ok())
    {
        PointForAllCXBNN(PointDomain,T,src,bx,srccomp,numcomp)
        {
            this->operator()(ivR, destcomp + nR) = srcR;
        } EndPointForAll
    }
}
#endif /*!BL_CRAY_BUG_DEFARG*/
#endif /*BL_USE_POINTLIB*/

template <class T>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src,
                  const Box&        srcbox,
                  int               srccomp,
                  const Box&        destbox,
                  int               destcomp,
                  int               numcomp)
{
    assert(srcbox.sameSize(destbox));
    assert(src.box().contains(srcbox));
    assert(domain.contains(destbox));
    assert(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    assert(destcomp >= 0 && destcomp+numcomp <= nvar);
    performCopy(src,srcbox,srccomp,destbox,destcomp,numcomp);
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src)
{
    assert(nvar == src.nvar);
    assert(domain.sameType(src.domain));
    Box overlap(domain);
    overlap &= src.domain;
    if (overlap.ok())
        performCopy(src,overlap,0,overlap,0,nvar);
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src,
                  const Box&        destbox)
{
    assert(nvar == src.nvar);
    assert(domain.contains(destbox));
    Box overlap(destbox);
    overlap &= src.domain;
    if (overlap.ok())
        performCopy(src,overlap,0,overlap,0,nvar);
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::copy  (const BaseFab<T>& src,
                   int               srccomp,
                   int               destcomp,
                   int               numcomp)
{
    assert(srccomp >= 0 && srccomp + numcomp <= src.nvar);
    assert(destcomp >= 0 && destcomp + numcomp <= nvar);
    Box overlap(domain);
    overlap &= src.domain;
    if (overlap.ok())
        performCopy(src,overlap,srccomp,overlap,destcomp,numcomp);
    return *this;
}

template <class T>
void
BaseFab<T>::getVal  (T*             data,
                     const IntVect& pos,
                     int            n,
                     int            numcomp) const
{
    int loc      = domain.index(pos);
    int size     = domain.numPts();
    const T* ptr = dptr;
    assert(n >= 0 && n + numcomp <= nvar);
    for (int k = 0; k < numcomp; k++)
        data[k] = ptr[loc+(n+k)*size];
}

template <class T>
void
BaseFab<T>::performSetVal (T         val,
                           const Box& bx,
                           int        ns,
                           int        num)
{
    assert(domain.contains(bx));
    assert(ns >= 0 && ns + num <= nvar);
    ForAllThisBNN(T,bx,ns,num)
    {
        thisR = val;
    } EndFor
}

template <class T>
void
BaseFab<T>::setComplement (T          x,
                           const Box& b,
                           int        ns,
                           int        num)
{
    BoxList b_lst = boxDiff(domain,b);
    for (BoxListIterator bli(b_lst); bli; ++bli)
        performSetVal(x,bli(),ns,num);
}

