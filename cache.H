#ifndef _CACHE_H_
#define _CACHE_H_

//
// $Id: cache.H,v 1.7 1997-10-03 23:37:33 car Exp $
//

#include <interface.H>
#include <boundary.H>

#if 0
class boundary_cache 
{
 public:
  virtual ~boundary_cache() { }
  virtual void run() const = 0;
};
#endif

#ifdef HG_USE_CACHE
class copy_cache 
{
 public:
  // data are public so that users may set up their own calls
  int nsets;
  Real *dptr, *sptr;
  int *dstart, *sstart;
#if (BL_SPACEDIM == 2)
  int *dstrid, *sstrid;
  int *nvals;
#else
  int *dstrid1, *dstrid2, *sstrid1, *sstrid2;
  int *nvals1, *nvals2;
#endif

  copy_cache(int Nsets, Real *Dptr, Real *Sptr);

  // sync cache
  copy_cache(MultiFab& r, const level_interface& interface,
	     amr_boundary bdy);

  // border cache (only w == 1 case currently supported)
  copy_cache(MultiFab& r, const level_interface& interface,
	     amr_boundary bdy, int w);

  ~copy_cache()
  {
    delete [] dstart;
  }

#if (BL_SPACEDIM == 2)
  void set(int i, int Dstart, int Sstart, int Dstrid, int Sstrid, int Nvals);
#else
  void set(int i, int Dstart, int Sstart, int Dstrid1, int Dstrid2,
	   int Sstrid1, int Sstrid2, int Nvals1, int Nvals2);
#endif

  void run() const;
};
#endif

#ifdef HG_CONSTANT
class unroll_cache 
{
 public:
  // data are public so that users may set up their own calls
  int nsets;
  Real *ptr;
  int *start;
#if (BL_SPACEDIM == 2)
  int *strid;
#else
  int *strid1, *strid2;
#endif
  int *nvals;

  unroll_cache(int Nsets, Real *Ptr);
  unroll_cache(MultiFab& r);

  ~unroll_cache() 
  {
    delete [] start;
  }

#if (BL_SPACEDIM == 2)
  void set(int i, int Start, int Strid, int Nvals);
#else
  void set(int i, int Start, int Strid1, int Strid2, int Nvals);
#endif
};
#endif

#endif
