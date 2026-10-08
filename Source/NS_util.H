#ifndef NS_util_H
#define NS_util_H

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

namespace amrex
{
    template <class T>
    Gpu::DeviceVector<T>
    convertToDeviceVector ( Vector<T> v)
    {
        int ncomp = v.size();
        Gpu::DeviceVector<T> v_d(ncomp);
#ifdef AMREX_USE_GPU
        Gpu::htod_memcpy
#else
    std::memcpy
#endif
      (v_d.data(), v.data(), sizeof(T)*ncomp);

        return v_d;
    }


}

#endif
