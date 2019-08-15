
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <GODUNOV_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3
#define XVEL 1
#define YVEL 2
#define ZVEL 3
#define RHO  4

#define ALL  999

module godunov_3d_module

  use amrex_fort_module, only : rt=>amrex_real
  
  implicit none

  private
  
  public :: extrap_vel_to_faces, fort_estdt, fort_maxchng_velmag, &
            fort_test_umac_rho, &
            adv_forcing, sync_adv_forcing, &
            convscalminmax, consscalminmax, &
            fort_sum_tf_gp, fort_sum_tf_gp_visc, fort_sum_tf_divu, &
            fort_sum_tf_divu_visc, update_tf, update_aofs_tf, &
            update_aofs_tf_gp

contains
           

 subroutine extrap_vel_to_faces(lo,hi,&
       u,u_lo,u_hi,&
       ubc, tfx,tfx_lo,tfx_hi, umac,umac_lo,umac_hi, &
       vbc, tfy,tfy_lo,tfy_hi, vmac,vmac_lo,vmac_hi, &
       wbc, tfz,tfz_lo,tfz_hi, wmac,wmac_lo,wmac_hi, &
       corner_couple, &
       dt, dx, use_forces_in_trans, ppm_type)  bind(C,name="extrap_vel_to_faces")

    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate

    implicit none
    real(rt), intent(in) :: dt, dx(SDIM)
    integer,  intent(in) ::  ubc(SDIM,2),vbc(SDIM,2),wbc(SDIM,2), use_forces_in_trans, ppm_type, corner_couple
    integer,  intent(in), dimension(3) :: lo,hi,u_lo,u_hi,&
                                          tfx_lo, tfx_hi, tfy_lo, tfy_hi, tfz_lo, tfz_hi, &
                                          umac_lo, umac_hi, vmac_lo, vmac_hi, wmac_lo, wmac_hi

    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),3)
    real(rt), intent(in) :: tfx(tfx_lo(1):tfx_hi(1),tfx_lo(2):tfx_hi(2),tfx_lo(3):tfx_hi(3))
    real(rt), intent(inout) :: umac(umac_lo(1):umac_hi(1),umac_lo(2):umac_hi(2),umac_lo(3):umac_hi(3)) ! result
    real(rt), intent(in) :: tfy(tfy_lo(1):tfy_hi(1),tfy_lo(2):tfy_hi(2),tfy_lo(3):tfy_hi(3))
    real(rt), intent(inout) :: vmac(vmac_lo(1):vmac_hi(1),vmac_lo(2):vmac_hi(2),vmac_lo(3):vmac_hi(3)) ! result
    real(rt), intent(in) :: tfz(tfz_lo(1):tfz_hi(1),tfz_lo(2):tfz_hi(2),tfz_lo(3):tfz_hi(3))
    real(rt), intent(inout) :: wmac(wmac_lo(1):wmac_hi(1),wmac_lo(2):wmac_hi(2),wmac_lo(3):wmac_hi(3)) ! result

    integer, dimension(3) :: wklo,wkhi,uwlo,uwhi,vwlo,vwhi,wwlo,wwhi, &
                             eblo,ebhi,ebxhi,ebyhi,ebzhi,g2lo,g2hi
    real(rt), dimension(:,:,:), pointer, contiguous :: xlo,xhi,sx,xedge,uad
    real(rt), dimension(:,:,:), pointer, contiguous :: ylo,yhi,sy,yedge,vad
    real(rt), dimension(:,:,:), pointer, contiguous :: zlo,zhi,sz,zedge,wad
    real(rt), dimension(:,:,:), pointer, contiguous :: Imx,Ipx,sedgex
    real(rt), dimension(:,:,:), pointer, contiguous :: Imy,Ipy,sedgey
    real(rt), dimension(:,:,:), pointer, contiguous :: Imz,Ipz,sedgez
    real(rt), dimension(:,:,:), pointer, contiguous :: sm,sp,dsvl
    real(rt), dimension(:,:,:), pointer, contiguous :: xylo,xzlo,yxlo,yzlo,zxlo,zylo
    real(rt), dimension(:,:,:), pointer, contiguous :: xyhi,xzhi,yxhi,yzhi,zxhi,zyhi
    integer velpred
    parameter( velpred = 1 )
    
    ! Works space requirements:
    !  on wk=grow(bx,1): xlo,xhi,sx,ylo,yhi,sy,sm,sp,Imx,Ipx,Imy,Ipy,    uad,vad - ec(wk)
    !  on g2=grow(bx,2): dsvl
    !
    !  Grow cells for u,v:
    !  ppm:
    !    ppm=1 or 2, need 3 grow cells
    !  else:
    !    slope_order = 1, need 1 grow cell
    !    slope_order = 2, need 2 grow cells
    !    else , need 3 grow cells
    wklo = lo - 1
    wkhi = hi + 1
    call amrex_allocate(xlo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(xhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(sx, wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(xedge,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))

    call amrex_allocate(ylo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(yhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(sy, wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(yedge,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))

    call amrex_allocate(zlo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(sz, wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zedge,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    
    call amrex_allocate(xylo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(xzlo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(yxlo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(yzlo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zxlo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zylo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    
    call amrex_allocate(xyhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(xzhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(yxhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(yzhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zxhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zyhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))    

    uwlo = wklo
    uwhi = wkhi
    uwhi(1) = wkhi(1) + 1

    vwlo = wklo
    vwhi = wkhi
    vwhi(2) = vwhi(2) + 1

    wwlo = wklo
    wwhi = wkhi
    wwhi(3) = wwhi(3) + 1

    call amrex_allocate(uad,uwlo(1),uwhi(1),uwlo(2),uwhi(2),uwlo(3),uwhi(3))
    call amrex_allocate(vad,vwlo(1),vwhi(1),vwlo(2),vwhi(2),vwlo(3),vwhi(3))
    call amrex_allocate(wad,wwlo(1),wwhi(1),wwlo(2),wwhi(2),wwlo(3),wwhi(3))

    ! Note: Intel unhappy to pass pointers through subroutines if not allocated
    !     We just allocate something small (and still promise not to use it)
    if (ppm_type .eq. 0) then
       eblo = lo
       ebhi = lo
       call amrex_allocate(Imx,   eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(Ipx,   eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(Imy,   eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(Ipy,   eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(Imz,   eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(Ipz,   eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(sm,    eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(sp,    eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(sedgex,eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(sedgey,eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(sedgez,eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(dsvl,  eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))       
    else
    !if (ppm_type .gt. 0) then
       if (ppm_type .eq. 2) then
          eblo = lo - 2
          ebhi = hi + 2
       else
          eblo = lo - 1
          ebhi = hi + 1
       endif

       ebxhi = ebhi
       ebxhi(1) = ebxhi(1) + 1

       ebyhi = ebhi
       ebyhi(2) = ebyhi(2) + 1

       ebzhi = ebhi
       ebzhi(3) = ebzhi(3) + 1

       call amrex_allocate(Imx,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
       call amrex_allocate(Ipx,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))

       call amrex_allocate(Imy,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
       call amrex_allocate(Ipy,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))

       call amrex_allocate(Imz,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
       call amrex_allocate(Ipz,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))

       call amrex_allocate(sm,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
       call amrex_allocate(sp,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))

       call amrex_allocate(sedgex,eblo(1),ebxhi(1),eblo(2),ebxhi(2),eblo(3),ebxhi(3))
       call amrex_allocate(sedgey,eblo(1),ebyhi(1),eblo(2),ebyhi(2),eblo(3),ebyhi(3))
       call amrex_allocate(sedgez,eblo(1),ebzhi(1),eblo(2),ebzhi(2),eblo(3),ebzhi(3))

       g2lo = lo - 2
       g2hi = hi + 2
       call amrex_allocate(dsvl,g2lo(1),g2hi(1),g2lo(2),g2hi(2),g2lo(3),g2hi(3))
    endif

    ! get velocities that resolve upwind directions on faces used to compute transverse derivatives (uad,vad)
    ! Note that this is done only in this "pre-mac" situation, to get velocities on faces that will be projected.
    ! When advecting the other states, we will use the projected velocities, not these approximate versions.
    ! These face-centered arrays for each direction, dir, are needed on surroundingNodes(gBox(lo,hi),dir), where
    ! gBox is Box(lo,hi) grown in all directions but dir.
    call transvel(lo, hi,&
         u(u_lo(1),u_lo(2),u_lo(3),1),u_lo,u_hi,&
         uad,uwlo,uwhi,&
         xhi,wklo,wkhi,&
         sx,wklo,wkhi,&
         ubc,&
         Imx,wklo,wkhi,&
         Ipx,wklo,wkhi,&
         sedgex,eblo,ebxhi,&
         u(u_lo(1),u_lo(2),u_lo(3),2),u_lo,u_hi,&
         vad,vwlo,vwhi,&
         yhi,wklo,wkhi,&
         sy,wklo,wkhi,&
         vbc,&
         Imy,wklo,wkhi,&
         Ipy,wklo,wkhi,&
         sedgey,eblo,ebyhi,&
         u(u_lo(1),u_lo(2),u_lo(3),3),u_lo,u_hi,&
         wad,wwlo,wwhi,&
         zhi,wklo,wkhi,&
         sz,wklo,wkhi,&
         wbc,&
         Imz,wklo,wkhi,&
         Ipz,wklo,wkhi,&
         sedgez,eblo,ebzhi,&
         dsvl,g2lo,g2hi,&
         sm,wklo,wkhi,&
         sp,wklo,wkhi,&
         tfx,tfx_lo,tfx_hi,&
         tfy,tfy_lo,tfy_hi,&
         tfz,tfz_lo,tfz_hi,&
         dt, dx, use_forces_in_trans, ppm_type)

    ! Note that the final edge velocites are are resolved using the average of the velocities predicted from both
    ! sides (including the transverse terms resolved using uad,vad from above).

    ! get velocity on x-face, predict from cc u (because velpred=1 && n=XVEL, only predicts to xfaces)
    call estate_premac(lo,hi,&
         u(u_lo(1),u_lo(2),u_lo(3),1),u_lo,u_hi,&
         tfx,tfx_lo,tfx_hi,&
         u(u_lo(1),u_lo(2),u_lo(3),1),u_lo,u_hi,&
         xlo,wklo,wkhi,&
         xhi,wklo,wkhi,&
         sx,wklo,wkhi,&
         uad,uwlo,uwhi,&
         xedge,wklo,wkhi,&
         uad,uwlo,uwhi,& ! Unused uedge since velpred = 1
         umac,umac_lo,umac_hi,&
         Imx,wklo,wkhi,&
         Ipx,wklo,wkhi,&
         sedgex,eblo,ebxhi,&
         u(u_lo(1),u_lo(2),u_lo(3),2),u_lo,u_hi,&
         ylo,wklo,wkhi,&
         yhi,wklo,wkhi,&
         sy,wklo,wkhi,&
         vad,vwlo,vwhi,&
         yedge,wklo,wkhi,&
         vad,vwlo,vwhi,& ! Unused uedge since velpred = 1
         vmac,vmac_lo,vmac_hi,&
         Imy,wklo,wkhi,&
         Ipy,wklo,wkhi,&
         sedgey,eblo,ebyhi,&
         u(u_lo(1),u_lo(2),u_lo(3),3),u_lo,u_hi,&
         zlo,wklo,wkhi,&
         zhi,wklo,wkhi,&
         sz,wklo,wkhi,&
         wad,wwlo,wwhi,&
         zedge,wklo,wkhi,&
         wad,wwlo,wwhi,& ! Unused uedge since velpred = 1
         wmac,wmac_lo,wmac_hi,&
         Imz,wklo,wkhi,&
         Ipz,wklo,wkhi,&
         sedgez,eblo,ebzhi,&
         dsvl,g2lo,g2hi,&
         sm,wklo,wkhi,&
         sp,wklo,wkhi,&
         xylo,wklo,wkhi,&
         xzlo,wklo,wkhi,&
         yxlo,wklo,wkhi,&
         yzlo,wklo,wkhi,&
         zxlo,wklo,wkhi,&
         zylo,wklo,wkhi,&
         xyhi,wklo,wkhi,&
         xzhi,wklo,wkhi,&
         yxhi,wklo,wkhi,&
         yzhi,wklo,wkhi,&
         zxhi,wklo,wkhi,&
         zyhi,wklo,wkhi,&
         corner_couple,&
         ubc, dt, dx, XVEL, 1, velpred, use_forces_in_trans, ppm_type)

    ! get velocity on y-face, predict from cc v (because velpred=1 && n=YVEL, only predicts to yfaces)
    call estate_premac(lo,hi,&
         u(u_lo(1),u_lo(2),u_lo(3),2),u_lo,u_hi,&
         tfy,tfy_lo,tfy_hi,&
         u(u_lo(1),u_lo(2),u_lo(3),1),u_lo,u_hi,&
         xlo,wklo,wkhi,&
         xhi,wklo,wkhi,&
         sx,wklo,wkhi,&
         uad,uwlo,uwhi,&
         xedge,wklo,wkhi,&
         uad,uwlo,uwhi,& ! Unused uedge since velpred = 1
         umac,umac_lo,umac_hi,&
         Imx,wklo,wkhi,&
         Ipx,wklo,wkhi,&
         sedgex,eblo,ebxhi,&
         u(u_lo(1),u_lo(2),u_lo(3),2),u_lo,u_hi,&
         ylo,wklo,wkhi,&
         yhi,wklo,wkhi,&
         sy,wklo,wkhi,&
         vad,vwlo,vwhi,&
         yedge,wklo,wkhi,&
         vad,vwlo,vwhi,& ! Unused uedge since velpred = 1
         vmac,vmac_lo,vmac_hi,&
         Imy,wklo,wkhi,&
         Ipy,wklo,wkhi,&
         sedgey,eblo,ebyhi,&
         u(u_lo(1),u_lo(2),u_lo(3),3),u_lo,u_hi,&
         zlo,wklo,wkhi,&
         zhi,wklo,wkhi,&
         sz,wklo,wkhi,&
         wad,wwlo,wwhi,&
         zedge,wklo,wkhi,&
         wad,wwlo,wwhi,& ! Unused uedge since velpred = 1
         wmac,wmac_lo,wmac_hi,&
         Imz,wklo,wkhi,&
         Ipz,wklo,wkhi,&
         sedgez,eblo,ebzhi,&
         dsvl,g2lo,g2hi,&
         sm,wklo,wkhi,&
         sp,wklo,wkhi,&
         xylo,wklo,wkhi,&
         xzlo,wklo,wkhi,&
         yxlo,wklo,wkhi,&
         yzlo,wklo,wkhi,&
         zxlo,wklo,wkhi,&
         zylo,wklo,wkhi,&
         xyhi,wklo,wkhi,&
         xzhi,wklo,wkhi,&
         yxhi,wklo,wkhi,&
         yzhi,wklo,wkhi,&
         zxhi,wklo,wkhi,&
         zyhi,wklo,wkhi,&
         corner_couple,&
         vbc, dt, dx, YVEL, 1, velpred, use_forces_in_trans, ppm_type)

   ! get velocity on z-face, predict from cc w (because velpred=1 && n=ZVEL, only predicts to zfaces)
    call estate_premac(lo,hi,&
         u(u_lo(1),u_lo(2),u_lo(3),3),u_lo,u_hi,&
         tfz,tfz_lo,tfz_hi,&
         u(u_lo(1),u_lo(2),u_lo(3),1),u_lo,u_hi,&
         xlo,wklo,wkhi,&
         xhi,wklo,wkhi,&
         sx,wklo,wkhi,&
         uad,uwlo,uwhi,&
         xedge,wklo,wkhi,&
         uad,uwlo,uwhi,& ! Unused uedge since velpred = 1
         umac,umac_lo,umac_hi,&
         Imx,wklo,wkhi,&
         Ipx,wklo,wkhi,&
         sedgex,eblo,ebxhi,&
         u(u_lo(1),u_lo(2),u_lo(3),2),u_lo,u_hi,&
         ylo,wklo,wkhi,&
         yhi,wklo,wkhi,&
         sy,wklo,wkhi,&
         vad,vwlo,vwhi,&
         yedge,wklo,wkhi,&
         vad,vwlo,vwhi,& ! Unused uedge since velpred = 1
         vmac,vmac_lo,vmac_hi,&
         Imy,wklo,wkhi,&
         Ipy,wklo,wkhi,&
         sedgey,eblo,ebyhi,&
         u(u_lo(1),u_lo(2),u_lo(3),3),u_lo,u_hi,&
         zlo,wklo,wkhi,&
         zhi,wklo,wkhi,&
         sz,wklo,wkhi,&
         wad,wwlo,wwhi,&
         zedge,wklo,wkhi,&
         wad,wwlo,wwhi,& ! Unused uedge since velpred = 1
         wmac,wmac_lo,wmac_hi,&
         Imz,wklo,wkhi,&
         Ipz,wklo,wkhi,&
         sedgez,eblo,ebzhi,&
         dsvl,g2lo,g2hi,&
         sm,wklo,wkhi,&
         sp,wklo,wkhi,&
         xylo,wklo,wkhi,&
         xzlo,wklo,wkhi,&
         yxlo,wklo,wkhi,&
         yzlo,wklo,wkhi,&
         zxlo,wklo,wkhi,&
         zylo,wklo,wkhi,&
         xyhi,wklo,wkhi,&
         xzhi,wklo,wkhi,&
         yxhi,wklo,wkhi,&
         yzhi,wklo,wkhi,&
         zxhi,wklo,wkhi,&
         zyhi,wklo,wkhi,&
         corner_couple,&
         wbc, dt, dx, ZVEL, 1, velpred, use_forces_in_trans, ppm_type)

    call amrex_deallocate(xedge)
    call amrex_deallocate(yedge)
    call amrex_deallocate(zedge)
    call amrex_deallocate(xylo)
    call amrex_deallocate(xzlo)
    call amrex_deallocate(yxlo)
    call amrex_deallocate(yzlo)
    call amrex_deallocate(zxlo)
    call amrex_deallocate(zylo)
    call amrex_deallocate(xyhi)
    call amrex_deallocate(xzhi)
    call amrex_deallocate(yxhi)
    call amrex_deallocate(yzhi)
    call amrex_deallocate(zxhi)
    call amrex_deallocate(zyhi)
    
    call amrex_deallocate(xlo)
    call amrex_deallocate(xhi)
    call amrex_deallocate(sx)
    call amrex_deallocate(ylo)
    call amrex_deallocate(yhi)
    call amrex_deallocate(sy)
    call amrex_deallocate(zlo)
    call amrex_deallocate(zhi)
    call amrex_deallocate(sz)
    call amrex_deallocate(uad)
    call amrex_deallocate(vad)
    call amrex_deallocate(wad)
    !if (ppm_type .gt. 0) then
       call amrex_deallocate(Imx)
       call amrex_deallocate(Ipx)
       call amrex_deallocate(Imy)
       call amrex_deallocate(Ipy)
       call amrex_deallocate(Imz)
       call amrex_deallocate(Ipz)
       call amrex_deallocate(sm)
       call amrex_deallocate(sp)
       call amrex_deallocate(sedgex)
       call amrex_deallocate(sedgey)
       call amrex_deallocate(sedgez)
       call amrex_deallocate(dsvl)
    !endif

  end subroutine extrap_vel_to_faces


  subroutine extrap_state_to_faces(lo,hi,&
       s,s_lo,s_hi,nc,              tf, tf_lo,tf_hi,              divu,divu_lo,divu_hi,&
       umac,umac_lo,umac_hi,        xstate,xstate_lo,xstate_hi,&
       vmac,vmac_lo,vmac_hi,        ystate,ystate_lo,ystate_hi,&
       wmac,wmac_lo,wmac_hi,        zstate,zstate_lo,zstate_hi,&
       corner_couple, &
       dt, dx, bc, state_ind, use_forces_in_trans, ppm_type, iconserv)  bind(C,name="extrap_state_to_faces")
  
    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
  
    implicit none
    integer, intent(in) ::  nc, bc(SDIM,2,nc), state_ind, use_forces_in_trans, ppm_type, iconserv(nc), corner_couple
    integer, dimension(3), intent(in) :: lo,hi,s_lo,s_hi,tf_lo,tf_hi,&
                              divu_lo,divu_hi,xstate_lo,xstate_hi,ystate_lo,ystate_hi,zstate_lo,zstate_hi, &
                              umac_lo,umac_hi,vmac_lo,vmac_hi,wmac_lo,wmac_hi
  
    real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc)
    real(rt), intent(in) :: tf(tf_lo(1):tf_hi(1),tf_lo(2):tf_hi(2),tf_lo(3):tf_hi(3),nc)
    real(rt), intent(in) :: divu(divu_lo(1):divu_hi(1),divu_lo(2):divu_hi(2),divu_lo(3):divu_hi(3))
    real(rt), intent(in) :: umac(umac_lo(1):umac_hi(1),umac_lo(2):umac_hi(2),umac_lo(3):umac_hi(3))
    real(rt), intent(inout) :: xstate(xstate_lo(1):xstate_hi(1),xstate_lo(2):xstate_hi(2),xstate_lo(3):xstate_hi(3),nc) ! result
    real(rt), intent(in) :: vmac(vmac_lo(1):vmac_hi(1),vmac_lo(2):vmac_hi(2),vmac_lo(3):vmac_hi(3))
    real(rt), intent(inout) :: ystate(ystate_lo(1):ystate_hi(1),ystate_lo(2):ystate_hi(2),ystate_lo(3):ystate_hi(3),nc) ! result
    real(rt), intent(in) :: wmac(wmac_lo(1):wmac_hi(1),wmac_lo(2):wmac_hi(2),wmac_lo(3):wmac_hi(3))
    real(rt), intent(inout) :: zstate(zstate_lo(1):zstate_hi(1),zstate_lo(2):zstate_hi(2),zstate_lo(3):zstate_hi(3),nc) ! result
  
    real(rt), intent(in) :: dt, dx(SDIM)
  
    integer, dimension(3) :: wklo,wkhi,eblo,ebhi,ebxhi,ebyhi,ebzhi,g2lo,g2hi
    real(rt), dimension(:,:,:), pointer, contiguous :: xlo,xhi,sx,xedge
    real(rt), dimension(:,:,:), pointer, contiguous :: ylo,yhi,sy,yedge
    real(rt), dimension(:,:,:), pointer, contiguous :: zlo,zhi,sz,zedge
    real(rt), dimension(:,:,:), pointer, contiguous :: Imx,Ipx,sedgex
    real(rt), dimension(:,:,:), pointer, contiguous :: Imy,Ipy,sedgey
    real(rt), dimension(:,:,:), pointer, contiguous :: Imz,Ipz,sedgez
    real(rt), dimension(:,:,:), pointer, contiguous :: sm,sp,dsvl
    real(rt), dimension(:,:,:), pointer, contiguous :: xylo,xzlo,yxlo,yzlo,zxlo,zylo
    real(rt), dimension(:,:,:), pointer, contiguous :: xyhi,xzhi,yxhi,yzhi,zxhi,zyhi
  
    ! Works space requirements:
    !  on wk=grow(bx,1): xlo,xhi,sx,ylo,yhi,sy,sm,sp,Imx,Ipx,Imy,Ipy
    !  on g2=grow(bx,2): dsvl
    !
    !  Grow cells for s:
    !  ppm:
    !    ppm=1 or 2, need 3 grow cells
    !  else:
    !    slope_order = 1, need 1 grow cell
    !    slope_order = 2, need 2 grow cells
    !    else , need 3 grow cells
    wklo = lo - 1
    wkhi = hi + 1
    call amrex_allocate(xlo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(xhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(sx,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(xedge,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(ylo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(yhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(sy,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(yedge,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zlo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(sz,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zedge,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
  
    call amrex_allocate(xylo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(xzlo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(yxlo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(yzlo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zxlo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zylo,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    
    call amrex_allocate(xyhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(xzhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(yxhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(yzhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zxhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
    call amrex_allocate(zyhi,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))

    ! Note: Intel unhappy to pass pointers through subroutines if not allocated
    !     We just allocate something small (and still promise not to use it)
    if (ppm_type .eq. 0) then
       eblo = lo
       ebhi = lo
       call amrex_allocate(Imx,   eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(Ipx,   eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(Imy,   eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(Ipy,   eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(Imz,   eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(Ipz,   eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(sm,    eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(sp,    eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(sedgex,eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(sedgey,eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(sedgez,eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))
       call amrex_allocate(dsvl,  eblo(1),ebhi(1),eblo(2),ebhi(2),eblo(3),ebhi(3))       
    else
       !if (ppm_type .gt. 0) then
       if (ppm_type .eq. 2) then
          eblo = lo - 2
          ebhi = hi + 2
       else
          eblo = lo - 1
          ebhi = hi + 1
       endif
  
       ebxhi = ebhi
       ebxhi(1) = ebxhi(1) + 1
       ebyhi = ebhi
       ebyhi(2) = ebyhi(2) + 1
       ebzhi = ebhi
       ebzhi(3) = ebzhi(3) + 1
       call amrex_allocate(Imx,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
       call amrex_allocate(Ipx,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
       call amrex_allocate(Imy,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
       call amrex_allocate(Ipy,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
       call amrex_allocate(Imz,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
       call amrex_allocate(Ipz,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
       call amrex_allocate(sm,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
       call amrex_allocate(sp,wklo(1),wkhi(1),wklo(2),wkhi(2),wklo(3),wkhi(3))
       call amrex_allocate(sedgex,eblo(1),ebxhi(1),eblo(2),ebxhi(2),eblo(3),ebxhi(3))
       call amrex_allocate(sedgey,eblo(1),ebyhi(1),eblo(2),ebyhi(2),eblo(3),ebyhi(3))
       call amrex_allocate(sedgez,eblo(1),ebzhi(1),eblo(2),ebzhi(2),eblo(3),ebzhi(3))
       g2lo = lo - 2
       g2hi = hi + 2
       call amrex_allocate(dsvl,g2lo(1),g2hi(1),g2lo(2),g2hi(2),g2lo(3),g2hi(3))
    endif
  
    call estate_fpu(lo,hi,&
         s,s_lo,s_hi,&
         tf,tf_lo,tf_hi,&
         divu,divu_lo,divu_hi,&
         xlo,wklo,wkhi,&
         xhi,wklo,wkhi,&
         sx,wklo,wkhi,&
         xedge,wklo,wkhi,&
         umac,umac_lo,umac_hi,&
         xstate,xstate_lo,xstate_hi,&
         Imx,wklo,wkhi,&
         Ipx,wklo,wkhi,&
         sedgex,eblo,ebxhi,&
         ylo,wklo,wkhi,&
         yhi,wklo,wkhi,&
         sy,wklo,wkhi,&
         yedge,wklo,wkhi,&
         vmac,vmac_lo,vmac_hi,&
         ystate,ystate_lo,ystate_hi,&
         Imy,wklo,wkhi,&
         Ipy,wklo,wkhi,&
         sedgey,eblo,ebyhi,&
         zlo,wklo,wkhi,&
         zhi,wklo,wkhi,&
         sz,wklo,wkhi,&
         zedge,wklo,wkhi,&
         wmac,wmac_lo,wmac_hi,&
         zstate,zstate_lo,zstate_hi,&
         Imz,wklo,wkhi,&
         Ipz,wklo,wkhi,&
         sedgez,eblo,ebzhi,&
         dsvl,g2lo,g2hi,&
         sm,wklo,wkhi,&
         sp,wklo,wkhi,&
         xylo,wklo,wkhi,&
         xzlo,wklo,wkhi,&
         yxlo,wklo,wkhi,&
         yzlo,wklo,wkhi,&
         zxlo,wklo,wkhi,&
         zylo,wklo,wkhi,&
         xyhi,wklo,wkhi,&
         xzhi,wklo,wkhi,&
         yxhi,wklo,wkhi,&
         yzhi,wklo,wkhi,&
         zxhi,wklo,wkhi,&
         zyhi,wklo,wkhi,&
         corner_couple,&
         bc, dt, dx, state_ind, nc, use_forces_in_trans, iconserv, ppm_type)
    !
    
    call amrex_deallocate(xedge)
    call amrex_deallocate(yedge)
    call amrex_deallocate(zedge)
    call amrex_deallocate(xylo)
    call amrex_deallocate(xzlo)
    call amrex_deallocate(yxlo)
    call amrex_deallocate(yzlo)
    call amrex_deallocate(zxlo)
    call amrex_deallocate(zylo)
    call amrex_deallocate(xyhi)
    call amrex_deallocate(xzhi)
    call amrex_deallocate(yxhi)
    call amrex_deallocate(yzhi)
    call amrex_deallocate(zxhi)
    call amrex_deallocate(zyhi)
    
    call amrex_deallocate(xlo)
    call amrex_deallocate(xhi)
    call amrex_deallocate(sx)
    call amrex_deallocate(ylo)
    call amrex_deallocate(yhi)
    call amrex_deallocate(sy)
    call amrex_deallocate(zlo)
    call amrex_deallocate(zhi)
    call amrex_deallocate(sz)
    !if (ppm_type .gt. 0) then
       call amrex_deallocate(Imx)
       call amrex_deallocate(Ipx)
       call amrex_deallocate(Imy)
       call amrex_deallocate(Ipy)
       call amrex_deallocate(Imz)
       call amrex_deallocate(Ipz)
       call amrex_deallocate(sm)
       call amrex_deallocate(sp)
       call amrex_deallocate(sedgex)
       call amrex_deallocate(sedgey)
       call amrex_deallocate(sedgez)
       call amrex_deallocate(dsvl)
    !endif
  
  end subroutine extrap_state_to_faces

  subroutine fort_estdt ( &
          vel,DIMS(vel), &
          tforces,DIMS(tf), &
          rho,DIMS(rho), &
          lo,hi,dt,dx,cfl,u_max) bind(C,name="fort_estdt")
!c 
!c     ----------------------------------------------------------
!c     estimate the timestep for this grid and scale by CFL number
!c     This routine sets dt as dt = dt_est*cfl where
!c     dt_est is estimated from the actual velocities and their 
!c     total forcing
!c     ----------------------------------------------------------
!c
      implicit none
      integer i, j, k
      real(rt)  u, v, w
      real(rt)  small
      real(rt)  dt_start
      real(rt)  tforce1,tforce2,tforce3
      integer lo(SDIM), hi(SDIM)
      real(rt)  dt,dx(SDIM),cfl,u_max(SDIM)

      integer DIMDEC(vel)
      integer DIMDEC(rho)
      integer DIMDEC(tf)

      real(rt)  vel(DIMV(vel),SDIM)
      real(rt)  rho(DIMV(rho))
      real(rt)  tforces(DIMV(tf),SDIM), irho

      PARAMETER(small = 1.0D-8, dt_start = 1.0D+20)

      u       = zero
      v       = zero
      w       = zero
      tforce1 = zero
      tforce2 = zero
      tforce3 = zero

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               irho = 1.0d0/rho(i,j,k)
               u = max(u,abs(vel(i,j,k,1)))	
               v = max(v,abs(vel(i,j,k,2)))
               w = max(w,abs(vel(i,j,k,3)))
               tforce1 = max(tforce1,abs(tforces(i,j,k,1)*irho))
               tforce2 = max(tforce2,abs(tforces(i,j,k,2)*irho))
               tforce3 = max(tforce3,abs(tforces(i,j,k,3)*irho))
            end do
         end do
      end do

      u_max(1) = u
      u_max(2) = v
      u_max(3) = w

      dt = dt_start

      if (u .gt. small) dt = min(dt,dx(1)/u)
      if (v .gt. small) dt = min(dt,dx(2)/v)
      if (w .gt. small) dt = min(dt,dx(3)/w)

      if (tforce1 .gt. small) dt = min(dt,sqrt(two*dx(1)/tforce1))
      if (tforce2 .gt. small) dt = min(dt,sqrt(two*dx(2)/tforce2))
      if (tforce3 .gt. small) dt = min(dt,sqrt(two*dx(3)/tforce3))

      if (dt .eq. dt_start) dt = min(dx(1),dx(2),dx(3))

      dt = dt*cfl

      end subroutine fort_estdt

      subroutine fort_maxchng_velmag ( &
          old_vel,DIMS(old_vel), &
          new_vel,DIMS(new_vel), &
          lo,hi,max_change) bind(C,name="fort_maxchng_velmag")
!c 
!c     ----------------------------------------------------------
!c     Given the velocity field at the previous and current time steps
!c     (old_vel and new_vel, respectively), find the largest change in
!c     velocity magnitude between the two.
!c     ----------------------------------------------------------
!c
      implicit none
      real(rt)   old_velmag, new_velmag
      integer  i, j, k
      integer  lo(SDIM), hi(SDIM)
      real(rt)   max_change

      integer DIMDEC(old_vel)
      integer DIMDEC(new_vel)

      real(rt)  old_vel(DIMV(old_vel),SDIM)
      real(rt)  new_vel(DIMV(new_vel),SDIM)

      max_change = zero

      do k = lo(3), hi(3)
          do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                old_velmag = sqrt(old_vel(i,j,k,1)**2 + &
                                 old_vel(i,j,k,2)**2 + &
                                 old_vel(i,j,k,3)**2)
                new_velmag = sqrt(new_vel(i,j,k,1)**2 + &
                                 new_vel(i,j,k,2)**2 + &
                                 new_vel(i,j,k,3)**2)
                max_change = max(max_change, abs(new_velmag - old_velmag))
            end do
          end do
      end do

      end subroutine fort_maxchng_velmag

      subroutine fort_test_umac_rho ( &
          umac,DIMS(umac), &
          vmac,DIMS(vmac), &
          wmac,DIMS(wmac), &
          rho,DIMS(rho), &
          lo,hi,dt,dx,cflmax,u_max) &
          bind(C,name="fort_test_umac_rho")
!c
!c     This subroutine computes the extrema of the density
!c     and mac edge velocities
!c
      implicit none
      integer lo(SDIM),hi(SDIM)
      real(rt)  dt, dx(SDIM), u_max(SDIM), cflmax
      integer imin, imax, jmin, jmax, kmin, kmax
      integer i, j, k
      real(rt)  hx, hy, hz
      real(rt)  umax, vmax, wmax, rhomax
      real(rt)  umin, vmin, wmin, rhomin

      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(wmac)
      integer DIMDEC(rho)

      real(rt)  umac(DIMV(umac))
      real(rt)  vmac(DIMV(vmac))
      real(rt)  wmac(DIMV(wmac))
      real(rt)  rho(DIMV(rho))

      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)
      umax = -1.d200
      vmax = -1.d200
      wmax = -1.d200
      umin =  1.d200
      vmin =  1.d200
      wmin =  1.d200
      rhomax = -1.d200
      rhomin =  1.d200

      do k = kmin, kmax
         do j = jmin, jmax
            do i = imin, imax+1
               umax = max(umax,umac(i,j,k))
               umin = min(umin,umac(i,j,k))
            end do
         end do
      end do

      do k = kmin, kmax
         do j = jmin, jmax+1
            do i = imin, imax
               vmax = max(vmax,vmac(i,j,k))
               vmin = min(vmin,vmac(i,j,k))
            end do
         end do
      end do

      do k = kmin, kmax+1
         do j = jmin, jmax
            do i = imin, imax
               wmax = max(wmax,wmac(i,j,k))
               wmin = min(wmin,wmac(i,j,k))
            end do
         end do
      end do

      do k = kmin, kmax
         do j = jmin, jmax
            do i = imin, imax
               rhomax = max(rhomax,rho(i,j,k))
               rhomin = min(rhomin,rho(i,j,k))
            end do
         end do
      end do

      u_max(1) = max(abs(umax), abs(umin))
      u_max(2) = max(abs(vmax), abs(vmin))
      u_max(3) = max(abs(wmax), abs(wmin))
      cflmax   = dt*max(u_max(1)/hx,u_max(2)/hy,u_max(3)/hz)

      write(6,1000) umax,umin,u_max(1)
      write(6,1001) vmax,vmin,u_max(2)
      write(6,1002) wmax,wmin,u_max(3)
      write(6,1003) rhomax,rhomin

 1000 format('UMAC MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1001 format('VMAC MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1002 format('WMAC MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1003 format('RHO  MAX/MIN      ',e21.14,2x,e21.14)

      end subroutine fort_test_umac_rho

      subroutine transvel (lo, hi, & 
           u,u_lo,u_hi,&
           ulo,ulo_lo,ulo_hi,&
           uhi,uhi_lo,uhi_hi,&
           sx,sx_lo,sx_hi,&
           ubc,&
           Imx,Imx_lo,Imx_hi,&
           Ipx,Ipx_lo,Ipx_hi,&
           sedgex,sedgex_lo,sedgex_hi,&
           v,v_lo,v_hi,&
           vlo,vlo_lo,vlo_hi,&
           vhi,vhi_lo,vhi_hi,&
           sy,sy_lo,sy_hi,&
           vbc,&
           Imy,Imy_lo,Imy_hi,&
           Ipy,Ipy_lo,Ipy_hi,&
           sedgey,sedgey_lo,sedgey_hi,&
           w,w_lo,w_hi,&
           wlo,wlo_lo,wlo_hi,&
           whi,whi_lo,whi_hi,&
           sz,sz_lo,sz_hi,&
           wbc,&
           Imz,Imz_lo,Imz_hi,&
           Ipz,Ipz_lo,Ipz_hi,&
           sedgez,sedgez_lo,sedgez_hi,&
           dsvl,dsvl_lo,dsvl_hi,&
           sm,sm_lo,sm_hi,&
           sp,sp_lo,sp_hi,&
           tfx,tfx_lo,tfx_hi,&
           tfy,tfy_lo,tfy_hi,&
           tfz,tfz_lo,tfz_hi,&
           dt, dx, use_minion, ppm_type)

!c
!c     This subroutine computes the advective velocities used in
!c     the transverse derivatives of the Godunov box
!c
      implicit none
      integer, intent(in) ::  ubc(SDIM,2),vbc(SDIM,2),wbc(SDIM,2), use_minion, ppm_type
      integer, dimension(3), intent(in) :: lo,hi,&
           u_lo,u_hi,ulo_lo,ulo_hi,uhi_lo,uhi_hi,sx_lo,sx_hi,&
           Imx_lo,Imx_hi,Ipx_lo,Ipx_hi,sedgex_lo,sedgex_hi,&
           v_lo,v_hi,vlo_lo,vlo_hi,vhi_lo,vhi_hi,sy_lo,sy_hi,&
           Imy_lo,Imy_hi,Ipy_lo,Ipy_hi,sedgey_lo,sedgey_hi,&
           w_lo,w_hi,wlo_lo,wlo_hi,whi_lo,whi_hi,sz_lo,sz_hi,&
           Imz_lo,Imz_hi,Ipz_lo,Ipz_hi,sedgez_lo,sedgez_hi,&
           dsvl_lo,dsvl_hi,sm_lo,sm_hi,sp_lo,sp_hi, &
           tfx_lo,tfx_hi,tfy_lo,tfy_hi,tfz_lo,tfz_hi
          

      real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
      real(rt), intent(inout) :: ulo(ulo_lo(1):ulo_hi(1),ulo_lo(2):ulo_hi(2),ulo_lo(3):ulo_hi(3))
      real(rt), intent(inout) :: uhi(uhi_lo(1):uhi_hi(1),uhi_lo(2):uhi_hi(2),uhi_lo(3):uhi_hi(3))
      real(rt), intent(inout) :: sx(sx_lo(1):sx_hi(1),sx_lo(2):sx_hi(2),sx_lo(3):sx_hi(3))
      real(rt), intent(inout) :: Imx(Imx_lo(1):Imx_hi(1),Imx_lo(2):Imx_hi(2),Imx_lo(3):Imx_hi(3))
      real(rt), intent(inout) :: Ipx(Ipx_lo(1):Ipx_hi(1),Ipx_lo(2):Ipx_hi(2),Ipx_lo(3):Ipx_hi(3))
      real(rt), intent(inout) :: sedgex(sedgex_lo(1):sedgex_hi(1),sedgex_lo(2):sedgex_hi(2),sedgex_lo(3):sedgex_hi(3))
      real(rt), intent(in) :: tfx(tfx_lo(1):tfx_hi(1),tfx_lo(2):tfx_hi(2),tfx_lo(3):tfx_hi(3))

      real(rt), intent(in) :: v(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
      real(rt), intent(inout) :: vlo(vlo_lo(1):vlo_hi(1),vlo_lo(2):vlo_hi(2),vlo_lo(3):vlo_hi(3))
      real(rt), intent(inout) :: vhi(vhi_lo(1):vhi_hi(1),vhi_lo(2):vhi_hi(2),vhi_lo(3):vhi_hi(3))
      real(rt), intent(inout) :: sy(sy_lo(1):sy_hi(1),sy_lo(2):sy_hi(2),sy_lo(3):sy_hi(3))
      real(rt), intent(inout) :: Imy(Imy_lo(1):Imy_hi(1),Imy_lo(2):Imy_hi(2),Imy_lo(3):Imy_hi(3))
      real(rt), intent(inout) :: Ipy(Ipy_lo(1):Ipy_hi(1),Ipy_lo(2):Ipy_hi(2),Ipy_lo(3):Ipy_hi(3))
      real(rt), intent(inout) :: sedgey(sedgey_lo(1):sedgey_hi(1),sedgey_lo(2):sedgey_hi(2),sedgey_lo(3):sedgey_hi(3))
      real(rt), intent(in) :: tfy(tfy_lo(1):tfy_hi(1),tfy_lo(2):tfy_hi(2),tfy_lo(3):tfy_hi(3))

      real(rt), intent(in) :: w(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
      real(rt), intent(inout) :: wlo(wlo_lo(1):wlo_hi(1),wlo_lo(2):wlo_hi(2),wlo_lo(3):wlo_hi(3))
      real(rt), intent(inout) :: whi(whi_lo(1):whi_hi(1),whi_lo(2):whi_hi(2),whi_lo(3):whi_hi(3))
      real(rt), intent(inout) :: sz(sz_lo(1):sz_hi(1),sz_lo(2):sz_hi(2),sz_lo(3):sz_hi(3))
      real(rt), intent(inout) :: Imz(Imz_lo(1):Imz_hi(1),Imz_lo(2):Imz_hi(2),Imz_lo(3):Imz_hi(3))
      real(rt), intent(inout) :: Ipz(Ipz_lo(1):Ipz_hi(1),Ipz_lo(2):Ipz_hi(2),Ipz_lo(3):Ipz_hi(3))
      real(rt), intent(inout) :: sedgez(sedgez_lo(1):sedgez_hi(1),sedgez_lo(2):sedgez_hi(2),sedgez_lo(3):sedgez_hi(3))
      real(rt), intent(in) :: tfz(tfz_lo(1):tfz_hi(1),tfz_lo(2):tfz_hi(2),tfz_lo(3):tfz_hi(3))

      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2),dsvl_lo(3):dsvl_hi(3))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2),sp_lo(3):sp_hi(3))

      integer :: i,j, k, imin,jmin,kmin,imax,jmax,kmax
      real(rt) :: hx, hy, hz, dt, dth, dthx, dthy, dthz, dx(SDIM)
      real(rt) :: eps,eps_for_bc, val, tst, dt3
      logical :: ltm
      parameter( eps        = 1.0D-6 )
      parameter( eps_for_bc = 1.0D-10 )

      dth  = half*dt
      dt3  = dt / 3.d0
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      dthz = half*dt / dx(3)
      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

!c
!c     =============== THE SCREWY ORDER is to maximize comparability
!c     with the old fortran
!c     --------------------------------------------------------------
!c     compute the x transverse velocities
!c     --------------------------------------------------------------
!c     --------------------------------------------------------------
!c     compute the y transverse velocities
!c     --------------------------------------------------------------
!c     --------------------------------------------------------------
!c     compute the z transverse velocities
!c     --------------------------------------------------------------
!c
      if (ppm_type .gt. 0) then
         call ppm(lo,hi,&
              u,u_lo,u_hi,&
              u,u_lo,u_hi,&
              v,v_lo,v_hi,&
              w,w_lo,w_hi,&
              Ipx,Ipx_lo,Ipx_hi,&
              Imx,Imx_lo,Imx_hi,&
              Ipy,Ipy_lo,Ipy_hi,&
              Imy,Imy_lo,Imy_hi,&
              Ipz,Ipz_lo,Ipz_hi,&
              Imz,Imz_lo,Imz_hi,&
              sm,sm_lo,sm_hi,&
              sp,sp_lo,sp_hi,&
              dsvl,dsvl_lo,dsvl_hi,&
              sedgex,sedgex_lo,sedgex_hi,&
              sedgey,sedgey_lo,sedgey_hi,&
              sedgez,sedgez_lo,sedgez_hi,&
              dx, dt, ubc, eps_for_bc, ppm_type)
      else
         call slopes( XVEL, &
                          lo,hi,&
                          u,u_lo,u_hi,&
                          sx,sx_lo,sx_hi,&
                          sy,sy_lo,sy_hi,&
                          sz,sz_lo,sz_hi,&
                          ubc)

      end if

      if (ppm_type .gt. 0) then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  ulo(i,j,k) = Ipx(i-1,j,k)
                  uhi(i,j,k) = Imx(i  ,j,k)
               end do
            end do
         end do
      else
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  ulo(i,j,k) = u(i-1,j,k) + (half  - dthx*u(i-1,j,k))*sx(i-1,j,k)
                  uhi(i,j,k) = u(i,  j,k) + (-half - dthx*u(i,  j,k))*sx(i,  j,k)
               end do
            end do
         end do
      end if

      if (use_minion .eq. 1 )then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  ulo(i,j,k) = ulo(i,j,k) + dth*tfx(i-1,j,k)
                  uhi(i,j,k) = uhi(i,j,k) + dth*tfx(i,  j,k)
               end do
            end do
         end do
      end if

      if (ppm_type .gt. 0) then
         call ppm(lo,hi,&
              v,v_lo,v_hi,&
              u,u_lo,u_hi,&
              v,v_lo,v_hi,&
              w,w_lo,w_hi,&
              Ipx,Ipx_lo,Ipx_hi,&
              Imx,Imx_lo,Imx_hi,&
              Ipy,Ipy_lo,Ipy_hi,&
              Imy,Imy_lo,Imy_hi,&
              Ipz,Ipz_lo,Ipz_hi,&
              Imz,Imz_lo,Imz_hi,&
              sm,sm_lo,sm_hi,&
              sp,sp_lo,sp_hi,&
              dsvl,dsvl_lo,dsvl_hi,&
              sedgex,sedgex_lo,sedgex_hi,&
              sedgey,sedgey_lo,sedgey_hi,&
              sedgez,sedgez_lo,sedgez_hi,&
              dx, dt, vbc, eps_for_bc, ppm_type)
      else
         call slopes( YVEL, &
                          lo,hi,&
                          v,v_lo,v_hi,&
                          sx,sx_lo,sx_hi,&
                          sy,sy_lo,sy_hi,&
                          sz,sz_lo,sz_hi,&
                          vbc)

      end if

      if (ppm_type .gt. 0) then
         do k = kmin-1,kmax+1
            do j = jmin,  jmax+1
               do i = imin-1,imax+1
                  vlo(i,j,k) = Ipy(i,j-1,k)
                  vhi(i,j,k) = Imy(i,j  ,k)
               end do
            end do
         end do
      else
         do k = kmin-1,kmax+1
            do j = jmin,  jmax+1
               do i = imin-1,imax+1
                  vlo(i,j,k) = v(i,j-1,k) + (half  - dthy*v(i,j-1,k))*sy(i,j-1,k)
                  vhi(i,j,k) = v(i,j,  k) + (-half - dthy*v(i,j,  k))*sy(i,j,  k)
               end do
            end do
         end do
      end if

      if (use_minion .eq. 1 )then
         do k = kmin-1,kmax+1
            do j = jmin,    jmax+1
               do i = imin-1,  imax+1
                  vlo(i,j,k) = vlo(i,j,k) + dth*tfy(i,j-1,k)
                  vhi(i,j,k) = vhi(i,j,k) + dth*tfy(i,j,  k)
               end do
            end do
         end do
      end if

      if (ppm_type .gt. 0) then
         call ppm(lo,hi,&
              w,w_lo,w_hi,&
              u,u_lo,u_hi,&
              v,v_lo,v_hi,&
              w,w_lo,w_hi,&
              Ipx,Ipx_lo,Ipx_hi,&
              Imx,Imx_lo,Imx_hi,&
              Ipy,Ipy_lo,Ipy_hi,&
              Imy,Imy_lo,Imy_hi,&
              Ipz,Ipz_lo,Ipz_hi,&
              Imz,Imz_lo,Imz_hi,&
              sm,sm_lo,sm_hi,&
              sp,sp_lo,sp_hi,&
              dsvl,dsvl_lo,dsvl_hi,&
              sedgex,sedgex_lo,sedgex_hi,&
              sedgey,sedgey_lo,sedgey_hi,&
              sedgez,sedgez_lo,sedgez_hi,&
              dx, dt, wbc, eps_for_bc, ppm_type)
      else
         call slopes( ZVEL, &
                          lo,hi,&
                          w,w_lo,w_hi,&
                          sx,sx_lo,sx_hi,&
                          sy,sy_lo,sy_hi,&
                          sz,sz_lo,sz_hi,&
                          wbc)

      end if

      if (ppm_type .gt. 0) then
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  wlo(i,j,k) = Ipz(i,j,k-1)
                  whi(i,j,k) = Imz(i,j,k  )
               end do
            end do
         end do
      else
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  wlo(i,j,k) = w(i,j,k-1) + (half  - dthz*w(i,j,k-1))*sz(i,j,k-1)
                  whi(i,j,k) = w(i,j,k  ) + (-half - dthz*w(i,j,k  ))*sz(i,j,k  )
               end do
            end do
         end do
      end if

      if (use_minion .eq. 1 )then
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,  imax+1
                  wlo(i,j,k) = wlo(i,j,k) + dth*tfz(i,j,k-1)
                  whi(i,j,k) = whi(i,j,k) + dth*tfz(i,j,k )
               end do
            end do
         end do
      end if


      call trans_xbc(lo,hi,&
           u,u_lo,u_hi,&
           ulo,ulo_lo,ulo_hi,&
           uhi,uhi_lo,uhi_hi,&
           ulo,ulo_lo,ulo_hi,&
           XVEL, ubc, eps_for_bc,.false.,.false.)

      call trans_ybc(lo,hi,&
           v,v_lo,v_hi,&
           vlo,vlo_lo,vlo_hi,&
           vhi,vhi_lo,vhi_hi,&
           vlo,vlo_lo,vlo_hi,&
           YVEL, vbc, eps_for_bc,.false.,.false.)

      call trans_zbc(lo,hi,&
           w,w_lo,w_hi,&
           wlo,wlo_lo,wlo_hi,&
           whi,whi_lo,whi_hi,&
           wlo,wlo_lo,wlo_hi,&
           ZVEL, wbc, eps_for_bc,.false.,.false.)

      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               tst = ulo(i,j,k)+uhi(i,j,k)
               val = merge(ulo(i,j,k),uhi(i,j,k),tst .ge. 0.0d0)
               ltm = &
                   ( (ulo(i,j,k) .le. 0.0d0) .and. &
                   (uhi(i,j,k) .ge. 0.0d0) ) .or. &
                   (abs(tst)   .lt. eps )
               ulo(i,j,k) = merge(0.0d0,val,ltm)
            end do
         end do
      end do

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               tst = vlo(i,j,k)+vhi(i,j,k)
               val = merge(vlo(i,j,k),vhi(i,j,k),tst .ge. 0.0d0)
               ltm = &
                   ( (vlo(i,j,k) .le. 0.0d0) .and.  &
                   (vhi(i,j,k) .ge. 0.0d0) ) .or. &
                   (abs(tst)   .lt. eps )
               vlo(i,j,k) = merge(0.0d0,val,ltm)
            end do
         end do
      end do
      
      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               tst = wlo(i,j,k)+whi(i,j,k)
               val = merge(wlo(i,j,k),whi(i,j,k),tst .ge. 0.0d0)
               ltm = &
                   ( (wlo(i,j,k) .le. 0.0d0) .and. &
                   (whi(i,j,k) .ge. 0.0d0) ) .or. &
                   (abs(tst)   .lt. eps )
               wlo(i,j,k) = merge(0.0d0,val,ltm)
            end do
         end do
      end do

      end subroutine transvel

    subroutine estate_premac(lo,hi,&
         s,s_lo,s_hi,&
         tf,tf_lo,tf_hi,&
         u,u_lo,u_hi,&
         xlo,xlo_lo,xlo_hi,&
         xhi,xhi_lo,xhi_hi,&
         sx,sx_lo,sx_hi,&
         uad,uad_lo,uad_hi,&
         xedge,xedge_lo,xedge_hi,&
         uedge,uedge_lo,uedge_hi,&
         xstate,xstate_lo,xstate_hi,&
         Imx,Imx_lo,Imx_hi,&
         Ipx,Ipx_lo,Ipx_hi,&
         sedgex,sedgex_lo,sedgex_hi,&
         v,v_lo,v_hi,&
         ylo,ylo_lo,ylo_hi,&
         yhi,yhi_lo,yhi_hi,&
         sy,sy_lo,sy_hi,&
         vad,vad_lo,vad_hi,&
         yedge,yedge_lo,yedge_hi,&
         vedge,vedge_lo,vedge_hi,&
         ystate,ystate_lo,ystate_hi,&
         Imy,Imy_lo,Imy_hi,&
         Ipy,Ipy_lo,Ipy_hi,&
         sedgey,sedgey_lo,sedgey_hi,&
         w,w_lo,w_hi,&
         zlo,zlo_lo,zlo_hi,&
         zhi,zhi_lo,zhi_hi,&
         sz,sz_lo,sz_hi,&
         wad,wad_lo,wad_hi,&
         zedge,zedge_lo,zedge_hi,&
         wedge,wedge_lo,wedge_hi,&
         zstate,zstate_lo,zstate_hi,&
         Imz,Imz_lo,Imz_hi,&
         Ipz,Ipz_lo,Ipz_hi,&
         sedgez,sedgez_lo,sedgez_hi,&
         dsvl,dsvl_lo,dsvl_hi,&
         sm,sm_lo,sm_hi,&
         sp,sp_lo,sp_hi,&
         xylo,xylo_lo,xylo_hi,&
         xzlo,xzlo_lo,xzlo_hi,&
         yxlo,yxlo_lo,yxlo_hi,&
         yzlo,yzlo_lo,yzlo_hi,&
         zxlo,zxlo_lo,zxlo_hi,&
         zylo,zylo_lo,zylo_hi,&
         xyhi,xyhi_lo,xyhi_hi,&
         xzhi,xzhi_lo,xzhi_hi,&
         yxhi,yxhi_lo,yxhi_hi,&
         yzhi,yzhi_lo,yzhi_hi,&
         zxhi,zxhi_lo,zxhi_hi,&
         zyhi,zyhi_lo,zyhi_hi,&
         corner_couple, &
         bc, dt, dx, n, nc, velpred, use_minion, ppm_type)
         
      ! Here is what this routine does
      ! 1. Trace values of state from cell-centers to faces using cell-centered velocities
      !     and excluding transverse corrections.  This produces data in xlo,xho.  Enforce
      !     boundary conditions on faces, then resolve the upwind value from xlo,xhi using
      !     uad and put this into xedge.  If uad is small, xedge=(xlo+xhi)/2
      ! 1a. Repeat for yedge and zedge
      ! 2A. If "new" corner coupling
      !    xylo = xlo - dt/3 * vad_avg * grady(state)
      !       where vad_avg = (vad_{j} + vad_{j+1})/2, and grady is gradient in y using yedge values
      !    xyhi = xhi - dt/3 * vad_avg * grady(state)
      !    Apply BCs
      !    Upwind with uad, answer into xylo
      !    Do the same with xzlo, xzhi based on wad, apply BCs, and upwind into xzlo
      !    Repeat the other way...that is yxlo = ylo - dt/3 * uad_avg * gradx(state), etc, and zx, zy
      !
      !    Now, compute stxlo,stxhi as traced states with transverse corrections computed using
      !       the mixed terms computed above, and the corresponding ad velocities
      !    
      !    Finally, result xstate is upwinded between stxlo,stxhi using UFACE, where
      !      if velpred!=1:
      !          UFACE = uedge
      !      else
      !          UFACE = stxlo + stxhi
      ! 2B. As in 2A, but no corner corrections. The transverse terms are simple extensions
      !       to 2D form, with one twist....
      !      vbar = (vad_{j} + vad_{j+1})/2
      !      if (vbar < 0)
      !         grady = (s_{j+1} - s_{j})/dy
      !      else
      !         grady = (s_{j} - s_{j-1})/dy
      !      Then corresponding transverse correction is - vbar*grady*dt/2
      !
      implicit none

      integer, intent(in) :: velpred, use_minion, ppm_type, bc(SDIM,2), n, nc
      real(rt), intent(in) :: dt, dx(SDIM)

      integer, dimension(3), intent(in) :: s_lo,s_hi,tf_lo,tf_hi,&
           u_lo,u_hi,xlo_lo,xlo_hi,xhi_lo,xhi_hi,sx_lo,sx_hi,uad_lo,uad_hi,&
           uedge_lo,uedge_hi,xedge_lo,xedge_hi,xstate_lo,xstate_hi,Imx_lo,Imx_hi,Ipx_lo,Ipx_hi,sedgex_lo,sedgex_hi,&
           v_lo,v_hi,ylo_lo,ylo_hi,yhi_lo,yhi_hi,sy_lo,sy_hi,vad_lo,vad_hi,&
           vedge_lo,vedge_hi,yedge_lo,yedge_hi,ystate_lo,ystate_hi,Imy_lo,Imy_hi,Ipy_lo,Ipy_hi,sedgey_lo,sedgey_hi,&
           w_lo,w_hi,zlo_lo,zlo_hi,zhi_lo,zhi_hi,sz_lo,sz_hi,wad_lo,wad_hi,&
           wedge_lo,wedge_hi,zedge_lo,zedge_hi,zstate_lo,zstate_hi,Imz_lo,Imz_hi,Ipz_lo,Ipz_hi,sedgez_lo,sedgez_hi,&
           dsvl_lo,dsvl_hi,sm_lo,sm_hi,sp_lo,sp_hi,lo,hi
           
      integer, dimension(3), intent(in) :: xylo_lo,xylo_hi,xzlo_lo,xzlo_hi,yxlo_lo,yxlo_hi, &
                                           yzlo_lo,yzlo_hi,zxlo_lo,zxlo_hi,zylo_lo,zylo_hi, &
                                           xyhi_lo,xyhi_hi,xzhi_lo,xzhi_hi,yxhi_lo,yxhi_hi, &
                                           yzhi_lo,yzhi_hi,zxhi_lo,zxhi_hi,zyhi_lo,zyhi_hi

      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc)
      real(rt), intent(in) :: tf(tf_lo(1):tf_hi(1),tf_lo(2):tf_hi(2),tf_lo(3):tf_hi(3),nc)
      real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
      real(rt), intent(inout) :: xlo(xlo_lo(1):xlo_hi(1),xlo_lo(2):xlo_hi(2),xlo_lo(3):xlo_hi(3))
      real(rt), intent(inout) :: xhi(xhi_lo(1):xhi_hi(1),xhi_lo(2):xhi_hi(2),xhi_lo(3):xhi_hi(3))
      real(rt), intent(inout) :: sx(sx_lo(1):sx_hi(1),sx_lo(2):sx_hi(2),sx_lo(3):sx_hi(3))
      real(rt), intent(in) :: uad(uad_lo(1):uad_hi(1),uad_lo(2):uad_hi(2),uad_lo(3):uad_hi(3))
      real(rt), intent(in) :: uedge(uedge_lo(1):uedge_hi(1),uedge_lo(2):uedge_hi(2),uedge_lo(3):uedge_hi(3))
      real(rt), intent(inout) :: xedge(xedge_lo(1):xedge_hi(1),xedge_lo(2):xedge_hi(2),xedge_lo(3):xedge_hi(3))
      real(rt), intent(inout) :: xstate(xstate_lo(1):xstate_hi(1),xstate_lo(2):xstate_hi(2),xstate_lo(3):xstate_hi(3),nc) ! result
      real(rt), intent(inout) :: Imx(Imx_lo(1):Imx_hi(1),Imx_lo(2):Imx_hi(2),Imx_lo(3):Imx_hi(3))
      real(rt), intent(inout) :: Ipx(Ipx_lo(1):Ipx_hi(1),Ipx_lo(2):Ipx_hi(2),Ipx_lo(3):Ipx_hi(3))
      real(rt), intent(inout) :: sedgex(sedgex_lo(1):sedgex_hi(1),sedgex_lo(2):sedgex_hi(2),sedgex_lo(3):sedgex_hi(3))
      real(rt), intent(in) :: v(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
      real(rt), intent(inout) :: ylo(ylo_lo(1):ylo_hi(1),ylo_lo(2):ylo_hi(2),ylo_lo(3):ylo_hi(3))
      real(rt), intent(inout) :: yhi(yhi_lo(1):yhi_hi(1),yhi_lo(2):yhi_hi(2),yhi_lo(3):yhi_hi(3))
      real(rt), intent(inout) :: sy(sy_lo(1):sy_hi(1),sy_lo(2):sy_hi(2),sy_lo(3):sy_hi(3))
      real(rt), intent(in) :: vad(vad_lo(1):vad_hi(1),vad_lo(2):vad_hi(2),vad_lo(3):vad_hi(3))
      real(rt), intent(in) :: vedge(vedge_lo(1):vedge_hi(1),vedge_lo(2):vedge_hi(2),vedge_lo(3):vedge_hi(3))
      real(rt), intent(inout) :: yedge(yedge_lo(1):yedge_hi(1),yedge_lo(2):yedge_hi(2),yedge_lo(3):yedge_hi(3))
      real(rt), intent(inout) :: ystate(ystate_lo(1):ystate_hi(1),ystate_lo(2):ystate_hi(2),ystate_lo(3):ystate_hi(3),nc) ! result
      real(rt), intent(inout) :: Imy(Imy_lo(1):Imy_hi(1),Imy_lo(2):Imy_hi(2),Imy_lo(3):Imy_hi(3))
      real(rt), intent(inout) :: Ipy(Ipy_lo(1):Ipy_hi(1),Ipy_lo(2):Ipy_hi(2),Ipy_lo(3):Ipy_hi(3))
      real(rt), intent(inout) :: sedgey(sedgey_lo(1):sedgey_hi(1),sedgey_lo(2):sedgey_hi(2),sedgey_lo(3):sedgey_hi(3))
      real(rt), intent(in) :: w(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
      real(rt), intent(inout) :: zlo(zlo_lo(1):zlo_hi(1),zlo_lo(2):zlo_hi(2),zlo_lo(3):zlo_hi(3))
      real(rt), intent(inout) :: zhi(zhi_lo(1):zhi_hi(1),zhi_lo(2):zhi_hi(2),zhi_lo(3):zhi_hi(3))
      real(rt), intent(inout) :: sz(sz_lo(1):sz_hi(1),sz_lo(2):sz_hi(2),sz_lo(3):sz_hi(3))
      real(rt), intent(in) :: wad(wad_lo(1):wad_hi(1),wad_lo(2):wad_hi(2),wad_lo(3):wad_hi(3))
      real(rt), intent(in) :: wedge(wedge_lo(1):wedge_hi(1),wedge_lo(2):wedge_hi(2),wedge_lo(3):wedge_hi(3))
      real(rt), intent(inout) :: zedge(zedge_lo(1):zedge_hi(1),zedge_lo(2):zedge_hi(2),zedge_lo(3):zedge_hi(3))
      real(rt), intent(inout) :: zstate(zstate_lo(1):zstate_hi(1),zstate_lo(2):zstate_hi(2),zstate_lo(3):zstate_hi(3),nc) ! result
      real(rt), intent(inout) :: Imz(Imz_lo(1):Imz_hi(1),Imz_lo(2):Imz_hi(2),Imz_lo(3):Imz_hi(3))
      real(rt), intent(inout) :: Ipz(Ipz_lo(1):Ipz_hi(1),Ipz_lo(2):Ipz_hi(2),Ipz_lo(3):Ipz_hi(3))
      real(rt), intent(inout) :: sedgez(sedgez_lo(1):sedgez_hi(1),sedgez_lo(2):sedgez_hi(2),sedgez_lo(3):sedgez_hi(3))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2),sp_lo(3):sp_hi(3))
      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2),dsvl_lo(3):dsvl_hi(3))
      
      real(rt), intent(inout) :: xylo(xylo_lo(1):xylo_hi(1),xylo_lo(2):xylo_hi(2),xylo_lo(3):xylo_hi(3))
      real(rt), intent(inout) :: xzlo(xzlo_lo(1):xzlo_hi(1),xzlo_lo(2):xzlo_hi(2),xzlo_lo(3):xzlo_hi(3))
      real(rt), intent(inout) :: yxlo(yxlo_lo(1):yxlo_hi(1),yxlo_lo(2):yxlo_hi(2),yxlo_lo(3):yxlo_hi(3))
      real(rt), intent(inout) :: yzlo(yzlo_lo(1):yzlo_hi(1),yzlo_lo(2):yzlo_hi(2),yzlo_lo(3):yzlo_hi(3))
      real(rt), intent(inout) :: zxlo(zxlo_lo(1):zxlo_hi(1),zxlo_lo(2):zxlo_hi(2),zxlo_lo(3):zxlo_hi(3))
      real(rt), intent(inout) :: zylo(zylo_lo(1):zylo_hi(1),zylo_lo(2):zylo_hi(2),zylo_lo(3):zylo_hi(3))
      
      real(rt), intent(inout) :: xyhi(xyhi_lo(1):xyhi_hi(1),xyhi_lo(2):xyhi_hi(2),xyhi_lo(3):xyhi_hi(3))
      real(rt), intent(inout) :: xzhi(xzhi_lo(1):xzhi_hi(1),xzhi_lo(2):xzhi_hi(2),xzhi_lo(3):xzhi_hi(3))
      real(rt), intent(inout) :: yxhi(yxhi_lo(1):yxhi_hi(1),yxhi_lo(2):yxhi_hi(2),yxhi_lo(3):yxhi_hi(3))
      real(rt), intent(inout) :: yzhi(yzhi_lo(1):yzhi_hi(1),yzhi_lo(2):yzhi_hi(2),yzhi_lo(3):yzhi_hi(3))
      real(rt), intent(inout) :: zxhi(zxhi_lo(1):zxhi_hi(1),zxhi_lo(2):zxhi_hi(2),zxhi_lo(3):zxhi_hi(3))
      real(rt), intent(inout) :: zyhi(zyhi_lo(1):zyhi_hi(1),zyhi_lo(2):zyhi_hi(2),zyhi_lo(3):zyhi_hi(3))
      
     
      real(rt) :: stxlo(lo(1)-2:hi(1)+2)
      real(rt) :: stxhi(lo(1)-2:hi(1)+2)
      real(rt) :: stylo(lo(2)-2:hi(2)+2)
      real(rt) :: styhi(lo(2)-2:hi(2)+2)
      real(rt) :: stzlo(lo(3)-2:hi(3)+2)
      real(rt) :: stzhi(lo(3)-2:hi(3)+2)
      real(rt) ::  hx, hy, hz, dth, dthx, dthy, dthz, ihx, ihy, ihz
      real(rt) ::  dt3, dt3x, dt3y, dt3z, dt4, dt4x, dt4y, dt4z
      real(rt) ::  dt6, dt6x, dt6y, dt6z
      real(rt) ::  tr1,tr2,ubar,vbar,wbar,stx,sty,stz,fu,fv,fw

      real(rt) ::  eps, eps_for_bc
      integer  :: i,j,k,L,imin,jmin,kmin,imax,jmax,kmax, inc, corner_couple
      logical  :: ltx,lty,ltz
      parameter( eps        = 1.0D-6 )
      parameter( eps_for_bc = 1.0D-10 )         

      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      dthz = half*dt / dx(3)
      dt3  = dt / 3.0d0
      dt3x = dt3 / dx(1)
      dt3y = dt3 / dx(2)
      dt3z = dt3 / dx(3)
      dt4  = dt / 4.0d0
      dt4x = dt4 / dx(1)
      dt4y = dt4 / dx(2)
      dt4z = dt4 / dx(3)
      dt6  = dt / 6.0d0
      dt6x = dt6 / dx(1)
      dt6y = dt6 / dx(2)
      dt6z = dt6 / dx(3)
      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      ihx  = 1.0d0/dx(1)
      ihy  = 1.0d0/dx(2)
      ihz  = 1.0d0/dx(3)
      
      
      do L=1,nc
      
!c
!c     compute the slopes
!c

      if (ppm_type .gt. 0) then
         call ppm(lo,hi,&
              s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
              u,u_lo,u_hi,&
              v,v_lo,v_hi,&
              w,w_lo,w_hi,&
              Ipx,Ipx_lo,Ipx_hi,&
              Imx,Imx_lo,Imx_hi,&
              Ipy,Ipy_lo,Ipy_hi,&
              Imy,Imy_lo,Imy_hi,&
              Ipz,Ipz_lo,Ipz_hi,&
              Imz,Imz_lo,Imz_hi,&
              sm,sm_lo,sm_hi,&
              sp,sp_lo,sp_hi,&
              dsvl,dsvl_lo,dsvl_hi,&
              sedgex,sedgex_lo,sedgex_hi,&
              sedgey,sedgey_lo,sedgey_hi,&
              sedgez,sedgez_lo,sedgez_hi,&
              dx, dt, bc, eps_for_bc, ppm_type)
      else
         call slopes( ALL, &
                          lo,hi,&
                          s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
                          sx,sx_lo,sx_hi,&
                          sy,sy_lo,sy_hi,&
                          sz,sz_lo,sz_hi,&
                          bc)

      end if

!c
!c     trace the state to the cell edges
!c
      if (ppm_type .gt. 0) then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  xlo(i,j,k) = Ipx(i-1,j,k)
                  xhi(i,j,k) = Imx(i  ,j,k)
               end do
            end do
         end do
      else
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  xlo(i,j,k) = s(i-1,j,k,L) + (half  - dthx*u(i-1,j,k))*sx(i-1,j,k)
                  xhi(i,j,k) = s(i,  j,k,L) + (-half - dthx*u(i,  j,k))*sx(i,  j,k)
               end do
            end do
         end do
      end if

      if (use_minion.eq.1)then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  xlo(i,j,k) = xlo(i,j,k) + dth*tf(i-1,j,k,L)
                  xhi(i,j,k) = xhi(i,j,k) + dth*tf(i,  j,k,L)
               end do
            end do
         end do
      end if

      
      call trans_xbc(lo,hi,&
           s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
           xlo,xlo_lo,xlo_hi,&
           xhi,xhi_lo,xhi_hi,&
           uad,uad_lo,uad_hi,&
           n+L-1, bc, eps_for_bc,.false.,.false.)
      

      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               fu  = merge(0.0d0,one,abs(uad(i,j,k)).lt.eps)
               stx = merge(xlo(i,j,k),xhi(i,j,k),uad(i,j,k) .ge. 0.0d0)
               xedge(i,j,k) = fu*stx + (one - fu)*half*(xhi(i,j,k)+xlo(i,j,k))
            end do
         end do
      end do

      if (ppm_type .gt. 0) then
         do k = kmin-1,kmax+1
            do j = jmin,  jmax+1
               do i = imin-1,imax+1
                  ylo(i,j,k) = Ipy(i,j-1,k)
                  yhi(i,j,k) = Imy(i,j  ,k)
               end do
            end do
         end do
      else
         do k = kmin-1,kmax+1
            do j = jmin,  jmax+1
               do i = imin-1,imax+1
                  ylo(i,j,k) = s(i,j-1,k,L) + (half  - dthy*v(i,j-1,k))*sy(i,j-1,k)
                  yhi(i,j,k) = s(i,j, k,L)  + (-half - dthy*v(i,j,  k))*sy(i,j, k)
               end do
            end do
         end do
      end if

      if (use_minion.eq.1)then
         do k = kmin-1,kmax+1
            do j = jmin, jmax+1
               do i = imin-1,  imax+1
                  ylo(i,j,k) = ylo(i,j,k) + dth*tf(i,j-1,k,L)
                  yhi(i,j,k) = yhi(i,j,k) + dth*tf(i,j,  k,L)
               end do
            end do
         end do
      end if

      call trans_ybc(lo,hi,&
              s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
              ylo,ylo_lo,ylo_hi,&
              yhi,yhi_lo,yhi_hi,&
              vad,vad_lo,vad_hi,&
              n+L-1, bc, eps_for_bc,.false.,.false.)
      
      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               fv  = merge(0.0d0,one,abs(vad(i,j,k)).lt.eps)
               sty = merge(ylo(i,j,k),yhi(i,j,k),vad(i,j,k) .ge. 0.0d0)
               yedge(i,j,k) = fv*sty + (one - fv)*half*(yhi(i,j,k)+ylo(i,j,k))
            end do
         end do
      end do

      if (ppm_type .gt. 0) then
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  zlo(i,j,k) = Ipz(i,j,k-1)
                  zhi(i,j,k) = Imz(i,j,k  )
               end do
            end do
         end do
      else
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  zlo(i,j,k) = s(i,j,k-1,L) + (half  - dthz*w(i,j,k-1))*sz(i,j,k-1)
                  zhi(i,j,k) = s(i,j,k  ,L) + (-half - dthz*w(i,j,k  ))*sz(i,j,k  )
               end do
            end do
         end do
      end if

      if (use_minion.eq.1)then
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,  imax+1
                  zlo(i,j,k) = zlo(i,j,k) + dth*tf(i,j,k-1,L)
                  zhi(i,j,k) = zhi(i,j,k) + dth*tf(i,j,k,L)
               end do
            end do
         end do
      end if

      call trans_zbc(lo,hi,&
              s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
              zlo,zlo_lo,zlo_hi,&
              zhi,zhi_lo,zhi_hi,&
              wad,wad_lo,wad_hi,&
              n+L-1, bc, eps_for_bc,.false.,.false.)

      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               fw  = merge(0.0d0,one,abs(wad(i,j,k)).lt.eps)
               stz = merge(zlo(i,j,k),zhi(i,j,k),wad(i,j,k) .ge. 0.0d0)
               zedge(i,j,k) = fw*stz + (one-fw)*half*(zhi(i,j,k)+zlo(i,j,k))
            end do
         end do
      end do

      if (corner_couple .ne. 0) then

!c
!c     NEW CORNER COUPLING CODE
!c

!c
!c     compute the corner-coupled terms:
!c     xylo/hi, xzlo/hi, yxlo/hi, yzlo/hi, zxlo/hi, zylo/hi
!c

!c     loop over appropriate xy faces
      do k=kmin-1,kmax+1
         do j=jmin,jmax
            do i=imin,imax+1
               xylo(i,j,k) = xlo(i,j,k) &
                   - dt6y*(vad(i-1,j+1,k)+vad(i-1,j,k)) &
                   *(yedge(i-1,j+1,k)-yedge(i-1,j,k))
               xyhi(i,j,k) = xhi(i,j,k) &
                   - dt6y*(vad(i  ,j+1,k)+vad(i  ,j,k)) &
                   *(yedge(i  ,j+1,k)-yedge(i  ,j,k))
            end do
         end do
      end do

!c     boundary conditions
     call trans_xbc(lo,hi,&
              s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
              xylo,xylo_lo,xylo_hi,&
              xyhi,xyhi_lo,xyhi_hi,&
              uad,uad_lo,uad_hi,&
              n+L-1, bc, eps_for_bc,.true.,.false.)

!c     upwind
      do k=kmin-1,kmax+1
         do j=jmin,jmax
            do i=imin,imax+1
               fu  = merge(0.0d0,one,abs(uad(i,j,k)).lt.eps)
               stx = merge(xylo(i,j,k),xyhi(i,j,k),uad(i,j,k) .ge. 0.0d0)
               xylo(i,j,k) = fu*stx + (one - fu)*half*(xyhi(i,j,k)+xylo(i,j,k))
            end do
         end do
      end do

!c     loop over appropriate xz faces
      do k=kmin,kmax
         do j=jmin-1,jmax+1
            do i=imin,imax+1
               xzlo(i,j,k) = xlo(i,j,k) &
                   - dt6z*(wad(i-1,j,k+1)+wad(i-1,j,k)) &
                   *(zedge(i-1,j,k+1)-zedge(i-1,j,k))
               xzhi(i,j,k) = xhi(i,j,k) &
                   - dt6z*(wad(i  ,j,k+1)+wad(i  ,j,k)) &
                   *(zedge(i  ,j,k+1)-zedge(i  ,j,k))
            end do
         end do
      end do

!c     boundary conditions
     call trans_xbc(lo,hi,&
              s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
              xzlo,xzlo_lo,xzlo_hi,&
              xzhi,xzhi_lo,xzhi_hi,&
              uad,uad_lo,uad_hi,&
              n+L-1, bc, eps_for_bc,.false.,.true.)
              
!c     upwind
      do k=kmin,kmax
         do j=jmin-1,jmax+1
            do i=imin,imax+1
               fu  = merge(0.0d0,one,abs(uad(i,j,k)).lt.eps)
               stx = merge(xzlo(i,j,k),xzhi(i,j,k),uad(i,j,k) .ge. 0.0d0)
               xzlo(i,j,k) = fu*stx + (one - fu)*half*(xzhi(i,j,k)+xzlo(i,j,k))
            end do
         end do
      end do

!c     loop over appropriate yx faces
      do k=kmin-1,kmax+1
         do j=jmin,jmax+1
            do i=imin,imax
               yxlo(i,j,k) = ylo(i,j,k) &
                   - dt6x*(uad(i+1,j-1,k)+uad(i,j-1,k)) &
                   *(xedge(i+1,j-1,k)-xedge(i,j-1,k))
               yxhi(i,j,k) = yhi(i,j,k) &
                   - dt6x*(uad(i+1,j  ,k)+uad(i,j  ,k)) &
                   *(xedge(i+1,j  ,k)-xedge(i,j  ,k))
            end do
         end do
      end do

!c     boundary conditions
     call trans_ybc(lo,hi,&
              s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
              yxlo,yxlo_lo,yxlo_hi,&
              yxhi,yxhi_lo,yxhi_hi,&
              vad,vad_lo,vad_hi,&
              n+L-1, bc, eps_for_bc,.true.,.false.)
              
!c     upwind
      do k=kmin-1,kmax+1
         do j=jmin,jmax+1
            do i=imin,imax
               fv  = merge(0.0d0,one,abs(vad(i,j,k)).lt.eps)
               sty = merge(yxlo(i,j,k),yxhi(i,j,k),vad(i,j,k) .ge. 0.0d0)
               yxlo(i,j,k) = fv*sty + (one - fv)*half*(yxhi(i,j,k)+yxlo(i,j,k))
            end do
         end do
      end do

!c     loop over appropriate yz faces
      do k=kmin,kmax
         do j=jmin,jmax+1
            do i=imin-1,imax+1
               yzlo(i,j,k) = ylo(i,j,k) &
                   - dt6z*(wad(i,j-1,k+1)+wad(i,j-1,k)) &
                   *(zedge(i,j-1,k+1)-zedge(i,j-1,k))
               yzhi(i,j,k) = yhi(i,j,k) &
                   - dt6z*(wad(i,j  ,k+1)+wad(i,j  ,k)) &
                   *(zedge(i,j  ,k+1)-zedge(i,j  ,k))
            end do
         end do
      end do

!c     boundary conditions
     call trans_ybc(lo,hi,&
              s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
              yzlo,yzlo_lo,yzlo_hi,&
              yzhi,yzhi_lo,yzhi_hi,&
              vad,vad_lo,vad_hi,&
              n+L-1, bc, eps_for_bc,.false.,.true.)
              
!c     upwind
      do k=kmin,kmax
         do j=jmin,jmax+1
            do i=imin-1,imax+1
               fv  = merge(0.0d0,one,abs(vad(i,j,k)).lt.eps)
               sty = merge(yzlo(i,j,k),yzhi(i,j,k),vad(i,j,k) .ge. 0.0d0)
               yzlo(i,j,k) = fv*sty + (one - fv)*half*(yzhi(i,j,k)+yzlo(i,j,k))
            end do
         end do
      end do

!c     loop over appropriate zx faces
      do k=kmin,kmax+1
         do j=jmin-1,jmax+1
            do i=imin,imax
               zxlo(i,j,k) = zlo(i,j,k) &
                   - dt6x*(uad(i+1,j,k-1)+uad(i,j,k-1)) &
                   *(xedge(i+1,j,k-1)-xedge(i,j,k-1))
               zxhi(i,j,k) = zhi(i,j,k) &
                   - dt6x*(uad(i+1,j,k  )+uad(i,j,k  )) &
                   *(xedge(i+1,j,k  )-xedge(i,j,k  ))
            end do
         end do
      end do

!c     boundary conditions
     call trans_zbc(lo,hi,&
              s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
              zxlo,zxlo_lo,zxlo_hi,&
              zxhi,zxhi_lo,zxhi_hi,&
              wad,wad_lo,wad_hi,&
              n+L-1, bc, eps_for_bc,.true.,.false.)

!c     upwind
      do k=kmin,kmax+1
         do j=jmin-1,jmax+1
            do i=imin,imax
               fw  = merge(0.0d0,one,abs(wad(i,j,k)).lt.eps)
               stz = merge(zxlo(i,j,k),zxhi(i,j,k),wad(i,j,k) .ge. 0.0d0)
               zxlo(i,j,k) = fw*stz + (one-fw)*half*(zxhi(i,j,k)+zxlo(i,j,k))
            end do
         end do
      end do

!c     loop over appropriate zy faces
      do k=kmin,kmax+1
         do j=jmin,jmax
            do i=imin-1,imax+1
               zylo(i,j,k) = zlo(i,j,k) &
                   - dt6y*(vad(i,j+1,k-1)+vad(i,j,k-1)) &
                   *(yedge(i,j+1,k-1)-yedge(i,j,k-1))
               zyhi(i,j,k) = zhi(i,j,k) &
                   - dt6y*(vad(i,j+1,k  )+vad(i,j,k  ))  &
                   *(yedge(i,j+1,k  )-yedge(i,j,k  ))
            end do
         end do
      end do

!c     boundary conditions
     call trans_zbc(lo,hi,&
              s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
              zylo,zylo_lo,zylo_hi,&
              zyhi,zyhi_lo,zyhi_hi,&
              wad,wad_lo,wad_hi,&
              n+L-1, bc, eps_for_bc,.false.,.true.)
              
!c     upwind
      do k=kmin,kmax+1
         do j=jmin,jmax
            do i=imin-1,imax+1
               fw  = merge(0.0d0,one,abs(wad(i,j,k)).lt.eps)
               stz = merge(zylo(i,j,k),zyhi(i,j,k),wad(i,j,k) .ge. 0.0d0)
               zylo(i,j,k) = fw*stz + (one-fw)*half*(zyhi(i,j,k)+zylo(i,j,k))
            end do
         end do
      end do

!c
!c     compute the xedge states
!c
      if ((velpred.ne.1) .or. ((n+L-1).eq.XVEL)) then
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax+1

                  stxlo(i) = xlo(i,j,k) &
                      - dt4y*(vad(i-1,j+1,k  )+vad(i-1,j,k))* &
                      (yzlo(i-1,j+1,k  )-yzlo(i-1,j,k)) &
                      - dt4z*(wad(i-1,j  ,k+1)+wad(i-1,j,k))* &
                      (zylo(i-1,j  ,k+1)-zylo(i-1,j,k))
                  stxhi(i) = xhi(i,j,k) &
                      - (dt4/hy)*(vad(i  ,j+1,k  )+vad(i  ,j,k))* &
                      (yzlo(i  ,j+1,k  )-yzlo(i  ,j,k)) &
                      - (dt4/hz)*(wad(i  ,j  ,k+1)+wad(i  ,j,k))* &
                      (zylo(i  ,j  ,k+1)-zylo(i  ,j,k))

                  if (use_minion.eq.0) then
                     stxlo(i) = stxlo(i) + dth*tf(i-1,j,k,L)
                     stxhi(i) = stxhi(i) + dth*tf(i,  j,k,L)
                  end if

               end do

               if (bc(1,1).eq.EXT_DIR .and. velpred.eq.1) then
                  stxhi(imin) = s(imin-1,j,k,L)
                  stxlo(imin) = s(imin-1,j,k,L)
               else if (bc(1,1).eq.EXT_DIR .and. uad(imin,j,k).ge.0.0d0) then
                  stxhi(imin) = s(imin-1,j,k,L)
                  stxlo(imin) = s(imin-1,j,k,L)
               else if (bc(1,1).eq.EXT_DIR .and. uad(imin,j,k).lt.0.0d0) then
                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP) then
                  if ((n+L-1).eq.XVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                        stxhi(imin) = MIN(stxhi(imin),0.0d0)
#endif
                        stxlo(imin) = stxhi(imin)
                     else
                        if (uad(imin,j,k).ge.0.0d0) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                           stxhi(imin) = MIN(stxhi(imin),0.0d0)
#endif
                           stxlo(imin) = stxhi(imin)
                        endif
                     endif
                  else
                     stxlo(imin) = stxhi(imin)
#ifdef NOVERTICALINFLOW
!c     Hack for no vertical inflow velocity
                     if ((n+L-1).eq.ZVEL)then
                        stxlo(imin) = 0.0d0
                     endif
#endif
                  endif
               else if (bc(1,1).eq.REFLECT_EVEN) then
                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1).eq.REFLECT_ODD) then
                  stxhi(imin) = 0.0d0
                  stxlo(imin) = 0.0d0
               end if
               if (bc(1,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stxlo(imax+1) = s(imax+1,j,k,L)
                  stxhi(imax+1) = s(imax+1,j,k,L)
               else if (bc(1,2).eq.EXT_DIR .and. uad(imax+1,j,k).le.0.0d0) then
                  stxlo(imax+1) = s(imax+1,j,k,L)
                  stxhi(imax+1) = s(imax+1,j,k,L)
               else if (bc(1,2).eq.EXT_DIR .and. uad(imax+1,j,k).gt.0.0d0) then
                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP) then
                  if ((n+L-1).eq.XVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                        stxlo(imax+1) = MAX(stxlo(imax+1),0.0d0)
#endif
                        stxhi(imax+1) = stxlo(imax+1)
                     else
                        if (uad(imax+1,j,k).le.0.0d0) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                           stxlo(imax+1) = MAX(stxlo(imax+1),0.0d0)
#endif
                           stxhi(imax+1) = stxlo(imax+1)
                        endif
                     endif
                  else
                     stxhi(imax+1) = stxlo(imax+1)
#ifdef NOVERTICALINFLOW
!c     Hack for no vertical inflow velocity
                     if ((n+L-1).eq.ZVEL)then
                        stxhi(imax+1) = 0.0d0
                     endif
#endif
                  endif
               else if (bc(1,2).eq.REFLECT_EVEN) then
                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2).eq.REFLECT_ODD) then
                  stxlo(imax+1) = zero
                  stxhi(imax+1) = zero
               end if

               if ( velpred .eq. 1 ) then
                  do i = imin, imax+1
                     ltx = stxlo(i) .le. 0.0d0  .and.  stxhi(i) .ge. 0.0d0
                     ltx = ltx .or. (abs(stxlo(i)+stxhi(i)) .lt. eps)
                     stx = merge(stxlo(i),stxhi(i),(stxlo(i)+stxhi(i)) .ge. 0.0d0)
                     xstate(i,j,k,L) = merge(zero,stx,ltx)
                  end do
               else
                  do i = imin, imax+1
                     xstate(i,j,k,L) = merge(stxlo(i),stxhi(i),uedge(i,j,k) .ge. 0.0d0)
                     xstate(i,j,k,L) = merge(half*(stxlo(i)+stxhi(i)),xstate(i,j,k,L) &
                         ,abs(uedge(i,j,k)).lt.eps)
                  end do
               end if
            end do
         end do
      end if
!c
!c     compute the yedge states
!c
      if ((velpred.ne.1) .or. ((n+L-1).eq.YVEL)) then
         do k = kmin,kmax
            do i = imin,imax
               do j = jmin,jmax+1

                  stylo(j) = ylo(i,j,k) &
                      - dt4x*(uad(i+1,j-1,k  )+uad(i,j-1,k))*  &
                      (xzlo(i+1,j-1,k  )-xzlo(i,j-1,k)) &
                      - dt4z*(wad(i  ,j-1,k+1)+wad(i,j-1,k))* &
                      (zxlo(i  ,j-1,k+1)-zxlo(i,j-1,k))
                  styhi(j) = yhi(i,j,k) &
                      - dt4x*(uad(i+1,j  ,k  )+uad(i,j  ,k))* &
                      (xzlo(i+1,j  ,k  )-xzlo(i,j  ,k)) &
                      - dt4z*(wad(i  ,j  ,k+1)+wad(i,j  ,k))* &
                      (zxlo(i  ,j  ,k+1)-zxlo(i,j  ,k))

                  if (use_minion.eq.0) then
                     stylo(j) = stylo(j) + dth*tf(i,j-1,k,L)
                     styhi(j) = styhi(j) + dth*tf(i,j,  k,L)
                  end if

               end do

               if (bc(2,1).eq.EXT_DIR .and. velpred.eq.1) then
                  styhi(jmin) = s(i,jmin-1,k,L)
                  stylo(jmin) = s(i,jmin-1,k,L)
               else if (bc(2,1).eq.EXT_DIR .and. vad(i,jmin,k).ge.0.0d0) then
                  styhi(jmin) = s(i,jmin-1,k,L)
                  stylo(jmin) = s(i,jmin-1,k,L)
               else if (bc(2,1).eq.EXT_DIR .and. vad(i,jmin,k).lt.0.0d0) then
                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP) then
                  if ((n+L-1).eq.YVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                        styhi(jmin) = MIN(styhi(jmin),0.0d0)
#endif
                        stylo(jmin) = styhi(jmin)
                     else
                        if (vad(i,jmin,k).ge.0.0d0) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                           styhi(jmin) = MIN(styhi(jmin),0.0d0)
#endif
                           stylo(jmin) = styhi(jmin)
                        endif
                     endif
                  else
                     stylo(jmin) = styhi(jmin)
#ifdef NOVERTICALINFLOW
!c     Hack for no vertical inflow velocity
                     if ((n+L-1).eq.ZVEL)then
                        stylo(jmin) = 0.0d0
                     endif
#endif
                  endif
               else if (bc(2,1).eq.REFLECT_EVEN) then
                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1).eq.REFLECT_ODD) then
                  styhi(jmin) = 0.0d0
                  stylo(jmin) = 0.0d0
               end if
               
               if (bc(2,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stylo(jmax+1) = s(i,jmax+1,k,L)
                  styhi(jmax+1) = s(i,jmax+1,k,L)
               else if (bc(2,2).eq.EXT_DIR .and. vad(i,jmax+1,k).le.0.0d0) then
                  stylo(jmax+1) = s(i,jmax+1,k,L)
                  styhi(jmax+1) = s(i,jmax+1,k,L)
               else if (bc(2,2).eq.EXT_DIR .and. vad(i,jmax+1,k).gt.0.0d0) then
                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP) then
                  if ((n+L-1).eq.YVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                        stylo(jmax+1) = MAX(stylo(jmax+1),0.0d0)
#endif
                        styhi(jmax+1) = stylo(jmax+1)
                     else
                        if (vad(i,jmax+1,k).le.0.0d0) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                           stylo(jmax+1) = MAX(stylo(jmax+1),0.0d0)
#endif
                           styhi(jmax+1) = stylo(jmax+1)
                        endif
                     endif
                  else
                     styhi(jmax+1) = stylo(jmax+1)
#ifdef NOVERTICALINFLOW
!c     Hack for no vertical inflow velocity
                     if ((n+L-1).eq.ZVEL)then
                        styhi(jmax+1) = 0.0d0
                     endif
#endif
                  endif
               else if (bc(2,2).eq.REFLECT_EVEN) then
                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2).eq.REFLECT_ODD) then
                  stylo(jmax+1) = 0.0d0
                  styhi(jmax+1) = 0.0d0
               end if

               if ( velpred .eq. 1 ) then
                  do j = jmin, jmax+1
                     lty = stylo(j) .le. zero  .and.  styhi(j) .ge. 0.0d0
                     lty = lty .or. (abs(stylo(j)+styhi(j)) .lt. eps)
                     sty = merge(stylo(j),styhi(j),(stylo(j)+styhi(j)) .ge. 0.0d0)
                     ystate(i,j,k,L) = merge(0.0d0,sty,lty)
                  end do
               else
                  do j=jmin,jmax+1
                     ystate(i,j,k,L) = merge(stylo(j),styhi(j),vedge(i,j,k) .ge. 0.0d0)
                     ystate(i,j,k,L) = merge(half*(stylo(j)+styhi(j)),ystate(i,j,k,L), &
                         abs(vedge(i,j,k)).lt.eps)
                  end do
               end if
            end do
         end do
      end if
!c
!c     compute the zedge states
!c
      if ((velpred.ne.1) .or. ((n+L-1).eq.ZVEL)) then
         do j = jmin,jmax
            do i = imin,imax
               do k = kmin,kmax+1

                  stzlo(k) = zlo(i,j,k) &
                      - dt4x*(uad(i+1,j  ,k-1)+uad(i,j,k-1)) &
                      *(xylo(i+1,j  ,k-1)-xylo(i,j,k-1)) &
                      - dt4y*(vad(i  ,j+1,k-1)+vad(i,j,k-1)) &
                      *(yxlo(i  ,j+1,k-1)-yxlo(i,j,k-1))
                  
                  stzhi(k) = zhi(i,j,k) &
                      - dt4x*(uad(i+1,j  ,k  )+uad(i,j,k  )) &
                      *(xylo(i+1,j  ,k  )-xylo(i,j,k  )) &
                      - dt4y*(vad(i  ,j+1,k  )+vad(i,j,k  )) &
                      *(yxlo(i  ,j+1,k  )-yxlo(i,j,k  ))
                  
                  if (use_minion.eq.0) then
                     stzlo(k) = stzlo(k) + dth*tf(i,j,k-1,L)
                     stzhi(k) = stzhi(k) + dth*tf(i,j,k,L)
                  end if

               end do

               if (bc(3,1).eq.EXT_DIR .and. velpred.eq.1) then
                  stzlo(kmin) = s(i,j,kmin-1,L)
                  stzhi(kmin) = s(i,j,kmin-1,L)
               else if (bc(3,1).eq.EXT_DIR .and. wad(i,j,kmin).ge.0.0d0) then
                  stzlo(kmin) = s(i,j,kmin-1,L)
                  stzhi(kmin) = s(i,j,kmin-1,L)
               else if (bc(3,1).eq.EXT_DIR .and. wad(i,j,kmin).lt.0.0d0) then
                  stzlo(kmin) = stzhi(kmin)
               else if (bc(3,1).eq.FOEXTRAP.or.bc(3,1).eq.HOEXTRAP) then
                  if ((n+L-1).eq.ZVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWZINFLOW
!c     prevent backflow
                        stzhi(kmin) = MIN(stzhi(kmin),0.0d0)
#endif
                        stzlo(kmin) = stzhi(kmin)
                     else
                        if (wad(i,j,kmin).ge.0.0d0) then
#ifndef ALLOWZINFLOW
!c     prevent backflow
                           stzhi(kmin) = MIN(stzhi(kmin),0.0d0)
#endif
                           stzlo(kmin) = stzhi(kmin)
                        endif
                     endif
                  else
                     stzlo(kmin) = stzhi(kmin)
                  endif
               else if (bc(3,1).eq.REFLECT_EVEN) then
                  stzlo(kmin) = stzhi(kmin)
               else if (bc(3,1).eq.REFLECT_ODD) then
                  stzlo(kmin) = 0.0d0
                  stzhi(kmin) = 0.0d0
               end if
               if (bc(3,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stzlo(kmax+1) = s(i,j,kmax+1,L)
                  stzhi(kmax+1) = s(i,j,kmax+1,L)
               else if (bc(3,2).eq.EXT_DIR .and. wad(i,j,kmax+1).le.0.0d0) then
                  stzlo(kmax+1) = s(i,j,kmax+1,L)
                  stzhi(kmax+1) = s(i,j,kmax+1,L)
               else if (bc(3,2).eq.EXT_DIR .and. wad(i,j,kmax+1).gt.0.0d0) then
                  stzhi(kmax+1) = stzlo(kmax+1)
               else if (bc(3,2).eq.FOEXTRAP.or.bc(3,2).eq.HOEXTRAP) then
                  if ((n+L-1).eq.ZVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWZINFLOW
!c     prevent backflow
                        stzlo(kmax+1) = MAX(stzlo(kmax+1),0.0d0)
#endif
                        stzhi(kmax+1) = stzlo(kmax+1)
                     else
                        if (wad(i,j,kmax+1).le.zero) then
#ifndef ALLOWZINFLOW
!c     prevent backflow
                           stzlo(kmax+1) = MAX(stzlo(kmax+1),0.0d0)
#endif
                           stzhi(kmax+1) = stzlo(kmax+1)
                        endif
                     endif
                  else
                     stzhi(kmax+1) = stzlo(kmax+1)
                  endif
               else if (bc(3,2).eq.REFLECT_EVEN) then
                  stzhi(kmax+1) = stzlo(kmax+1)
               else if (bc(3,2).eq.REFLECT_ODD) then
                  stzlo(kmax+1) = 0.0d0
                  stzhi(kmax+1) = 0.0d0
               end if

               if ( velpred .eq. 1 ) then
                  do k = kmin,kmax+1
                     ltz = stzlo(k) .le. zero  .and.  stzhi(k) .ge. 0.0d0
                     ltz = ltz .or. (abs(stzlo(k)+stzhi(k)) .lt. eps)
                     stz = merge(stzlo(k),stzhi(k),(stzlo(k)+stzhi(k)) .ge. 0.0d0)
                     zstate(i,j,k,L) = merge(0.0d0,stz,ltz)
                  end do
               else
                  do k = kmin,kmax+1
                     zstate(i,j,k,L) = merge(stzlo(k),stzhi(k),wedge(i,j,k) .ge. 0.0d0)
                     zstate(i,j,k,L) = merge(half*(stzlo(k)+stzhi(k)),zstate(i,j,k,L), &
                         abs(wedge(i,j,k)).lt.eps)
                  end do
               end if
            end do
         end do
      end if

      else

!c    
!c     ORIGINAL NON-CORNER COUPLING CODE
!c
!c
!c     compute the xedge states
!c
      if ((velpred.ne.1) .or. ((n+L-1).eq.XVEL)) then

         do k = kmin,kmax
            do j = jmin,jmax

               do i = imin-1,imax+1
                  if (vad(i,j,k)*vad(i,j+1,k).lt.0.d0) then
                      vbar = 0.5d0*(vad(i,j,k)+vad(i,j+1,k))
                      if (vbar.lt.0.d0) then
                          inc = 1
                      else
                          inc = 0
                      endif
                      tr1 = vbar*(s(i,j+inc,k,L)-s(i,j+inc-1,k,L))*ihy
                  else
                      tr1 = half* &
                      (vad(i,j+1,k) + vad(i,j,k)) * &
                      (yedge(i,j+1,k) - yedge(i,j,k)) *ihy
                  endif

                  if (wad(i,j,k)*wad(i,j,k+1).lt.0.d0) then
                      wbar = 0.5d0*(wad(i,j,k)+wad(i,j,k+1))
                      if (wbar.lt.0.d0) then
                          inc = 1
                      else
                          inc = 0
                      endif
                      tr2 = wbar*(s(i,j,k+inc,L)-s(i,j,k+inc-1,L))*ihz
                  else
                      tr2 = half* &
                      (wad(i,j,k+1) + wad(i,j,k)) * &
                      (zedge(i,j,k+1) - zedge(i,j,k)) *ihz
                  endif

                  if (ppm_type .gt. 0) then
                     stxlo(i+1)= Ipx(i,j,k) &
                         - dth*tr1 - dth*tr2 &
                         + dth*tf(i,j,k,L)
                     stxhi(i  )= Imx(i,j,k) &
                         - dth*tr1 - dth*tr2 &
                         + dth*tf(i,j,k,L)
                  else
                     stxlo(i+1)= s(i,j,k,L) + (half-dthx*u(i,j,k))*sx(i,j,k) &
                         - dth*tr1 - dth*tr2 &
                         + dth*tf(i,j,k,L)
                     stxhi(i  )= s(i,j,k,L) - (half+dthx*u(i,j,k))*sx(i,j,k) &
                         - dth*tr1 - dth*tr2 &
                         + dth*tf(i,j,k,L)
                  end if
               end do

               if (bc(1,1).eq.EXT_DIR .and. velpred.eq.1) then
                  stxhi(imin) = s(imin-1,j,k,L)
                  stxlo(imin) = s(imin-1,j,k,L)
               else if (bc(1,1).eq.EXT_DIR .and. uad(imin,j,k).ge.0.0d0) then
                  stxhi(imin) = s(imin-1,j,k,L)
                  stxlo(imin) = s(imin-1,j,k,L)
               else if (bc(1,1).eq.EXT_DIR .and. uad(imin,j,k).lt.0.0d0) then
                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP) then
                  if ((n+L-1).eq.XVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                        stxhi(imin) = MIN(stxhi(imin),0.0d0)
#endif
                        stxlo(imin) = stxhi(imin)
                     else
                        if (uad(imin,j,k).ge.zero) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                           stxhi(imin) = MIN(stxhi(imin),0.0d0)
#endif
                           stxlo(imin) = stxhi(imin)
                        endif
                     endif
                  else
                     stxlo(imin) = stxhi(imin)
#ifdef NOVERTICALINFLOW
!c     Hack for no vertical inflow velocity
                     if ((n+L-1).eq.ZVEL)then
                        stxlo(imin) = 0.0d0
                     endif
#endif
                  endif
               else if (bc(1,1).eq.REFLECT_EVEN) then
                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1).eq.REFLECT_ODD) then
                  stxhi(imin) = 0.0d0
                  stxlo(imin) = 0.0d0
               end if
               if (bc(1,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stxlo(imax+1) = s(imax+1,j,k,L)
                  stxhi(imax+1) = s(imax+1,j,k,L)
               else if (bc(1,2).eq.EXT_DIR .and. uad(imax+1,j,k).le.0.0d0) then
                  stxlo(imax+1) = s(imax+1,j,k,L)
                  stxhi(imax+1) = s(imax+1,j,k,L)
               else if (bc(1,2).eq.EXT_DIR .and. uad(imax+1,j,k).gt.0.0d0) then
                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP) then
                  if ((n+L-1).eq.XVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                        stxlo(imax+1) = MAX(stxlo(imax+1),0.0d0)
#endif
                        stxhi(imax+1) = stxlo(imax+1)
                     else
                        if (uad(imax+1,j,k).le.zero) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                           stxlo(imax+1) = MAX(stxlo(imax+1),0.0d0)
#endif
                           stxhi(imax+1) = stxlo(imax+1)
                        endif
                     endif
                  else
                     stxhi(imax+1) = stxlo(imax+1)
#ifdef NOVERTICALINFLOW
!c     Hack for no vertical inflow velocity
                     if ((n+L-1).eq.ZVEL)then
                        stxhi(imax+1) = 0.0d0
                     endif
#endif
                  endif
               else if (bc(1,2).eq.REFLECT_EVEN) then
                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2).eq.REFLECT_ODD) then
                  stxlo(imax+1) = 0.0d0
                  stxhi(imax+1) = 0.0d0
               end if

               if ( velpred .eq. 1 ) then
                  do i = imin, imax+1
                     ltx = stxlo(i) .le. 0.0d0  .and.  stxhi(i) .ge. 0.0d0
                     ltx = ltx .or. (abs(stxlo(i)+stxhi(i)) .lt. eps)
                     stx = merge(stxlo(i),stxhi(i),(stxlo(i)+stxhi(i)) .ge. 0.0d0)
                     xstate(i,j,k,L) = merge(0.0d0,stx,ltx)
                  end do
               else
                  do i = imin, imax+1
                     xstate(i,j,k,L) = merge(stxlo(i),stxhi(i),uedge(i,j,k) .ge. 0.0d0)
                     xstate(i,j,k,L) = merge(half*(stxlo(i)+stxhi(i)),xstate(i,j,k,L) &
                         ,abs(uedge(i,j,k)).lt.eps)
                  end do
               end if
            end do
         end do
      end if
!c
!c     compute the yedge states
!c
      if ((velpred.ne.1) .or. ((n+L-1).eq.YVEL)) then
         do k = kmin,kmax
            do i = imin,imax
               
               do j = jmin-1,jmax+1
                  if (uad(i,j,k)*uad(i+1,j,k).lt.0.d0) then
                      ubar = 0.5d0*(uad(i,j,k)+uad(i+1,j,k))
                      if (ubar.lt.0.d0) then
                          inc = 1
                      else
                          inc = 0
                      endif
                      tr1 = ubar*(s(i+inc,j,k,L)-s(i+inc-1,j,k,L))*ihx
                  else
                      tr1 = half* &
                      (uad(i+1,j,k) + uad(i,j,k)) * &
                      (xedge(i+1,j,k) - xedge(i,j,k)) *ihx
                  endif

                  if (wad(i,j,k)*wad(i,j,k+1).lt.0.d0) then
                      wbar = 0.5d0*(wad(i,j,k)+wad(i,j,k+1))
                      if (wbar.lt.0.d0) then
                          inc = 1
                      else
                          inc = 0
                      endif
                      tr2 = wbar*(s(i,j,k+inc,L)-s(i,j,k+inc-1,L))*ihz
                  else
                      tr2 = half* &
                      (wad(i,j,k+1) + wad(i,j,k)) * &
                      (zedge(i,j,k+1) - zedge(i,j,k)) *ihz
                  endif

                  if (ppm_type .gt. 0) then
                     stylo(j+1)= Ipy(i,j,k) &
                         - dth*tr1 - dth*tr2 &
                         + dth*tf(i,j,k,L)
                     styhi(j)  = Imy(i,j,k) &
                         - dth*tr1 - dth*tr2 &
                         + dth*tf(i,j,k,L)
                  else
                     stylo(j+1)= s(i,j,k,L) + (half-dthy*v(i,j,k))*sy(i,j,k) &
                         - dth*tr1 - dth*tr2 &
                         + dth*tf(i,j,k,L)
                     styhi(j)  = s(i,j,k,L) - (half+dthy*v(i,j,k))*sy(i,j,k) &
                         - dth*tr1 - dth*tr2 &
                         + dth*tf(i,j,k,L)
                  end if
               end do

               if (bc(2,1).eq.EXT_DIR .and. velpred.eq.1) then
                  styhi(jmin) = s(i,jmin-1,k,L)
                  stylo(jmin) = s(i,jmin-1,k,L)
               else if (bc(2,1).eq.EXT_DIR .and. vad(i,jmin,k).ge.0.0d0) then
                  styhi(jmin) = s(i,jmin-1,k,L)
                  stylo(jmin) = s(i,jmin-1,k,L)
               else if (bc(2,1).eq.EXT_DIR .and. vad(i,jmin,k).lt.0.0d0) then
                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP) then
                  if ((n+L-1).eq.YVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                        styhi(jmin) = MIN(styhi(jmin),0.0d0)
#endif
                        stylo(jmin) = styhi(jmin)
                     else
                        if (vad(i,jmin,k).ge.zero) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                           styhi(jmin) = MIN(styhi(jmin),0.0d0)
#endif
                           stylo(jmin) = styhi(jmin)
                        endif
                     endif
                  else
                     stylo(jmin) = styhi(jmin)
#ifdef NOVERTICALINFLOW
!c     Hack for no vertical inflow velocity
                     if ((n+L-1).eq.ZVEL)then
                        stylo(jmin) = 0.0d0
                     endif
#endif
                  endif
               else if (bc(2,1).eq.REFLECT_EVEN) then
                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1).eq.REFLECT_ODD) then
                  styhi(jmin) = 0.0d0
                  stylo(jmin) = 0.0d0
               end if
               
               if (bc(2,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stylo(jmax+1) = s(i,jmax+1,k,L)
                  styhi(jmax+1) = s(i,jmax+1,k,L)
               else if (bc(2,2).eq.EXT_DIR .and. vad(i,jmax+1,k).le.0.0d0) then
                  stylo(jmax+1) = s(i,jmax+1,k,L)
                  styhi(jmax+1) = s(i,jmax+1,k,L)
               else if (bc(2,2).eq.EXT_DIR .and. vad(i,jmax+1,k).gt.0.0d0) then
                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP) then
                  if ((n+L-1).eq.YVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                        stylo(jmax+1) = MAX(stylo(jmax+1),0.0d0)
#endif
                        styhi(jmax+1) = stylo(jmax+1)
                     else
                        if (vad(i,jmax+1,k).le.zero) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                           stylo(jmax+1) = MAX(stylo(jmax+1),0.0d0)
#endif
                           styhi(jmax+1) = stylo(jmax+1)
                        endif
                     endif
                  else
                     styhi(jmax+1) = stylo(jmax+1)
#ifdef NOVERTICALINFLOW
!c     Hack for no vertical inflow velocity
                     if ((n+L-1).eq.ZVEL)then
                        styhi(jmax+1) = 0.0d0
                     endif
#endif
                  endif
               else if (bc(2,2).eq.REFLECT_EVEN) then
                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2).eq.REFLECT_ODD) then
                  stylo(jmax+1) = 0.0d0
                  styhi(jmax+1) = 0.0d0
               end if

               if ( velpred .eq. 1 ) then
                  do j = jmin, jmax+1
                     lty = stylo(j) .le. zero  .and.  styhi(j) .ge. 0.0d0
                     lty = lty .or. (abs(stylo(j)+styhi(j)) .lt. eps)
                     sty = merge(stylo(j),styhi(j),(stylo(j)+styhi(j)) .ge. 0.0d0)
                     ystate(i,j,k,L) = merge(0.0d0,sty,lty)
                  end do
               else
                  do j=jmin,jmax+1
                     ystate(i,j,k,L) = merge(stylo(j),styhi(j),vedge(i,j,k) .ge. 0.0d0)
                     ystate(i,j,k,L) = merge(half*(stylo(j)+styhi(j)),ystate(i,j,k,L), &
                         abs(vedge(i,j,k)).lt.eps)
                  end do
               end if
            end do
         end do
      end if
!c
!c     compute the zedge states
!c
      if ((velpred.ne.1) .or. ((n+L-1).eq.ZVEL)) then
         do j = jmin,jmax
            do i = imin,imax
               
               do k = kmin-1,kmax+1
                  if (uad(i,j,k)*uad(i+1,j,k).lt.0.d0) then
                      ubar = 0.5d0*(uad(i,j,k)+uad(i+1,j,k))
                      if (ubar.lt.0.d0) then
                          inc = 1
                      else
                          inc = 0
                      endif
                      tr1 = ubar*(s(i+inc,j,k,L)-s(i+inc-1,j,k,L))*ihx
                  else
                      tr1 = half* &
                      (uad(i+1,j,k) + uad(i,j,k)) * &
                      (xedge(i+1,j,k) - xedge(i,j,k)) *ihx
                  endif

                  if (vad(i,j,k)*vad(i,j+1,k).lt.0.d0) then
                      vbar = 0.5d0*(vad(i,j,k)+vad(i,j+1,k))
                      if (vbar.lt.0.d0) then
                          inc = 1
                      else
                          inc = 0
                      endif
                      tr2 = vbar*(s(i,j+inc,k,L)-s(i,j+inc-1,k,L))*ihy
                  else
                      tr2 = half* &
                      (vad(i,j+1,k) + vad(i,j,k)) * &
                      (yedge(i,j+1,k) - yedge(i,j,k)) *ihy
                  endif

                  if (ppm_type .gt. 0) then
                     stzlo(k+1)= Ipz(i,j,k) &
                         - dth*tr1 - dth*tr2 &
                         + dth*tf(i,j,k,L)
                     stzhi(k)  = Imz(i,j,k) &
                         - dth*tr1 - dth*tr2 &
                         + dth*tf(i,j,k,L)
                  else
                     stzlo(k+1)= s(i,j,k,L) + (half-dthz*w(i,j,k))*sz(i,j,k) &
                         - dth*tr1 - dth*tr2 &
                         + dth*tf(i,j,k,L) 
                     stzhi(k)  = s(i,j,k,L) - (half+dthz*w(i,j,k))*sz(i,j,k) &
                         - dth*tr1 - dth*tr2 &
                         + dth*tf(i,j,k,L)
                  end if
               end do

               if (bc(3,1).eq.EXT_DIR .and. velpred.eq.1) then
                  stzlo(kmin) = s(i,j,kmin-1,L)
                  stzhi(kmin) = s(i,j,kmin-1,L)
               else if (bc(3,1).eq.EXT_DIR .and. wad(i,j,kmin).ge.0.0d0) then
                  stzlo(kmin) = s(i,j,kmin-1,L)
                  stzhi(kmin) = s(i,j,kmin-1,L)
               else if (bc(3,1).eq.EXT_DIR .and. wad(i,j,kmin).lt.0.0d0) then
                  stzlo(kmin) = stzhi(kmin)
               else if (bc(3,1).eq.FOEXTRAP.or.bc(3,1).eq.HOEXTRAP) then
                  if ((n+L-1).eq.ZVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWZINFLOW
!c     prevent backflow
                        stzhi(kmin) = MIN(stzhi(kmin),0.0d0)
#endif
                        stzlo(kmin) = stzhi(kmin)
                     else
                        if (wad(i,j,kmin).ge.0.0d0) then
#ifndef ALLOWZINFLOW
!c     prevent backflow
                           stzhi(kmin) = MIN(stzhi(kmin),0.0d0)
#endif
                           stzlo(kmin) = stzhi(kmin)
                        endif
                     endif
                  else
                     stzlo(kmin) = stzhi(kmin)
                  endif
               else if (bc(3,1).eq.REFLECT_EVEN) then
                  stzlo(kmin) = stzhi(kmin)
               else if (bc(3,1).eq.REFLECT_ODD) then
                  stzlo(kmin) = 0.0d0
                  stzhi(kmin) = 0.0d0
               end if
               if (bc(3,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stzlo(kmax+1) = s(i,j,kmax+1,L)
                  stzhi(kmax+1) = s(i,j,kmax+1,L)
               else if (bc(3,2).eq.EXT_DIR .and. wad(i,j,kmax+1).le.0.0d0) then
                  stzlo(kmax+1) = s(i,j,kmax+1,L)
                  stzhi(kmax+1) = s(i,j,kmax+1,L)
               else if (bc(3,2).eq.EXT_DIR .and. wad(i,j,kmax+1).gt.0.0d0) then
                  stzhi(kmax+1) = stzlo(kmax+1)
               else if (bc(3,2).eq.FOEXTRAP.or.bc(3,2).eq.HOEXTRAP) then
                  if ((n+L-1).eq.ZVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWZINFLOW
!c     prevent backflow
                        stzlo(kmax+1) = MAX(stzlo(kmax+1),0.0d0)
#endif
                        stzhi(kmax+1) = stzlo(kmax+1)
                     else
                        if (wad(i,j,kmax+1).le.0.0d0) then
#ifndef ALLOWZINFLOW
!c     prevent backflow
                           stzlo(kmax+1) = MAX(stzlo(kmax+1),0.0d0)
#endif
                           stzhi(kmax+1) = stzlo(kmax+1)
                        endif
                     endif
                  else
                     stzhi(kmax+1) = stzlo(kmax+1)
                  endif
               else if (bc(3,2).eq.REFLECT_EVEN) then
                  stzhi(kmax+1) = stzlo(kmax+1)
               else if (bc(3,2).eq.REFLECT_ODD) then
                  stzlo(kmax+1) = 0.0d0
                  stzhi(kmax+1) = 0.0d0
               end if

               if ( velpred .eq. 1 ) then
                  do k = kmin,kmax+1
                     ltz = stzlo(k) .le. zero  .and.  stzhi(k) .ge. 0.0d0
                     ltz = ltz .or. (abs(stzlo(k)+stzhi(k)) .lt. eps)
                     stz = merge(stzlo(k),stzhi(k),(stzlo(k)+stzhi(k)) .ge. 0.0d0)
                     zstate(i,j,k,L) = merge(zero,stz,ltz)
                  end do
               else
                  do k = kmin,kmax+1
                     zstate(i,j,k,L) = merge(stzlo(k),stzhi(k),wedge(i,j,k) .ge. 0.0d0)
                     zstate(i,j,k,L) = merge(half*(stzlo(k)+stzhi(k)),zstate(i,j,k,L), &
                         abs(wedge(i,j,k)).lt.eps)
                  end do
               end if
            end do
         end do
      end if

      end if
      
      end do

      end subroutine estate_premac

      subroutine estate_fpu(lo,hi,&
         s,s_lo,s_hi,&
         tf,tf_lo,tf_hi,&
         divu,divu_lo,divu_hi,&
         xlo,xlo_lo,xlo_hi,&
         xhi,xhi_lo,xhi_hi,&
         sx,sx_lo,sx_hi,&
         xedge,xedge_lo,xedge_hi,&
         uedge,uedge_lo,uedge_hi,&
         xstate,xstate_lo,xstate_hi,&
         Imx,Imx_lo,Imx_hi,&
         Ipx,Ipx_lo,Ipx_hi,&
         sedgex,sedgex_lo,sedgex_hi,&
         ylo,ylo_lo,ylo_hi,&
         yhi,yhi_lo,yhi_hi,&
         sy,sy_lo,sy_hi,&
         yedge,yedge_lo,yedge_hi,&
         vedge,vedge_lo,vedge_hi,&
         ystate,ystate_lo,ystate_hi,&
         Imy,Imy_lo,Imy_hi,&
         Ipy,Ipy_lo,Ipy_hi,&
         sedgey,sedgey_lo,sedgey_hi,&
         zlo,zlo_lo,zlo_hi,&
         zhi,zhi_lo,zhi_hi,&
         sz,sz_lo,sz_hi,&
         zedge,zedge_lo,zedge_hi,&
         wedge,wedge_lo,wedge_hi,&
         zstate,zstate_lo,zstate_hi,&
         Imz,Imz_lo,Imz_hi,&
         Ipz,Ipz_lo,Ipz_hi,&
         sedgez,sedgez_lo,sedgez_hi,&
         dsvl,dsvl_lo,dsvl_hi,&
         sm,sm_lo,sm_hi,&
         sp,sp_lo,sp_hi,&
         xylo,xylo_lo,xylo_hi,&
         xzlo,xzlo_lo,xzlo_hi,&
         yxlo,yxlo_lo,yxlo_hi,&
         yzlo,yzlo_lo,yzlo_hi,&
         zxlo,zxlo_lo,zxlo_hi,&
         zylo,zylo_lo,zylo_hi,&
         xyhi,xyhi_lo,xyhi_hi,&
         xzhi,xzhi_lo,xzhi_hi,&
         yxhi,yxhi_lo,yxhi_hi,&
         yzhi,yzhi_lo,yzhi_hi,&
         zxhi,zxhi_lo,zxhi_hi,&
         zyhi,zyhi_lo,zyhi_hi,&
         corner_couple, &
         bc, dt, dx, n, nc, use_minion, iconserv, ppm_type)
!c
!c     This subroutine computes edges states, right now it uses
!c     a lot of memory, but there becomes a trade off between
!c     simplicity-efficiency in the new way of computing states
!c     and complexity in the old way.  By eliminating loops over
!c     state components though, the new way uses much less memory.
!c
!c     This routine differs from the default ESTATE function above in that
!c     it assumes that the edge velocities are valid in a grow cell outside
!c     the box, and no *ad (unprojected) velocities are used.  This routine
!c     will fail if the UMAC coming in hasn't been "fillpatched"
!c

      implicit none

      integer, intent(in) :: nc, use_minion, iconserv(nc), ppm_type, bc(SDIM,2,nc), n
      real(rt), intent(in) :: dt, dx(SDIM)

      integer, dimension(3), intent(in) :: s_lo,s_hi,tf_lo,tf_hi,divu_lo,divu_hi,&
           xlo_lo,xlo_hi,xhi_lo,xhi_hi,sx_lo,sx_hi,&
           uedge_lo,uedge_hi,xedge_lo,xedge_hi,xstate_lo,xstate_hi,Imx_lo,Imx_hi,Ipx_lo,Ipx_hi,sedgex_lo,sedgex_hi,&
           ylo_lo,ylo_hi,yhi_lo,yhi_hi,sy_lo,sy_hi,&
           vedge_lo,vedge_hi,yedge_lo,yedge_hi,ystate_lo,ystate_hi,Imy_lo,Imy_hi,Ipy_lo,Ipy_hi,sedgey_lo,sedgey_hi,&
           zlo_lo,zlo_hi,zhi_lo,zhi_hi,sz_lo,sz_hi,&
           wedge_lo,wedge_hi,zedge_lo,zedge_hi,zstate_lo,zstate_hi,Imz_lo,Imz_hi,Ipz_lo,Ipz_hi,sedgez_lo,sedgez_hi,&
           dsvl_lo,dsvl_hi,sm_lo,sm_hi,sp_lo,sp_hi,lo,hi

      integer, dimension(3), intent(in) :: xylo_lo,xylo_hi,xzlo_lo,xzlo_hi,yxlo_lo,yxlo_hi, &
                                           yzlo_lo,yzlo_hi,zxlo_lo,zxlo_hi,zylo_lo,zylo_hi, &
                                           xyhi_lo,xyhi_hi,xzhi_lo,xzhi_hi,yxhi_lo,yxhi_hi, &
                                           yzhi_lo,yzhi_hi,zxhi_lo,zxhi_hi,zyhi_lo,zyhi_hi

      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc)
      real(rt), intent(in) :: tf(tf_lo(1):tf_hi(1),tf_lo(2):tf_hi(2),tf_lo(3):tf_hi(3),nc)
      real(rt), intent(in) :: divu(divu_lo(1):divu_hi(1),divu_lo(2):divu_hi(2),divu_lo(3):divu_hi(3))
      real(rt), intent(inout) :: xlo(xlo_lo(1):xlo_hi(1),xlo_lo(2):xlo_hi(2),xlo_lo(3):xlo_hi(3))
      real(rt), intent(inout) :: xhi(xhi_lo(1):xhi_hi(1),xhi_lo(2):xhi_hi(2),xhi_lo(3):xhi_hi(3))
      real(rt), intent(inout) :: sx(sx_lo(1):sx_hi(1),sx_lo(2):sx_hi(2),sx_lo(3):sx_hi(3))
      real(rt), intent(in) :: uedge(uedge_lo(1):uedge_hi(1),uedge_lo(2):uedge_hi(2),uedge_lo(3):uedge_hi(3))
      real(rt), intent(inout) :: xedge(xedge_lo(1):xedge_hi(1),xedge_lo(2):xedge_hi(2),xedge_lo(3):xedge_hi(3))
      real(rt), intent(inout) :: xstate(xstate_lo(1):xstate_hi(1),xstate_lo(2):xstate_hi(2),xstate_lo(3):xstate_hi(3),nc) ! result
      real(rt), intent(inout) :: Imx(Imx_lo(1):Imx_hi(1),Imx_lo(2):Imx_hi(2),Imx_lo(3):Imx_hi(3))
      real(rt), intent(inout) :: Ipx(Ipx_lo(1):Ipx_hi(1),Ipx_lo(2):Ipx_hi(2),Ipx_lo(3):Ipx_hi(3))
      real(rt), intent(inout) :: sedgex(sedgex_lo(1):sedgex_hi(1),sedgex_lo(2):sedgex_hi(2),sedgex_lo(3):sedgex_hi(3))
      real(rt), intent(inout) :: ylo(ylo_lo(1):ylo_hi(1),ylo_lo(2):ylo_hi(2),ylo_lo(3):ylo_hi(3))
      real(rt), intent(inout) :: yhi(yhi_lo(1):yhi_hi(1),yhi_lo(2):yhi_hi(2),yhi_lo(3):yhi_hi(3))
      real(rt), intent(inout) :: sy(sy_lo(1):sy_hi(1),sy_lo(2):sy_hi(2),sy_lo(3):sy_hi(3))
      real(rt), intent(in) :: vedge(vedge_lo(1):vedge_hi(1),vedge_lo(2):vedge_hi(2),vedge_lo(3):vedge_hi(3))
      real(rt), intent(inout) :: yedge(yedge_lo(1):yedge_hi(1),yedge_lo(2):yedge_hi(2),yedge_lo(3):yedge_hi(3))
      real(rt), intent(inout) :: ystate(ystate_lo(1):ystate_hi(1),ystate_lo(2):ystate_hi(2),ystate_lo(3):ystate_hi(3),nc) ! result
      real(rt), intent(inout) :: Imy(Imy_lo(1):Imy_hi(1),Imy_lo(2):Imy_hi(2),Imy_lo(3):Imy_hi(3))
      real(rt), intent(inout) :: Ipy(Ipy_lo(1):Ipy_hi(1),Ipy_lo(2):Ipy_hi(2),Ipy_lo(3):Ipy_hi(3))
      real(rt), intent(inout) :: sedgey(sedgey_lo(1):sedgey_hi(1),sedgey_lo(2):sedgey_hi(2),sedgey_lo(3):sedgey_hi(3))
      
      real(rt), intent(inout) :: zlo(zlo_lo(1):zlo_hi(1),zlo_lo(2):zlo_hi(2),zlo_lo(3):zlo_hi(3))
      real(rt), intent(inout) :: zhi(zhi_lo(1):zhi_hi(1),zhi_lo(2):zhi_hi(2),zhi_lo(3):zhi_hi(3))
      real(rt), intent(inout) :: sz(sz_lo(1):sz_hi(1),sz_lo(2):sz_hi(2),sz_lo(3):sz_hi(3))
      real(rt), intent(in) :: wedge(wedge_lo(1):wedge_hi(1),wedge_lo(2):wedge_hi(2),wedge_lo(3):wedge_hi(3))
      real(rt), intent(inout) :: zedge(zedge_lo(1):zedge_hi(1),zedge_lo(2):zedge_hi(2),zedge_lo(3):zedge_hi(3))
      real(rt), intent(inout) :: zstate(zstate_lo(1):zstate_hi(1),zstate_lo(2):zstate_hi(2),zstate_lo(3):zstate_hi(3),nc) ! result
      real(rt), intent(inout) :: Imz(Imz_lo(1):Imz_hi(1),Imz_lo(2):Imz_hi(2),Imz_lo(3):Imz_hi(3))
      real(rt), intent(inout) :: Ipz(Ipz_lo(1):Ipz_hi(1),Ipz_lo(2):Ipz_hi(2),Ipz_lo(3):Ipz_hi(3))
      real(rt), intent(inout) :: sedgez(sedgez_lo(1):sedgez_hi(1),sedgez_lo(2):sedgez_hi(2),sedgez_lo(3):sedgez_hi(3))
      
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2),sp_lo(3):sp_hi(3))
      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2),dsvl_lo(3):dsvl_hi(3))

      real(rt), intent(inout) :: xylo(xylo_lo(1):xylo_hi(1),xylo_lo(2):xylo_hi(2),xylo_lo(3):xylo_hi(3))
      real(rt), intent(inout) :: xzlo(xzlo_lo(1):xzlo_hi(1),xzlo_lo(2):xzlo_hi(2),xzlo_lo(3):xzlo_hi(3))
      real(rt), intent(inout) :: yxlo(yxlo_lo(1):yxlo_hi(1),yxlo_lo(2):yxlo_hi(2),yxlo_lo(3):yxlo_hi(3))
      real(rt), intent(inout) :: yzlo(yzlo_lo(1):yzlo_hi(1),yzlo_lo(2):yzlo_hi(2),yzlo_lo(3):yzlo_hi(3))
      real(rt), intent(inout) :: zxlo(zxlo_lo(1):zxlo_hi(1),zxlo_lo(2):zxlo_hi(2),zxlo_lo(3):zxlo_hi(3))
      real(rt), intent(inout) :: zylo(zylo_lo(1):zylo_hi(1),zylo_lo(2):zylo_hi(2),zylo_lo(3):zylo_hi(3))

      real(rt), intent(inout) :: xyhi(xyhi_lo(1):xyhi_hi(1),xyhi_lo(2):xyhi_hi(2),xyhi_lo(3):xyhi_hi(3))
      real(rt), intent(inout) :: xzhi(xzhi_lo(1):xzhi_hi(1),xzhi_lo(2):xzhi_hi(2),xzhi_lo(3):xzhi_hi(3))
      real(rt), intent(inout) :: yxhi(yxhi_lo(1):yxhi_hi(1),yxhi_lo(2):yxhi_hi(2),yxhi_lo(3):yxhi_hi(3))
      real(rt), intent(inout) :: yzhi(yzhi_lo(1):yzhi_hi(1),yzhi_lo(2):yzhi_hi(2),yzhi_lo(3):yzhi_hi(3))
      real(rt), intent(inout) :: zxhi(zxhi_lo(1):zxhi_hi(1),zxhi_lo(2):zxhi_hi(2),zxhi_lo(3):zxhi_hi(3))
      real(rt), intent(inout) :: zyhi(zyhi_lo(1):zyhi_hi(1),zyhi_lo(2):zyhi_hi(2),zyhi_lo(3):zyhi_hi(3))
      
      real(rt) :: stxlo(lo(1)-2:hi(1)+2)
      real(rt) :: stxhi(lo(1)-2:hi(1)+2)
      real(rt) :: stylo(lo(2)-2:hi(2)+2)
      real(rt) :: styhi(lo(2)-2:hi(2)+2)
      real(rt) :: stzlo(lo(3)-2:hi(3)+2)
      real(rt) :: stzhi(lo(3)-2:hi(3)+2)
      real(rt) :: hx, hy, hz, dth, dthx, dthy, dthz, ihx, ihy, ihz
      real(rt) :: tr,tr1,tr2,ubar,vbar,wbar,stx,sty,stz,fu,fv,fw,eps,eps_for_bc,st
      real(rt) ::  dt3, dt3x, dt3y, dt3z, dt4, dt4x, dt4y, dt4z
      real(rt) ::  dt6, dt6x, dt6y, dt6z
      integer  :: i,j,k,L,imin,jmin,kmin,imax,jmax,kmax,inc,corner_couple
      parameter( eps        = 1.0D-6 )
      parameter( eps_for_bc = 1.0D-10 )

      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      dthz = half*dt / dx(3)
      dt3  = dt / 3.0d0
      dt3x = dt3 / dx(1)
      dt3y = dt3 / dx(2)
      dt3z = dt3 / dx(3)
      dt4  = dt / 4.0d0
      dt4x = dt4 / dx(1)
      dt4y = dt4 / dx(2)
      dt4z = dt4 / dx(3)
      dt6  = dt / 6.0d0
      dt6x = dt6 / dx(1)
      dt6y = dt6 / dx(2)
      dt6z = dt6 / dx(3)
      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      ihx  = 1.0d0/dx(1)
      ihy  = 1.0d0/dx(2)
      ihz  = 1.0d0/dx(3)

      do L=1,nc
      
!c
!c     compute the slopes
!c
      if (ppm_type .gt. 0) then
         call ppm_fpu(lo, hi,&
              s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
              uedge,uedge_lo,uedge_hi,&
              vedge,vedge_lo,vedge_hi,&
              wedge,wedge_lo,wedge_hi,&
              Ipx,Ipx_lo,Ipx_hi,&
              Imx,Imx_lo,Imx_hi,&
              Ipy,Ipy_lo,Ipy_hi,&
              Imy,Imy_lo,Imy_hi,&
              Ipz,Ipz_lo,Ipz_hi,&
              Imz,Imz_lo,Imz_hi,&
              sm,sm_lo,sm_hi,&
              sp,sp_lo,sp_hi,&
              dsvl,dsvl_lo,dsvl_hi,&
              sedgex,sedgex_lo,sedgex_hi,&
              sedgey,sedgey_lo,sedgey_hi,&
              sedgez,sedgez_lo,sedgez_hi,&
              dx, dt, bc(1,1,L), eps_for_bc, ppm_type) 
      else
         call slopes(ALL, &
                          lo,hi,&
                          s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
                          sx,sx_lo,sx_hi,&
                          sy,sy_lo,sy_hi,&
                          sz,sz_lo,sz_hi,&
                          bc(1,1,L))
      end if
!c
!c     trace the state to the cell edges
!c
      if (ppm_type .gt. 0) then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  xlo(i,j,k) = Ipx(i-1,j,k)
                  xhi(i,j,k) = Imx(i  ,j,k)
               end do
            end do
         end do
      else
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  xlo(i,j,k) = s(i-1,j,k,L) + (half  - dthx*uedge(i,j,k))*sx(i-1,j,k)
                  xhi(i,j,k) = s(i,  j,k,L) + (-half - dthx*uedge(i,j,k))*sx(i,  j,k)
               end do
            end do
         end do
      end if

      if (use_minion.eq.1)then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  xlo(i,j,k) = xlo(i,j,k) + dth*tf(i-1,j,k,L)
                  xhi(i,j,k) = xhi(i,j,k) + dth*tf(i,  j,k,L)
               end do
            end do
         end do
         if (iconserv(L) .eq. 1) then
           do k = kmin-1,kmax+1
             do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  xlo(i,j,k) = xlo(i,j,k) - dth*s(i-1,j,k,L)*divu(i-1,j,k)
                  xhi(i,j,k) = xhi(i,j,k) - dth*s(i  ,j,k,L)*divu(i,  j,k)
               end do
             end do
           end do
         end if
      end if

      call trans_xbc(lo,hi,&
           s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
           xlo,xlo_lo,xlo_hi,&
           xhi,xhi_lo,xhi_hi,&
           uedge,uedge_lo,uedge_hi,&
           n+L-1, bc(1,1,L), eps_for_bc,.false.,.false.)

      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               fu  = merge(zero,one,abs(uedge(i,j,k)).lt.eps)
               stx = merge(xlo(i,j,k),xhi(i,j,k),uedge(i,j,k) .ge. 0.0d0)
               xedge(i,j,k) = fu*stx + (one - fu)*half*(xhi(i,j,k)+xlo(i,j,k))
            end do
         end do
      end do

      if (ppm_type .gt. 0) then
         do k = kmin-1,kmax+1
            do j = jmin,  jmax+1
               do i = imin-1,imax+1
                  ylo(i,j,k) = Ipy(i,j-1,k)
                  yhi(i,j,k) = Imy(i,j  ,k)
               end do
            end do
         end do
      else
         do k = kmin-1,kmax+1
            do j = jmin,  jmax+1
               do i = imin-1,imax+1
                  ylo(i,j,k) = s(i,j-1,k,L) + (half  - dthy*vedge(i,j,k))*sy(i,j-1,k)
                  yhi(i,j,k) = s(i,j, k,L)  + (-half - dthy*vedge(i,j,k))*sy(i,j, k)
               end do
            end do
         end do
      end if

      if (use_minion.eq.1)then
         do k = kmin-1,kmax+1
            do j = jmin, jmax+1
               do i = imin-1,  imax+1
                  ylo(i,j,k) = ylo(i,j,k) + dth*tf(i,j-1,k,L)
                  yhi(i,j,k) = yhi(i,j,k) + dth*tf(i,j,  k,L)
               end do
            end do
         end do
         if (iconserv(L) .eq. 1) then
           do k = kmin-1,kmax+1
             do j = jmin, jmax+1
               do i = imin-1,  imax+1
                  ylo(i,j,k) = ylo(i,j,k) - dth*s(i,j-1,k,L)*divu(i,j-1,k)
                  yhi(i,j,k) = yhi(i,j,k) - dth*s(i,j  ,k,L)*divu(i,j,  k)
               end do
             end do
           end do
         end if
      end if

      call trans_ybc(lo,hi,&
           s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
           ylo,ylo_lo,ylo_hi,&
           yhi,yhi_lo,yhi_hi,&
           vedge,vedge_lo,vedge_hi,&
           n+L-1, bc(1,1,L), eps_for_bc,.false.,.false.)

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               fv  = merge(zero,one,abs(vedge(i,j,k)).lt.eps)
               sty = merge(ylo(i,j,k),yhi(i,j,k),vedge(i,j,k) .ge. 0.0d0)
               yedge(i,j,k) = fv*sty + (one - fv)*half*(yhi(i,j,k)+ylo(i,j,k))
            end do
         end do
      end do

      if (ppm_type .gt. 0) then
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  zlo(i,j,k) = Ipz(i,j,k-1)
                  zhi(i,j,k) = Imz(i,j,k  )
               end do
            end do
         end do
      else
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  zlo(i,j,k) = s(i,j,k-1,L) + (half  - dthz*wedge(i,j,k))*sz(i,j,k-1)
                  zhi(i,j,k) = s(i,j,k  ,L) + (-half - dthz*wedge(i,j,k))*sz(i,j,k  )
               end do
            end do
         end do
      end if

      if (use_minion.eq.1)then
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,  imax+1
                  zlo(i,j,k) = zlo(i,j,k) + dth*tf(i,j,k-1,L)
                  zhi(i,j,k) = zhi(i,j,k) + dth*tf(i,j,k,L)
               end do
            end do
         end do
         if (iconserv(L) .eq. 1) then
           do k = kmin,kmax+1
             do j = jmin-1,jmax+1
               do i = imin-1,  imax+1
                  zlo(i,j,k) = zlo(i,j,k) - dth*s(i,j,k-1,L)*divu(i,j,k-1)
                  zhi(i,j,k) = zhi(i,j,k) - dth*s(i,j,k  ,L)*divu(i,j,k  )
               end do
             end do
           end do
         end if
      end if

      call trans_zbc(lo,hi,&
           s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
           zlo,zlo_lo,zlo_hi,&
           zhi,zhi_lo,zhi_hi,&
           wedge,wedge_lo,wedge_hi,&
           n+L-1, bc(1,1,L), eps_for_bc,.false.,.false.)
      
      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               fw  = merge(zero,one,abs(wedge(i,j,k)).lt.eps)
               stz = merge(zlo(i,j,k),zhi(i,j,k),wedge(i,j,k) .ge. 0.0d0)
               zedge(i,j,k) = fw*stz + (one-fw)*half*(zhi(i,j,k)+zlo(i,j,k))
            end do
         end do
      end do

      if (corner_couple .ne. 0) then

!c
!c     NEW CORNER COUPLING CODE
!c

!c
!c     compute the corner-coupled terms:
!c     xylo/hi, xzlo/hi, yxlo/hi, yzlo/hi, zxlo/hi, zylo/hi
!c

!c     loop over appropriate xy faces
      if (iconserv(L).eq.1) then
         do k=kmin-1,kmax+1
            do j=jmin,jmax
               do i=imin,imax+1
                  xylo(i,j,k) = xlo(i,j,k) &
                      - dt3y*(yedge(i-1,j+1,k)*vedge(i-1,j+1,k) &
                      - yedge(i-1,j,k)*vedge(i-1,j,k)) &
                      - dt3*s(i-1,j,k,L)*divu(i-1,j,k) &
                      + dt3y*s(i-1,j,k,L)*(vedge(i-1,j+1,k)-vedge(i-1,j,k))
                 xyhi(i,j,k) = xhi(i,j,k) &
                      - dt3y*(yedge(i  ,j+1,k)*vedge(i  ,j+1,k) &
                      - yedge(i  ,j,k)*vedge(i  ,j,k)) &
                      - dt3*s(i,j,k,L)*divu(i,j,k) &
                      + dt3y*s(i,j,k,L)*(vedge(i,j+1,k)-vedge(i,j,k))
               end do
            end do
         end do
      else
         do k=kmin-1,kmax+1
            do j=jmin,jmax
               do i=imin,imax+1
                  xylo(i,j,k) = xlo(i,j,k) &
                      - dt6y*(vedge(i-1,j+1,k)+vedge(i-1,j,k)) &
                      *(yedge(i-1,j+1,k)-yedge(i-1,j,k))
                  xyhi(i,j,k) = xhi(i,j,k) &
                      - dt6y*(vedge(i  ,j+1,k)+vedge(i  ,j,k)) &
                      *(yedge(i  ,j+1,k)-yedge(i  ,j,k))
               end do
            end do
         end do
      end if

!c     boundary conditions

     call trans_xbc(lo,hi,&
           s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
           xylo,xylo_lo,xylo_hi,&
           xyhi,xyhi_lo,xyhi_hi,&
           uedge,uedge_lo,uedge_hi,&
           n+L-1, bc(1,1,L), eps_for_bc,.true.,.false.)

!c     upwind
      do k=kmin-1,kmax+1
         do j=jmin,jmax
            do i=imin,imax+1
               fu  = merge(zero,one,abs(uedge(i,j,k)).lt.eps)
               stx = merge(xylo(i,j,k),xyhi(i,j,k),uedge(i,j,k) .ge. 0.0d0)
               xylo(i,j,k) = fu*stx + (one - fu)*half*(xyhi(i,j,k)+xylo(i,j,k))
            end do
         end do
      end do

!c     loop over appropriate xz faces
      if (iconserv(L).eq.1) then
         do k=kmin,kmax
            do j=jmin-1,jmax+1
               do i=imin,imax+1
                  xzlo(i,j,k) = xlo(i,j,k) &
                      - dt3z*(zedge(i-1,j,k+1)*wedge(i-1,j,k+1) &
                      - zedge(i-1,j,k)*wedge(i-1,j,k)) &
                      - dt3*s(i-1,j,k,L)*divu(i-1,j,k) &
                      + dt3z*s(i-1,j,k,L)*(wedge(i-1,j,k+1)-wedge(i-1,j,k))
                 xzhi(i,j,k) = xhi(i,j,k) &
                      - dt3z*(zedge(i  ,j,k+1)*wedge(i  ,j,k+1) &
                      - zedge(i  ,j,k)*wedge(i  ,j,k)) &
                      - dt3*s(i,j,k,L)*divu(i,j,k) &
                      + dt3z*s(i,j,k,L)*(wedge(i,j,k+1)-wedge(i,j,k))
               end do
            end do
         end do
      else
         do k=kmin,kmax
            do j=jmin-1,jmax+1
               do i=imin,imax+1
                  xzlo(i,j,k) = xlo(i,j,k) &
                      - dt6z*(wedge(i-1,j,k+1)+wedge(i-1,j,k)) &
                      *(zedge(i-1,j,k+1)-zedge(i-1,j,k))
                  xzhi(i,j,k) = xhi(i,j,k) &
                      - dt6z*(wedge(i  ,j,k+1)+wedge(i  ,j,k)) &
                      *(zedge(i  ,j,k+1)-zedge(i  ,j,k))
               end do
            end do
         end do
      end if

!c     boundary conditions
     call trans_xbc(lo,hi,&
           s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
           xzlo,xzlo_lo,xzlo_hi,&
           xzhi,xzhi_lo,xzhi_hi,&
           uedge,uedge_lo,uedge_hi,&
           n+L-1, bc(1,1,L), eps_for_bc,.false.,.true.)
           
!c     upwind
      do k=kmin,kmax
         do j=jmin-1,jmax+1
            do i=imin,imax+1
               fu  = merge(zero,one,abs(uedge(i,j,k)).lt.eps)
               stx = merge(xzlo(i,j,k),xzhi(i,j,k),uedge(i,j,k) .ge. 0.0d0)
               xzlo(i,j,k) = fu*stx + (one - fu)*half*(xzhi(i,j,k)+xzlo(i,j,k))
            end do
         end do
      end do

!c     loop over appropriate yx faces
      if (iconserv(L).eq.1) then
         do k=kmin-1,kmax+1
            do j=jmin,jmax+1
               do i=imin,imax
                  yxlo(i,j,k) = ylo(i,j,k) &
                      - dt3x*(xedge(i+1,j-1,k)*uedge(i+1,j-1,k) &
                      - xedge(i,j-1,k)*uedge(i,j-1,k)) &
                      - dt3*s(i,j-1,k,L)*divu(i,j-1,k) &
                      + dt3x*s(i,j-1,k,L)*(uedge(i+1,j-1,k)-uedge(i,j-1,k))
                 yxhi(i,j,k) = yhi(i,j,k) &
                      - dt3x*(xedge(i+1,j  ,k)*uedge(i+1,j  ,k) &
                      - xedge(i,j  ,k)*uedge(i,j  ,k)) &
                      - dt3*s(i,j,k,L)*divu(i,j,k) &
                      + dt3x*s(i,j,k,L)*(uedge(i+1,j,k)-uedge(i,j,k))
               end do
            end do
         end do
      else
         do k=kmin-1,kmax+1
            do j=jmin,jmax+1
               do i=imin,imax
                  yxlo(i,j,k) = ylo(i,j,k) &
                      - dt6x*(uedge(i+1,j-1,k)+uedge(i,j-1,k)) &
                      *(xedge(i+1,j-1,k)-xedge(i,j-1,k))
                  yxhi(i,j,k) = yhi(i,j,k) &
                      - dt6x*(uedge(i+1,j  ,k)+uedge(i,j  ,k)) &
                      *(xedge(i+1,j  ,k)-xedge(i,j  ,k))
                  
               end do
            end do
         end do
      end if

!c     boundary conditions

      call trans_ybc(lo,hi,&
           s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
           yxlo,yxlo_lo,yxlo_hi,&
           yxhi,yxhi_lo,yxhi_hi,&
           vedge,vedge_lo,vedge_hi,&
           n+L-1, bc(1,1,L), eps_for_bc,.true.,.false.)

!c     upwind
      do k=kmin-1,kmax+1
         do j=jmin,jmax+1
            do i=imin,imax
               fv  = merge(zero,one,abs(vedge(i,j,k)).lt.eps)
               sty = merge(yxlo(i,j,k),yxhi(i,j,k),vedge(i,j,k) .ge. 0.0d0)
               yxlo(i,j,k) = fv*sty + (one - fv)*half*(yxhi(i,j,k)+yxlo(i,j,k))
            end do
         end do
      end do

!c     loop over appropriate yz faces
      if (iconserv(L).eq.1) then
         do k=kmin,kmax
            do j=jmin,jmax+1
               do i=imin-1,imax+1
                  yzlo(i,j,k) = ylo(i,j,k) &
                      - dt3z*(zedge(i,j-1,k+1)*wedge(i,j-1,k+1) &
                      - zedge(i,j-1,k)*wedge(i,j-1,k)) &
                      - dt3*s(i,j-1,k,L)*divu(i,j-1,k) &
                      + dt3z*s(i,j-1,k,L)*(wedge(i,j-1,k+1)-wedge(i,j-1,k))
                 yzhi(i,j,k) = yhi(i,j,k) &
                      - dt3z*(zedge(i,j  ,k+1)*wedge(i,j  ,k+1) &
                      - zedge(i,j  ,k)*wedge(i,j  ,k)) &
                      - dt3*s(i,j,k,L)*divu(i,j,k) &
                      + dt3z*s(i,j,k,L)*(wedge(i,j,k+1)-wedge(i,j,k))
               end do
            end do
         end do
      else
         do k=kmin,kmax
            do j=jmin,jmax+1
               do i=imin-1,imax+1
                  yzlo(i,j,k) = ylo(i,j,k) &
                      - dt6z*(wedge(i,j-1,k+1)+wedge(i,j-1,k)) &
                      *(zedge(i,j-1,k+1)-zedge(i,j-1,k))
                  yzhi(i,j,k) = yhi(i,j,k) &
                      - dt6z*(wedge(i,j  ,k+1)+wedge(i,j  ,k)) &
                      *(zedge(i,j  ,k+1)-zedge(i,j  ,k))
               end do
            end do
         end do
      end if

!c     boundary conditions
      call trans_ybc(lo,hi,&
           s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
           yzlo,yzlo_lo,yzlo_hi,&
           yzhi,yzhi_lo,yzhi_hi,&
           vedge,vedge_lo,vedge_hi,&
           n+L-1, bc(1,1,L), eps_for_bc,.false.,.true.)
           
!c     upwind
      do k=kmin,kmax
         do j=jmin,jmax+1
            do i=imin-1,imax+1
               fv  = merge(zero,one,abs(vedge(i,j,k)).lt.eps)
               sty = merge(yzlo(i,j,k),yzhi(i,j,k),vedge(i,j,k) .ge. 0.0d0)
               yzlo(i,j,k) = fv*sty + (one - fv)*half*(yzhi(i,j,k)+yzlo(i,j,k))
            end do
         end do
      end do

!c     loop over appropriate zx faces
      if (iconserv(L).eq.1) then
         do k=kmin,kmax+1
            do j=jmin-1,jmax+1
               do i=imin,imax
                  zxlo(i,j,k) = zlo(i,j,k) &
                      - dt3x*(xedge(i+1,j,k-1)*uedge(i+1,j,k-1) &
                      - xedge(i,j,k-1)*uedge(i,j,k-1)) &
                      - dt3*s(i,j,k-1,L)*divu(i,j,k-1) &
                      + dt3x*s(i,j,k-1,L)*(uedge(i+1,j,k-1)-uedge(i,j,k-1))
                 zxhi(i,j,k) = zhi(i,j,k) &
                      - dt3x*(xedge(i+1,j,k  )*uedge(i+1,j,k  ) &
                      - xedge(i,j,k  )*uedge(i,j,k  )) &
                      - dt3*s(i,j,k,L)*divu(i,j,k) &
                      + dt3x*s(i,j,k,L)*(uedge(i+1,j,k)-uedge(i,j,k))
               end do
            end do
         end do
      else
         do k=kmin,kmax+1
            do j=jmin-1,jmax+1
               do i=imin,imax
                  zxlo(i,j,k) = zlo(i,j,k) &
                      - dt6x*(uedge(i+1,j,k-1)+uedge(i,j,k-1)) &
                      *(xedge(i+1,j,k-1)-xedge(i,j,k-1))
                  zxhi(i,j,k) = zhi(i,j,k) &
                      - dt6x*(uedge(i+1,j,k  )+uedge(i,j,k  )) &
                      *(xedge(i+1,j,k  )-xedge(i,j,k  ))
               end do
            end do
         end do
      end if

!c     boundary conditions
     call trans_zbc(lo,hi,&
           s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
           zxlo,zxlo_lo,zxlo_hi,&
           zxhi,zxhi_lo,zxhi_hi,&
           wedge,wedge_lo,wedge_hi,&
           n+L-1, bc(1,1,L), eps_for_bc,.true.,.false.)

!c     upwind
      do k=kmin,kmax+1
         do j=jmin-1,jmax+1
            do i=imin,imax
               fw  = merge(zero,one,abs(wedge(i,j,k)).lt.eps)
               stz = merge(zxlo(i,j,k),zxhi(i,j,k),wedge(i,j,k) .ge. 0.0d0)
               zxlo(i,j,k) = fw*stz + (one-fw)*half*(zxhi(i,j,k)+zxlo(i,j,k))
            end do
         end do
      end do

!c     loop over appropriate zy faces
      if (iconserv(L).eq.1) then
         do k=kmin,kmax+1
            do j=jmin,jmax
               do i=imin-1,imax+1
                  zylo(i,j,k) = zlo(i,j,k) &
                      - dt3y*(yedge(i,j+1,k-1)*vedge(i,j+1,k-1) &
                      - yedge(i,j,k-1)*vedge(i,j,k-1)) &
                      - dt3*s(i,j,k-1,L)*divu(i,j,k-1) &
                      + dt3y*s(i,j,k-1,L)*(vedge(i,j+1,k-1)-vedge(i,j,k-1))
                 zyhi(i,j,k) = zhi(i,j,k) &
                      - dt3y*(yedge(i,j+1,k  )*vedge(i,j+1,k  ) &
                      - yedge(i,j,k  )*vedge(i,j,k  )) &
                      - dt3*s(i,j,k,L)*divu(i,j,k) &
                      + dt3y*s(i,j,k,L)*(vedge(i,j+1,k)-vedge(i,j,k))
               end do
            end do
         end do
      else
         do k=kmin,kmax+1
            do j=jmin,jmax
               do i=imin-1,imax+1
                  zylo(i,j,k) = zlo(i,j,k) &
                      - dt6y*(vedge(i,j+1,k-1)+vedge(i,j,k-1)) &
                      *(yedge(i,j+1,k-1)-yedge(i,j,k-1))
                  zyhi(i,j,k) = zhi(i,j,k) &
                      - dt6y*(vedge(i,j+1,k  )+vedge(i,j,k  )) &
                      *(yedge(i,j+1,k  )-yedge(i,j,k  ))
               end do
            end do
         end do
      end if
         
!c     boundary conditions
     call trans_zbc(lo,hi,&
           s(s_lo(1),s_lo(2),s_lo(3),L),s_lo,s_hi,&
           zylo,zylo_lo,zylo_hi,&
           zyhi,zyhi_lo,zyhi_hi,&
           wedge,wedge_lo,wedge_hi,&
           n+L-1, bc(1,1,L), eps_for_bc,.false.,.true.)

!c     upwind
      do k=kmin,kmax+1
         do j=jmin,jmax
            do i=imin-1,imax+1
               fw  = merge(zero,one,abs(wedge(i,j,k)).lt.eps)
               stz = merge(zylo(i,j,k),zyhi(i,j,k),wedge(i,j,k) .ge. 0.0d0)
               zylo(i,j,k) = fw*stz + (one-fw)*half*(zyhi(i,j,k)+zylo(i,j,k))
            end do
         end do
      end do
      
!c
!c     compute the xedge states
!c
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax+1
                  
               if (iconserv(L).eq.1) then

                  stxlo(i) = xlo(i,j,k) &
                      - dthy*(yzlo(i-1,j+1,k  )*vedge(i-1,j+1,k  ) &
                      - yzlo(i-1,j,k)*vedge(i-1,j,k)) &
                      - dthz*(zylo(i-1,j  ,k+1)*wedge(i-1,j  ,k+1) &
                      - zylo(i-1,j,k)*wedge(i-1,j,k)) &
                      + dthy*s(i-1,j,k,L)*(vedge(i-1,j+1,k)-vedge(i-1,j,k)) &
                      + dthz*s(i-1,j,k,L)*(wedge(i-1,j,k+1)-wedge(i-1,j,k)) 
                  stxhi(i) = xhi(i,j,k) &
                      - dthy*(yzlo(i  ,j+1,k  )*vedge(i  ,j+1,  k) &
                      - yzlo(i  ,j,k)*vedge(i  ,j,k)) &
                      - dthz*(zylo(i  ,j  ,k+1)*wedge(i  ,j  ,k+1) &
                      - zylo(i  ,j,k)*wedge(i  ,j,k)) &
                      + dthy*s(i  ,j,k,L)*(vedge(i,j+1,k)-vedge(i,j,k)) &
                      + dthz*s(i  ,j,k,L)*(wedge(i,j,k+1)-wedge(i,j,k))
                
                  if (use_minion.eq.0) then
                     stxlo(i) = stxlo(i) - dth*s(i-1,j,k,L)*divu(i-1,j,k)
                     stxhi(i) = stxhi(i) - dth*s(i  ,j,k,L)*divu(i,  j,k)
                  end if

               else

                  stxlo(i) = xlo(i,j,k) &
                      - dt4y*(vedge(i-1,j+1,k  )+vedge(i-1,j,k))* &
                      (yzlo(i-1,j+1,k  )-yzlo(i-1,j,k)) &
                      - dt4z*(wedge(i-1,j  ,k+1)+wedge(i-1,j,k))* &
                      (zylo(i-1,j  ,k+1)-zylo(i-1,j,k))
                  stxhi(i) = xhi(i,j,k) &
                      - (dt4*ihy)*(vedge(i  ,j+1,k  )+vedge(i  ,j,k))* &
                      (yzlo(i  ,j+1,k  )-yzlo(i  ,j,k)) &
                      - (dt4*ihz)*(wedge(i  ,j  ,k+1)+wedge(i  ,j,k))* &
                      (zylo(i  ,j  ,k+1)-zylo(i  ,j,k))

               endif

               if (use_minion.eq.0) then
                  stxlo(i) = stxlo(i) + dth*tf(i-1,j,k,L)
                  stxhi(i) = stxhi(i) + dth*tf(i,  j,k,L)
               end if
               
            end do
            
            if (bc(1,1,L).eq.EXT_DIR .and. uedge(imin,j,k).ge.zero) then
               stxhi(imin) = s(imin-1,j,k,L)
               stxlo(imin) = s(imin-1,j,k,L)
            else if (bc(1,1,L).eq.EXT_DIR .and. uedge(imin,j,k).lt.zero) then
               stxlo(imin) = stxhi(imin)
            else if (bc(1,1,L).eq.FOEXTRAP.or.bc(1,1,L).eq.HOEXTRAP) then
               if (n.eq.XVEL) then
                  if (uedge(imin,j,k).ge.zero) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                     stxhi(imin) = MIN(stxhi(imin),zero)
#endif
                     stxlo(imin) = stxhi(imin)
                  endif
               else
                  stxlo(imin) = stxhi(imin)
               endif
            else if (bc(1,1,L).eq.REFLECT_EVEN) then
               stxlo(imin) = stxhi(imin)
            else if (bc(1,1,L).eq.REFLECT_ODD) then
               stxhi(imin) = zero
               stxlo(imin) = zero
            end if
            if (bc(1,2,L).eq.EXT_DIR .and. uedge(imax+1,j,k).le.zero) then
               stxlo(imax+1) = s(imax+1,j,k,L)
               stxhi(imax+1) = s(imax+1,j,k,L)
            else if (bc(1,2,L).eq.EXT_DIR .and. uedge(imax+1,j,k).gt.zero) then
               stxhi(imax+1) = stxlo(imax+1)
            else if (bc(1,2,L).eq.FOEXTRAP.or.bc(1,2,L).eq.HOEXTRAP) then
               if (n.eq.XVEL) then
                  if (uedge(imax+1,j,k).le.zero) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                     stxlo(imax+1) = MAX(stxlo(imax+1),zero)
#endif
                     stxhi(imax+1) = stxlo(imax+1)
                  endif
               else
                  stxhi(imax+1) = stxlo(imax+1)
               endif
            else if (bc(1,2,L).eq.REFLECT_EVEN) then
               stxhi(imax+1) = stxlo(imax+1)
            else if (bc(1,2,L).eq.REFLECT_ODD) then
               stxlo(imax+1) = zero
               stxhi(imax+1) = zero
            end if
            
            do i = imin, imax+1
               xstate(i,j,k,L) = merge(stxlo(i),stxhi(i),uedge(i,j,k) .ge. 0.0d0) 
               xstate(i,j,k,L) = merge(half*(stxlo(i)+stxhi(i)),xstate(i,j,k,L) &
                   ,abs(uedge(i,j,k)).lt.eps)
            end do
         end do
      end do
!c
!c     compute the yedge states
!c
      do k = kmin,kmax
         do i = imin,imax
            do j = jmin,jmax+1

               if (iconserv(L).eq.1) then

                  stylo(j) = ylo(i,j,k) &
                      - dthx*(xzlo(i+1,j-1,k  )*uedge(i+1,j-1,k  ) &
                      - xzlo(i,j-1,k)*uedge(i,j-1,k)) &
                      - dthz*(zxlo(i  ,j-1,k+1)*wedge(i  ,j-1,k+1) &
                      - zxlo(i,j-1,k)*wedge(i,j-1,k)) &
                      + dthx*s(i,j-1,k,L)*(uedge(i+1,j-1,k)-uedge(i,j-1,k)) &
                      + dthz*s(i,j-1,k,L)*(wedge(i,j-1,k+1)-wedge(i,j-1,k))
                  styhi(j) = yhi(i,j,k) &
                      - dthx*(xzlo(i+1,j  ,k  )*uedge(i+1,j  ,k  ) &
                      - xzlo(i,j  ,k)*uedge(i,j  ,k)) &
                      - dthz*(zxlo(i  ,j  ,k+1)*wedge(i  ,j  ,k+1) &
                      - zxlo(i,j  ,k)*wedge(i,j  ,k)) &
                      + dthx*s(i,j  ,k,L)*(uedge(i+1,j,k)-uedge(i,j,k)) &
                      + dthz*s(i,j  ,k,L)*(wedge(i,j,k+1)-wedge(i,j,k))
                  
                  if (use_minion .eq. 0) then
                     stylo(j) = stylo(j) - dth*s(i,j-1,k,L)*divu(i,j-1,k)
                     styhi(j) = styhi(j) - dth*s(i,j  ,k,L)*divu(i,j,  k)
                  end if

               else
                  
                  stylo(j) = ylo(i,j,k) &
                      - dt4x*(uedge(i+1,j-1,k  )+uedge(i,j-1,k))* &
                      (xzlo(i+1,j-1,k  )-xzlo(i,j-1,k)) &
                      - dt4z*(wedge(i  ,j-1,k+1)+wedge(i,j-1,k))* &
                      (zxlo(i  ,j-1,k+1)-zxlo(i,j-1,k))
                  styhi(j) = yhi(i,j,k) &
                      - dt4x*(uedge(i+1,j  ,k  )+uedge(i,j  ,k))* &
                      (xzlo(i+1,j  ,k  )-xzlo(i,j  ,k)) &
                      - dt4z*(wedge(i  ,j  ,k+1)+wedge(i,j  ,k))* &
                      (zxlo(i  ,j  ,k+1)-zxlo(i,j  ,k))

               endif

               if (use_minion.eq.0) then
                  stylo(j) = stylo(j) + dth*tf(i,j-1,k,L)
                  styhi(j) = styhi(j) + dth*tf(i,j,  k,L)
               end if
               
            end do

            if (bc(2,1,L).eq.EXT_DIR .and. vedge(i,jmin,k).ge.zero) then
               styhi(jmin) = s(i,jmin-1,k,L)
               stylo(jmin) = s(i,jmin-1,k,L)
            else if (bc(2,1,L).eq.EXT_DIR .and. vedge(i,jmin,k).lt.zero) then
               stylo(jmin) = styhi(jmin)
            else if (bc(2,1,L).eq.FOEXTRAP.or.bc(2,1,L).eq.HOEXTRAP) then
               if (n.eq.YVEL) then
                  if (vedge(i,jmin,k).ge.zero) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                     styhi(jmin) = MIN(styhi(jmin),zero)
#endif
                     stylo(jmin) = styhi(jmin)
                  endif
               else
                  stylo(jmin) = styhi(jmin)
               endif
            else if (bc(2,1,L).eq.REFLECT_EVEN) then
               stylo(jmin) = styhi(jmin)
            else if (bc(2,1,L).eq.REFLECT_ODD) then
               styhi(jmin) = zero
               stylo(jmin) = zero
            end if
            
            if (bc(2,2,L).eq.EXT_DIR .and. vedge(i,jmax+1,k).le.zero) then
               stylo(jmax+1) = s(i,jmax+1,k,L)
               styhi(jmax+1) = s(i,jmax+1,k,L)
            else if (bc(2,2,L).eq.EXT_DIR .and. vedge(i,jmax+1,k).le.zero) then
               styhi(jmax+1) = stylo(jmax+1)
            else if (bc(2,2,L).eq.FOEXTRAP.or.bc(2,2,L).eq.HOEXTRAP) then
               if (n.eq.YVEL) then
                  if (vedge(i,jmax+1,k).le.zero) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                     stylo(jmax+1) = MAX(stylo(jmax+1),zero)
#endif
                     styhi(jmax+1) = stylo(jmax+1)
                  endif
               else
                  styhi(jmax+1) = stylo(jmax+1)
               endif
            else if (bc(2,2,L).eq.REFLECT_EVEN) then
               styhi(jmax+1) = stylo(jmax+1)
            else if (bc(2,2,L).eq.REFLECT_ODD) then
               stylo(jmax+1) = zero
               styhi(jmax+1) = zero
            end if
            
            do j=jmin,jmax+1
               ystate(i,j,k,L) = merge(stylo(j),styhi(j),vedge(i,j,k) .ge. 0.0d0)
               ystate(i,j,k,L) = merge(half*(stylo(j)+styhi(j)),ystate(i,j,k,L), &
                   abs(vedge(i,j,k)).lt.eps)
            end do
         end do
      end do
!c     
!c     compute the zedge states
!c
      do j = jmin,jmax
         do i = imin,imax
            do k = kmin,kmax+1
                  
               if (iconserv(L).eq.1) then
                 
                  stzlo(k) = zlo(i,j,k) &
                      - dthx*(xylo(i+1,j  ,k-1)*uedge(i+1,j  ,k-1) &
                      - xylo(i,j,k-1)*uedge(i,j,k-1)) &
                      - dthy*(yxlo(i  ,j+1,k-1)*vedge(i  ,j+1,k-1) &
                      - yxlo(i,j,k-1)*vedge(i,j,k-1)) &
                      + dthx*s(i,j,k-1,L)*(uedge(i+1,j,k-1)-uedge(i,j,k-1)) &
                      + dthy*s(i,j,k-1,L)*(vedge(i,j+1,k-1)-vedge(i,j,k-1))
                  stzhi(k) = zhi(i,j,k) &
                      - dthx*(xylo(i+1,j  ,k  )*uedge(i+1,j  ,k  ) &
                      - xylo(i,j,k  )*uedge(i,j,k  )) &
                      - dthy*(yxlo(i  ,j+1,k  )*vedge(i  ,j+1,k  ) &
                      - yxlo(i,j,k  )*vedge(i,j,k  )) &
                      + dthx*s(i,j,k,L)*(uedge(i+1,j,k)-uedge(i,j,k)) &
                      + dthy*s(i,j,k,L)*(vedge(i,j+1,k)-vedge(i,j,k))

                  if (use_minion.eq.0) then
                     stzlo(k) = stzlo(k) - dth*s(i,j,k-1,L)*divu(i,j,k-1)
                     stzhi(k) = stzhi(k) - dth*s(i,j,k  ,L)*divu(i,j,k  )
                  end if

               else

                  stzlo(k) = zlo(i,j,k) &
                      - dt4x*(uedge(i+1,j  ,k-1)+uedge(i,j,k-1)) &
                      *(xylo(i+1,j  ,k-1)-xylo(i,j,k-1)) &
                      - dt4y*(vedge(i  ,j+1,k-1)+vedge(i,j,k-1)) &
                      *(yxlo(i  ,j+1,k-1)-yxlo(i,j,k-1))
                  
                  stzhi(k) = zhi(i,j,k) &
                      - dt4x*(uedge(i+1,j  ,k  )+uedge(i,j,k  )) &
                      *(xylo(i+1,j  ,k  )-xylo(i,j,k  )) &
                      - dt4y*(vedge(i  ,j+1,k  )+vedge(i,j,k  )) &
                      *(yxlo(i  ,j+1,k  )-yxlo(i,j,k  ))

               endif

               if (use_minion.eq.0) then
                  stzlo(k) = stzlo(k) + dth*tf(i,j,k-1,L)
                  stzhi(k) = stzhi(k) + dth*tf(i,j,k,L)
               end if

            end do

            if (bc(3,1,L).eq.EXT_DIR .and. wedge(i,j,kmin).ge.zero) then
               stzlo(kmin) = s(i,j,kmin-1,L)
               stzhi(kmin) = s(i,j,kmin-1,L)
            else if (bc(3,1,L).eq.EXT_DIR .and. wedge(i,j,kmin).lt.zero) then
               stzlo(kmin) = stzhi(kmin)
            else if (bc(3,1,L).eq.FOEXTRAP.or.bc(3,1,L).eq.HOEXTRAP) then
               if (n.eq.ZVEL) then
                  if (wedge(i,j,kmin).ge.zero) then
#ifndef ALLOWZINFLOW
!c     prevent backflow
                     stzhi(kmin) = MIN(stzhi(kmin),zero)
#endif
                     stzlo(kmin) = stzhi(kmin)
                  endif
               else
                  stzlo(kmin) = stzhi(kmin)
               endif
            else if (bc(3,1,L).eq.REFLECT_EVEN) then
               stzlo(kmin) = stzhi(kmin)
            else if (bc(3,1,L).eq.REFLECT_ODD) then
               stzlo(kmin) = zero
               stzhi(kmin) = zero
            end if
            if (bc(3,2,L).eq.EXT_DIR .and. wedge(i,j,kmax+1).le.zero) then
               stzlo(kmax+1) = s(i,j,kmax+1,L)
               stzhi(kmax+1) = s(i,j,kmax+1,L)
            else if (bc(3,2,L).eq.EXT_DIR .and. wedge(i,j,kmax+1).gt.zero) then
               stzhi(kmax+1) = stzlo(kmax+1)
            else if (bc(3,2,L).eq.FOEXTRAP.or.bc(3,2,L).eq.HOEXTRAP) then
               if (n.eq.ZVEL) then
                  if (wedge(i,j,kmax+1).le.zero) then
#ifndef ALLOWZINFLOW
!c     prevent backflow
                     stzlo(kmax+1) = MAX(stzlo(kmax+1),zero)
#endif
                     stzhi(kmax+1) = stzlo(kmax+1)
                  endif
               else
                  stzhi(kmax+1) = stzlo(kmax+1)
               endif
            else if (bc(3,2,L).eq.REFLECT_EVEN) then
               stzhi(kmax+1) = stzlo(kmax+1)
            else if (bc(3,2,L).eq.REFLECT_ODD) then
               stzlo(kmax+1) = zero
               stzhi(kmax+1) = zero
            end if
               
            do k = kmin,kmax+1
               zstate(i,j,k,L) = merge(stzlo(k),stzhi(k),wedge(i,j,k) .ge. 0.0d0)
               zstate(i,j,k,L) = merge(half*(stzlo(k)+stzhi(k)),zstate(i,j,k,L), &
                   abs(wedge(i,j,k)).lt.eps)
            end do
         end do
      end do
      
      else
!c    
!c     ORIGINAL NON-CORNER COUPLING CODE
!c
!c
!c     compute the xedge states
!c

      do k = kmin,kmax
            do j = jmin,jmax
               do i = imin-1,imax+1
                  if (iconserv(L).eq.1) then
                     tr = &
                         (vedge(i,j+1,k)*yedge(i,j+1,k) - vedge(i,j,k)*yedge(i,j,k))*ihy +  & 
                         (wedge(i,j,k+1)*zedge(i,j,k+1) - wedge(i,j,k)*zedge(i,j,k))*ihz   
                     st = -dth*tr + dth*(tf(i,j,k,L) - s(i,j,k,L)*divu(i,j,k)) &
                         + dth*s(i,j,k,L)*(vedge(i,j+1,k)-vedge(i,j,k))*ihy &
                         + dth*s(i,j,k,L)*(wedge(i,j,k+1)-wedge(i,j,k))*ihz
                  else
                     if (vedge(i,j,k)*vedge(i,j+1,k).le.0.d0) then
                        vbar = 0.5d0*(vedge(i,j,k)+vedge(i,j+1,k))
                        if (vbar.lt.0.d0) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr1 = vbar*(s(i,j+inc,k,L)-s(i,j+inc-1,k,L))*ihy
                     else
                        tr1 = half*(vedge(i,j+1,k) + vedge(i,j,k)) * &
                                    (yedge(i,j+1,k) -   yedge(i,j,k)  ) *ihy
                     endif
                     if (wedge(i,j,k)*wedge(i,j,k+1).lt.0.d0) then
                        wbar = 0.5d0*(wedge(i,j,k)+wedge(i,j,k+1))
                        if (wbar.lt.0.d0) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr2 = wbar*(s(i,j,k+inc,L)-s(i,j,k+inc-1,L))*ihz
                     else
                        tr2 = half*(wedge(i,j,k+1) + wedge(i,j,k)) * &
                                    (zedge(i,j,k+1) -   zedge(i,j,k)  ) *ihz
                     endif

                     st = -dth*(tr1 + tr2) + dth*tf(i,j,k,L)
                  endif

                  if (ppm_type .gt. 0) then
                     stxlo(i+1)= Ipx(i,j,k) + st
                     stxhi(i  )= Imx(i,j,k) + st
                  else
                     stxlo(i+1)= s(i,j,k,L) + (half-dthx*uedge(i+1,j,k))*sx(i,j,k) + st
                     stxhi(i  )= s(i,j,k,L) - (half+dthx*uedge(i  ,j,k))*sx(i,j,k) + st
                  end if

               end do

               if (bc(1,1,L).eq.EXT_DIR .and. uedge(imin,j,k).ge.zero) then
                  stxhi(imin) = s(imin-1,j,k,L)
                  stxlo(imin) = s(imin-1,j,k,L)
               else if (bc(1,1,L).eq.EXT_DIR .and. uedge(imin,j,k).lt.zero) then
                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1,L).eq.FOEXTRAP.or.bc(1,1,L).eq.HOEXTRAP) then
                  if (n.eq.XVEL) then
                     if (uedge(imin,j,k).ge.zero) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                        stxhi(imin) = MIN(stxhi(imin),zero)
#endif
                        stxlo(imin) = stxhi(imin)
                     endif
                  else
                     stxlo(imin) = stxhi(imin)
                  endif
               else if (bc(1,1,L).eq.REFLECT_EVEN) then
                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1,L).eq.REFLECT_ODD) then
                  stxhi(imin) = zero
                  stxlo(imin) = zero
               end if
               if (bc(1,2,L).eq.EXT_DIR .and. uedge(imax+1,j,k).le.zero) then
                  stxlo(imax+1) = s(imax+1,j,k,L)
                  stxhi(imax+1) = s(imax+1,j,k,L)
               else if (bc(1,2,L).eq.EXT_DIR .and. uedge(imax+1,j,k).gt.zero) then
                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2,L).eq.FOEXTRAP.or.bc(1,2,L).eq.HOEXTRAP) then
                  if (n.eq.XVEL) then
                     if (uedge(imax+1,j,k).le.zero) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                        stxlo(imax+1) = MAX(stxlo(imax+1),zero)
#endif
                        stxhi(imax+1) = stxlo(imax+1)
                     endif
                  else
                     stxhi(imax+1) = stxlo(imax+1)
                  endif
               else if (bc(1,2,L).eq.REFLECT_EVEN) then
                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2,L).eq.REFLECT_ODD) then
                  stxlo(imax+1) = zero
                  stxhi(imax+1) = zero
               end if

               do i = imin, imax+1
                  xstate(i,j,k,L) = merge(stxlo(i),stxhi(i),uedge(i,j,k) .ge. 0.0d0)
                  xstate(i,j,k,L) = merge(half*(stxlo(i)+stxhi(i)),xstate(i,j,k,L) &
                      ,abs(uedge(i,j,k)).lt.eps)
               end do
            end do
      end do

!c
!c     compute the yedge states
!c

      do k = kmin,kmax
            do i = imin,imax
               do j = jmin-1,jmax+1

                  if (iconserv(L).eq.1) then

                     tr = &
                         (uedge(i+1,j,k)*xedge(i+1,j,k) - uedge(i,j,k)*xedge(i,j,k))*ihx +    &
                         (wedge(i,j,k+1)*zedge(i,j,k+1) - wedge(i,j,k)*zedge(i,j,k))*ihz   

                     st = -dth*tr + dth*(tf(i,j,k,L) - s(i,j,k,L)*divu(i,j,k)) &
                         + dth*s(i,j,k,L)*(uedge(i+1,j,k)-uedge(i,j,k))*ihx &
                         + dth*s(i,j,k,L)*(wedge(i,j,k+1)-wedge(i,j,k))*ihz
                  else
                     if (uedge(i,j,k)*uedge(i+1,j,k).le.0.d0) then
                        ubar = 0.5d0*(uedge(i,j,k)+uedge(i+1,j,k))
                        if (ubar.lt.0.d0) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr1 = ubar*(s(i+inc,j,k,L)-s(i+inc-1,j,k,L))*ihx
                     else
                        tr1 = half*(uedge(i+1,j,k) + uedge(i,j,k)) * &
                                    (xedge(i+1,j,k) -   xedge(i,j,k)  ) *ihx
                     endif
                     if (wedge(i,j,k)*wedge(i,j,k+1).lt.0.d0) then
                        wbar = 0.5d0*(wedge(i,j,k)+wedge(i,j,k+1))
                        if (wbar.lt.0.d0) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr2 = wbar*(s(i,j,k+inc,L)-s(i,j,k+inc-1,L))*ihz
                     else
                        tr2 = half*(wedge(i,j,k+1) + wedge(i,j,k)) * &
                                    (zedge(i,j,k+1) -   zedge(i,j,k)  ) *ihz
                     endif

                     st = -dth*(tr1 + tr2) + dth*tf(i,j,k,L)
                  endif

                  if (ppm_type .gt. 0) then
                     stylo(j+1)= Ipy(i,j,k) + st
                     styhi(j  )= Imy(i,j,k) + st
                  else
                     stylo(j+1)= s(i,j,k,L) + (half-dthy*vedge(i,j+1,k))*sy(i,j,k) + st
                     styhi(j  )= s(i,j,k,L) - (half+dthy*vedge(i,j  ,k))*sy(i,j,k) + st
                  end if

               end do

               if (bc(2,1,L).eq.EXT_DIR .and. vedge(i,jmin,k).ge.zero) then
                  styhi(jmin) = s(i,jmin-1,k,L)
                  stylo(jmin) = s(i,jmin-1,k,L)
               else if (bc(2,1,L).eq.EXT_DIR .and. vedge(i,jmin,k).lt.zero) then
                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1,L).eq.FOEXTRAP.or.bc(2,1,L).eq.HOEXTRAP) then
                  if (n.eq.YVEL) then
                     if (vedge(i,jmin,k).ge.zero) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                        styhi(jmin) = MIN(styhi(jmin),zero)
#endif
                        stylo(jmin) = styhi(jmin)
                     endif
                  else
                     stylo(jmin) = styhi(jmin)
                  endif
               else if (bc(2,1,L).eq.REFLECT_EVEN) then
                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1,L).eq.REFLECT_ODD) then
                  styhi(jmin) = zero
                  stylo(jmin) = zero
               end if
               
               if (bc(2,2,L).eq.EXT_DIR .and. vedge(i,jmax+1,k).le.zero) then
                  stylo(jmax+1) = s(i,jmax+1,k,L)
                  styhi(jmax+1) = s(i,jmax+1,k,L)
               else if (bc(2,2,L).eq.EXT_DIR .and. vedge(i,jmax+1,k).le.zero) then
                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2,L).eq.FOEXTRAP.or.bc(2,2,L).eq.HOEXTRAP) then
                  if (n.eq.YVEL) then
                     if (vedge(i,jmax+1,k).le.zero) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                        stylo(jmax+1) = MAX(stylo(jmax+1),zero)
#endif
                        styhi(jmax+1) = stylo(jmax+1)
                     endif
                  else
                     styhi(jmax+1) = stylo(jmax+1)
                  endif
               else if (bc(2,2,L).eq.REFLECT_EVEN) then
                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2,L).eq.REFLECT_ODD) then
                  stylo(jmax+1) = zero
                  styhi(jmax+1) = zero
               end if

               do j=jmin,jmax+1
                  ystate(i,j,k,L) = merge(stylo(j),styhi(j),vedge(i,j,k) .ge. 0.0d0)
                  ystate(i,j,k,L) = merge(half*(stylo(j)+styhi(j)),ystate(i,j,k,L), &
                      abs(vedge(i,j,k)).lt.eps)
               end do
            end do
      end do

!c
!c     compute the zedge states
!c

      do j = jmin,jmax
            do i = imin,imax
               do k = kmin-1,kmax+1

                  if (iconserv(L).eq.1) then
                     tr = &
                         (uedge(i+1,j,k)*xedge(i+1,j,k) - uedge(i,j,k)*xedge(i,j,k))*ihx +  &  
                         (vedge(i,j+1,k)*yedge(i,j+1,k) - vedge(i,j,k)*yedge(i,j,k))*ihy   
                     
                     st = -dth*tr + dth*(tf(i,j,k,L) - s(i,j,k,L)*divu(i,j,k)) &
                         + dth*s(i,j,k,L)*(uedge(i+1,j,k)-uedge(i,j,k))*ihx &
                         + dth*s(i,j,k,L)*(vedge(i,j+1,k)-vedge(i,j,k))*ihy
                  else
                     if (uedge(i,j,k)*uedge(i+1,j,k).le.0.d0) then
                        ubar = 0.5d0*(uedge(i,j,k)+uedge(i+1,j,k))
                        if (ubar.lt.0.d0) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr1 = ubar*(s(i+inc,j,k,L)-s(i+inc-1,j,k,L))*ihx
                     else
                        tr1 = half*(uedge(i+1,j,k) + uedge(i,j,k)) * &
                            (xedge(i+1,j,k) - xedge(i,j,k)  ) *ihx
                     endif
                     if (vedge(i,j,k)*vedge(i,j+1,k).lt.0.d0) then
                        vbar = 0.5d0*(vedge(i,j,k)+vedge(i,j+1,k))
                        if (vbar.lt.0.d0) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr2 = vbar*(s(i,j+inc,k,L)-s(i,j+inc-1,k,L))*ihy
                     else
                        tr2 = half*(vedge(i,j+1,k) + vedge(i,j,k)) * &
                            (yedge(i,j+1,k) - yedge(i,j,k)  ) *ihy
                     endif

                     st = -dth*(tr1 + tr2) + dth*tf(i,j,k,L)
                  endif

                  if (ppm_type .gt. 0) then
                     stzlo(k+1)= Ipz(i,j,k) + st
                     stzhi(k  )= Imz(i,j,k) + st
                  else
                     stzlo(k+1)= s(i,j,k,L) + (half-dthz*wedge(i,j,k+1))*sz(i,j,k) + st
                     stzhi(k  )= s(i,j,k,L) - (half+dthz*wedge(i,j,k  ))*sz(i,j,k) + st
                  end if

               end do

               if (bc(3,1,L).eq.EXT_DIR .and. wedge(i,j,kmin).ge.zero) then
                  stzlo(kmin) = s(i,j,kmin-1,L)
                  stzhi(kmin) = s(i,j,kmin-1,L)
               else if (bc(3,1,L).eq.EXT_DIR .and. wedge(i,j,kmin).lt.zero) then
                  stzlo(kmin) = stzhi(kmin)
               else if (bc(3,1,L).eq.FOEXTRAP.or.bc(3,1,L).eq.HOEXTRAP) then
                  if (n.eq.ZVEL) then
                     if (wedge(i,j,kmin).ge.zero) then
#ifndef ALLOWZINFLOW
!c     prevent backflow
                        stzhi(kmin) = MIN(stzhi(kmin),zero)
#endif
                        stzlo(kmin) = stzhi(kmin)
                     endif
                  else
                     stzlo(kmin) = stzhi(kmin)
                  endif
               else if (bc(3,1,L).eq.REFLECT_EVEN) then
                  stzlo(kmin) = stzhi(kmin)
               else if (bc(3,1,L).eq.REFLECT_ODD) then
                  stzlo(kmin) = zero
                  stzhi(kmin) = zero
               end if
               if (bc(3,2,L).eq.EXT_DIR .and. wedge(i,j,kmax+1).le.zero) then
                  stzlo(kmax+1) = s(i,j,kmax+1,L)
                  stzhi(kmax+1) = s(i,j,kmax+1,L)
               else if (bc(3,2,L).eq.EXT_DIR .and. wedge(i,j,kmax+1).gt.zero) then
                  stzhi(kmax+1) = stzlo(kmax+1)
               else if (bc(3,2,L).eq.FOEXTRAP.or.bc(3,2,L).eq.HOEXTRAP) then
                  if (n.eq.ZVEL) then
                     if (wedge(i,j,kmax+1).le.zero) then
#ifndef ALLOWZINFLOW
!c     prevent backflow
                        stzlo(kmax+1) = MAX(stzlo(kmax+1),zero)
#endif
                        stzhi(kmax+1) = stzlo(kmax+1)
                     endif
                  else
                     stzhi(kmax+1) = stzlo(kmax+1)
                  endif
               else if (bc(3,2,L).eq.REFLECT_EVEN) then
                  stzhi(kmax+1) = stzlo(kmax+1)
               else if (bc(3,2,L).eq.REFLECT_ODD) then
                  stzlo(kmax+1) = zero
                  stzhi(kmax+1) = zero
               end if

               do k = kmin,kmax+1
                  zstate(i,j,k,L) = merge(stzlo(k),stzhi(k),wedge(i,j,k) .ge. 0.0d0)
                  zstate(i,j,k,L) = merge(half*(stzlo(k)+stzhi(k)),zstate(i,j,k,L), &
                      abs(wedge(i,j,k)).lt.eps)
               end do
            end do
      end do
      
      end if

      end do
      
      end subroutine estate_fpu

      subroutine trans_xbc(lo,hi,&
         s,s_lo,s_hi,&
         xlo,xlo_lo,xlo_hi,&
         xhi,xhi_lo,xhi_hi,&
         uad,uad_lo,uad_hi,&
         n, xbc, eps,ycouple,zcouple)
!c
!c     This subroutine processes boundary conditions on information
!c     traced to cell faces in the x direction.  This is used for
!c     computing velocities and edge states used in calculating
!c     transverse derivatives
!c
      implicit none
      
      integer, intent(in) :: n,xbc(SDIM,2)
      integer, dimension(3), intent(in) :: &
           s_lo,s_hi,xlo_lo,xlo_hi,xhi_lo,xhi_hi,uad_lo,uad_hi,lo,hi
           
      real(rt), intent(in)    :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
      real(rt), intent(inout) :: xlo(xlo_lo(1):xlo_hi(1),xlo_lo(2):xlo_hi(2),xlo_lo(3):xlo_hi(3))
      real(rt), intent(inout) :: xhi(xhi_lo(1):xhi_hi(1),xhi_lo(2):xhi_hi(2),xhi_lo(3):xhi_hi(3))
      real(rt), intent(in)    :: uad(uad_lo(1):uad_hi(1),uad_lo(2):uad_hi(2),uad_lo(3):uad_hi(3))
      real(rt), intent(in)    :: eps
      
      logical ycouple
      logical zcouple

      real(rt) ::   stx
      logical ltest
      integer j,k
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

!c     if we are applying bc's to intermediate corner-copuled terms
!c     the bounds are slightly different on the valid data
      if (ycouple) then
         jmin = jmin + 1
         jmax = jmax - 1
      end if

      if (zcouple) then
         kmin = kmin + 1
         kmax = kmax - 1
      end if
!c
!c     -------------- the lower x boundary
!c
      if (xbc(1,1).eq.EXT_DIR) then
         if ( n .eq. XVEL ) then
            do j = jmin-1,jmax+1
             do k = kmin-1,kmax+1
              if (uad(imin,j,k) .ge. 0.0d0) then
                  xlo(imin,j,k) = s(imin-1,j,k)
                  xhi(imin,j,k) = s(imin-1,j,k)
              else
                  xlo(imin,j,k) = xhi(imin,j,k)
              endif
             end do
            end do
         else
            do j = jmin-1,jmax+1
               do k = kmin-1,kmax+1
                  ltest = uad(imin,j,k).le.eps
                  stx   = merge(xhi(imin,j,k),s(imin-1,j,k),ltest)
                  xlo(imin,j,k) = stx
                  xhi(imin,j,k) = stx
               end do
            end do
         end if
      else if (xbc(1,1).eq.FOEXTRAP.or.xbc(1,1).eq.HOEXTRAP &
             .or.xbc(1,1).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            do k = kmin-1,kmax+1
               xlo(imin,j,k) = xhi(imin,j,k)
            end do
         end do
      else if (xbc(1,1).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            do k = kmin-1,kmax+1
               xhi(imin,j,k) = 0.0d0
               xlo(imin,j,k) = 0.0d0
            end do
         end do
      end if
!c
!c     -------------- the upper x boundary
!c
      if (xbc(1,2).eq.EXT_DIR) then
         if ( n .eq. XVEL ) then
            do j = jmin-1,jmax+1
             do k = kmin-1,kmax+1
               if (uad(imax+1,j,k) .le. 0.0d0) then
                  xlo(imax+1,j,k) = s(imax+1,j,k)
                  xhi(imax+1,j,k) = s(imax+1,j,k)
               else
                  xhi(imax+1,j,k) = xlo(imax+1,j,k)
               endif
             end do
            end do
         else
            do j = jmin-1,jmax+1
               do k = kmin-1,kmax+1
                  ltest = uad(imax+1,j,k).ge.-eps
                  stx   = merge(xlo(imax+1,j,k),s(imax+1,j,k),ltest)
                  xlo(imax+1,j,k) = stx
                  xhi(imax+1,j,k) = stx
               end do
            end do
         end if
      else if (xbc(1,2).eq.FOEXTRAP.or.xbc(1,2).eq.HOEXTRAP &
             .or.xbc(1,2).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            do k = kmin-1,kmax+1
               xhi(imax+1,j,k) = xlo(imax+1,j,k)
            end do
         end do
      else if (xbc(1,2).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            do k = kmin-1,kmax+1
               xhi(imax+1,j,k) = 0.0d0
               xlo(imax+1,j,k) = 0.0d0
            end do
         end do
      end if

      end subroutine trans_xbc

      subroutine trans_ybc( lo,hi,&
         s,s_lo,s_hi,&
         ylo,ylo_lo,ylo_hi,&
         yhi,yhi_lo,yhi_hi,&
         vad,vad_lo,vad_hi,&
         n, ybc, eps,xcouple,zcouple)
!c
!c     This subroutine processes boundary conditions on information
!c     traced to cell faces in the y direction.  This is used for
!c     computing velocities and edge states used in calculating
!c     transverse derivatives
!c

      implicit none
      
      integer, intent(in) :: n,ybc(SDIM,2)
      integer, dimension(3), intent(in) :: &
           s_lo,s_hi,ylo_lo,ylo_hi,yhi_lo,yhi_hi,vad_lo,vad_hi,lo,hi
           
      real(rt), intent(in)    :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
      real(rt), intent(inout) :: ylo(ylo_lo(1):ylo_hi(1),ylo_lo(2):ylo_hi(2),ylo_lo(3):ylo_hi(3))
      real(rt), intent(inout) :: yhi(yhi_lo(1):yhi_hi(1),yhi_lo(2):yhi_hi(2),yhi_lo(3):yhi_hi(3))
      real(rt), intent(in)    :: vad(vad_lo(1):vad_hi(1),vad_lo(2):vad_hi(2),vad_lo(3):vad_hi(3))
      real(rt), intent(in)    ::  eps
      
      logical xcouple
      logical zcouple

      real(rt) ::   sty
      logical ltest
      integer i,k
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

!c     if we are applying bc's to intermediate corner-copuled terms
!c     the bounds are slightly different on the valid data
      if (xcouple) then
         imin = imin + 1
         imax = imax - 1
      end if

      if (zcouple) then
         kmin = kmin + 1
         kmax = kmax - 1
      end if
!c
!c     -------------- the lower y boundary
!c
      if (ybc(2,1).eq.EXT_DIR) then
         if ( n .eq. YVEL ) then
            do i = imin-1,imax+1
             do k = kmin-1,kmax+1
              if (vad(i,jmin,k) .ge. zero) then
                  ylo(i,jmin,k) = s(i,jmin-1,k)
                  yhi(i,jmin,k) = s(i,jmin-1,k)
              else
                  ylo(i,jmin,k) = yhi(i,jmin,k)
              endif
             end do
            end do
         else
            do i = imin-1,imax+1
               do k = kmin-1,kmax+1
                  ltest = vad(i,jmin,k).le.eps
                  sty   = merge(yhi(i,jmin,k),s(i,jmin-1,k),ltest)
                  ylo(i,jmin,k) = sty
                  yhi(i,jmin,k) = sty
               end do
            end do
         end if
      else if (ybc(2,1).eq.FOEXTRAP.or.ybc(2,1).eq.HOEXTRAP &
             .or.ybc(2,1).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            do k = kmin-1,kmax+1
               ylo(i,jmin,k) = yhi(i,jmin,k)
            end do
         end do
      else if (ybc(2,1).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            do k = kmin-1,kmax+1
               yhi(i,jmin,k) = zero
               ylo(i,jmin,k) = zero
            end do
         end do
      end if
!c
!c     -------------- the upper y boundary
!c
      if (ybc(2,2).eq.EXT_DIR) then
         if ( n .eq. YVEL ) then
            do i = imin-1,imax+1
             do k = kmin-1,kmax+1
               if (vad(i,jmax+1,k) .le. zero) then
                  ylo(i,jmax+1,k) = s(i,jmax+1,k)
                  yhi(i,jmax+1,k) = s(i,jmax+1,k)
               else
                  yhi(i,jmax+1,k) = ylo(i,jmax+1,k)
               endif
             end do
            end do
         else
            do i = imin-1,imax+1
               do k = kmin-1,kmax+1
                  ltest = vad(i,jmax+1,k).ge.-eps
                  sty   = merge(ylo(i,jmax+1,k),s(i,jmax+1,k),ltest)
                  ylo(i,jmax+1,k) = sty
                  yhi(i,jmax+1,k) = sty
               end do
            end do
         end if
      else if (ybc(2,2).eq.FOEXTRAP.or.ybc(2,2).eq.HOEXTRAP &
             .or.ybc(2,2).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            do k = kmin-1,kmax+1
               yhi(i,jmax+1,k) = ylo(i,jmax+1,k)
            end do
         end do
      else if (ybc(2,2).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            do k = kmin-1,kmax+1
               ylo(i,jmax+1,k) = zero
               yhi(i,jmax+1,k) = zero
            end do
         end do
      end if

      end subroutine trans_ybc

      subroutine trans_zbc(lo,hi,&
         s,s_lo,s_hi,&
         zlo,zlo_lo,zlo_hi,&
         zhi,zhi_lo,zhi_hi,&
         wad,wad_lo,wad_hi,&
         n, zbc, eps,xcouple,ycouple)
!c
!c     This subroutine processes boundary conditions on information
!c     traced to cell faces in the z direction.  This is used for
!c     computing velocities and edge states used in calculating
!c     transverse derivatives
!c
      implicit none
      integer, intent(in) :: n,zbc(SDIM,2)
      integer, dimension(3), intent(in) :: &
           s_lo,s_hi,zlo_lo,zlo_hi,zhi_lo,zhi_hi,wad_lo,wad_hi,lo,hi

      real(rt), intent(in)    :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
      real(rt), intent(inout) :: zlo(zlo_lo(1):zlo_hi(1),zlo_lo(2):zlo_hi(2),zlo_lo(3):zlo_hi(3))
      real(rt), intent(inout) :: zhi(zhi_lo(1):zhi_hi(1),zhi_lo(2):zhi_hi(2),zhi_lo(3):zhi_hi(3))
      real(rt), intent(in)    :: wad(wad_lo(1):wad_hi(1),wad_lo(2):wad_hi(2),wad_lo(3):wad_hi(3))
      real(rt), intent(in)    ::  eps

      logical xcouple
      logical ycouple
      
      real(rt) :: stz
      logical ltest
      integer i,j
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

!c     if we are applying bc's to intermediate corner-copuled terms
!c     the bounds are slightly different on the valid data
      if (xcouple) then
         imin = imin + 1
         imax = imax - 1
      end if

      if (ycouple) then
         jmin = jmin + 1
         jmax = jmax - 1
      end if
!c
!c     -------------- the lower z boundary
!c
      if (zbc(3,1).eq.EXT_DIR) then
         if ( n .eq. ZVEL ) then
            do i = imin-1,imax+1
             do j = jmin-1,jmax+1
               if (wad(i,j,kmin) .ge. zero) then
                  zhi(i,j,kmin) = s(i,j,kmin-1)
                  zlo(i,j,kmin) = s(i,j,kmin-1)
               else
                  zlo(i,j,kmin) = zhi(i,j,kmin)
               endif
             end do
            end do
         else
            do i = imin-1,imax+1
               do j = jmin-1,jmax+1
                  ltest = wad(i,j,kmin).le.eps
                  stz   = merge(zhi(i,j,kmin),s(i,j,kmin-1),ltest)
                  zhi(i,j,kmin) = stz
                  zlo(i,j,kmin) = stz
               end do
            end do
         end if
      else if (zbc(3,1).eq.FOEXTRAP.or.zbc(3,1).eq.HOEXTRAP &
             .or.zbc(3,1).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            do j = jmin-1,jmax+1
               zlo(i,j,kmin) = zhi(i,j,kmin)
            end do
         end do
      else if (zbc(3,1).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            do j = jmin-1,jmax+1
               zhi(i,j,kmin) = zero
               zlo(i,j,kmin) = zero
            end do
         end do
      end if
!c
!c     -------------- the upper z boundary
!c
      if (zbc(3,2).eq.EXT_DIR) then
         if ( n .eq. ZVEL ) then
            do i = imin-1,imax+1
             do j = jmin-1,jmax+1
               if (wad(i,j,kmax+1) .le. zero) then
                  zlo(i,j,kmax+1) = s(i,j,kmax+1)
                  zhi(i,j,kmax+1) = s(i,j,kmax+1)
               else
                  zhi(i,j,kmax+1) = zlo(i,j,kmax+1)
               endif
             end do
            end do
         else
            do i = imin-1,imax+1
               do j = jmin-1,jmax+1
                  ltest = wad(i,j,kmax+1).ge.-eps
                  stz   = merge(zlo(i,j,kmax+1),s(i,j,kmax+1),ltest)
                  zhi(i,j,kmax+1) = stz
                  zlo(i,j,kmax+1) = stz
               end do
            end do
         end if
      else if (zbc(3,2).eq.FOEXTRAP.or.zbc(3,2).eq.HOEXTRAP &
             .or.zbc(3,2).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            do j = jmin-1,jmax+1
               zhi(i,j,kmax+1) = zlo(i,j,kmax+1)
            end do
         end do
      else if (zbc(3,2).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            do j = jmin-1,jmax+1
               zlo(i,j,kmax+1) = zero
               zhi(i,j,kmax+1) = zero
            end do
         end do
      end if

      end subroutine trans_zbc

      subroutine slopes(dir, lo,hi,&
         s,s_lo,s_hi,&
         slx,slx_lo,slx_hi,&
         sly,sly_lo,sly_hi,&
         slz,slz_lo,slz_hi,&
         bc)
!c 
!c     this subroutine computes first or forth order slopes of
!c     a 3D scalar field.
!c
!c     (dir) is used to eliminate calculating extra slopes in transvel
!c
!c     Boundary conditions on interior slopes are handled automatically
!c     by the ghost cells
!c
!c     Boundary conditions on EXT_DIR and HOEXTRAP slopes are implemented
!c     by setting them to zero outside of the domain and using a
!c     one-sided derivative from the interior
!c
      implicit none

#include <GODCOMM_F.H>

      integer :: dir
      integer, intent(in) :: bc(SDIM,2)
      integer, dimension(3), intent(in) :: lo,hi,s_lo,s_hi,slx_lo,slx_hi,sly_lo,sly_hi,slz_lo,slz_hi
      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
      real(rt), intent(inout) :: slx(slx_lo(1):slx_hi(1),slx_lo(2):slx_hi(2),slx_lo(3):slx_hi(3))
      real(rt), intent(inout) :: sly(sly_lo(1):sly_hi(1),sly_lo(2):sly_hi(2),sly_lo(3):sly_hi(3))
      real(rt), intent(inout) :: slz(slz_lo(1):slz_hi(1),slz_lo(2):slz_hi(2),slz_lo(3):slz_hi(3))
      real(rt) :: slxscr(lo(1)-2:hi(1)+2,4)
      real(rt) :: slyscr(lo(2)-2:hi(2)+2,4)
      real(rt) :: slzscr(lo(3)-2:hi(3)+2,4)
      
      integer imin,jmin,kmin,imax,jmax,kmax,i,j,k
      integer ng
      real(rt) dpls,dmin,ds
      real(rt) del,slim,sflg,sixteen15ths
      integer cen,lim,flag,fromm

      parameter( cen = 1, lim = 2, flag = 3, fromm = 4 )
      parameter( sixteen15ths = sixteen/fifteen )


!C
!C     Determine ng in a way that covers the case of tiling where
!C     (lo:hi) is only a portion of the box s is defined on.
!C
      ng = lo(1) - s_lo(1)
      if (slope_order .eq.1) then
         if (ng .lt. 1) then
            call bl_abort("FORT_SLOPES: too few bndry cells for first order")
         endif
         ng = 1
      else if (slope_order .eq. 2) then
         if (ng .lt. 2) then
            call bl_abort("SLOPE_2D: not enough bndry cells for 2nd order")
         endif
         ng = 2
      else
         if (ng .lt. 3) then
            call bl_abort("SLOPE_2D: not enough bndry cells for 4th order")
         end if
         ng = 3
      endif

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
!c
!c     Added to prevent underflow for small s values.
!c

!      do k = lo(3)-ng, hi(3)+ng
!         do j = lo(2)-ng, hi(2)+ng
!            do i = lo(1)-ng, hi(1)+ng
!               s(i,j,k) = merge(s(i,j,k), zero, abs(s(i,j,k)).gt.1.0D-20)
!            end do
!         end do
!      end do

!c
!c     COMPUTE 0TH order slopes
!c
      if (slope_order.eq.1) then
         do k = kmin-1, kmax+1
            do j = jmin-1, jmax+1 
               do i = imin-1, imax+1
                  slx(i,j,k) = zero
                  sly(i,j,k) = zero
                  slz(i,j,k) = zero
               end do
            end do
         end do
         return
      end if
!c
!c     COMPUTE 2ND order slopes
!c
      if (slope_order.eq.2) then
!c
!c     ------------------------ x slopes
!c
        if ( (dir.eq.XVEL) .or. (dir.eq.ALL) ) then
           if (use_unlimited_slopes) then
              do k = kmin-1,kmax+1
                 do j = jmin-1,jmax+1
                    do i = imin-1,imax+1
                       slx(i,j,k)= half*(s(i+1,j,k) - s(i-1,j,k))
                    end do
                    if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                       slx(imin-1,j,k) = zero
                       slx(imin  ,j,k) = (s(imin+1,j,k)+three*s(imin,j,k)-four*s(imin-1,j,k))*third
                    end if
                    if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                       slx(imax+1,j,k) = zero
                       slx(imax  ,j,k) = -(s(imax-1,j,k)+three*s(imax,j,k)-four*s(imax+1,j,k))*third
                    end if
                 end do
              end do
           else
              do k = kmin-1,kmax+1
                 do j = jmin-1,jmax+1
                    do i = imin-1,imax+1
                       del  = half*(s(i+1,j,k) - s(i-1,j,k))
                       dpls =  two*(s(i+1,j,k) - s(i  ,j,k))
                       dmin =  two*(s(i  ,j,k) - s(i-1,j,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       slx(i,j,k)= sflg*min(slim,abs(del))
                    end do
                    if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                       slx(imin-1,j,k) = zero
                       del  = (s(imin+1,j,k)+three*s(imin,j,k)-four*s(imin-1,j,k))*third
                       dpls = two*(s(imin+1,j,k) - s(imin  ,j,k))
                       dmin = two*(s(imin  ,j,k) - s(imin-1,j,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       slx(imin,j,k)= sflg*min(slim,abs(del))
                    end if
                    if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                       slx(imax+1,j,k) = zero
                       del  = -(s(imax-1,j,k)+three*s(imax,j,k)-four*s(imax+1,j,k))*third
                       dpls = two*(s(imax+1,j,k) - s(imax  ,j,k))
                       dmin = two*(s(imax  ,j,k) - s(imax-1,j,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       slx(imax,j,k)= sflg*min(slim,abs(del))
                    end if
                 end do
              end do
           end if
        end if
!c
!c     ------------------------ y slopes
!c
        if ( (dir.eq.YVEL) .or. (dir.eq.ALL) ) then
           if (use_unlimited_slopes) then
              do k = kmin-1,kmax+1
                 do i = imin-1,imax+1
                    do j = jmin-1,jmax+1
                       sly(i,j,k) = half*(s(i,j+1,k)-s(i,j-1,k))
                    end do
                    if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                       sly(i,jmin-1,k) = zero
                       sly(i,jmin  ,k) = (s(i,jmin+1,k)+three*s(i,jmin,k)-four*s(i,jmin-1,k))*third
                    end if
                    if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                       sly(i,jmax+1,k) = zero
                       sly(i,jmax  ,k) = -(s(i,jmax-1,k)+three*s(i,jmax,k)-four*s(i,jmax+1,k))*third
                    end if
                 end do
              end do
           else
              do k = kmin-1,kmax+1
                 do i = imin-1,imax+1
                    do j = jmin-1,jmax+1
                       del  = half*(s(i,j+1,k) - s(i,j-1,k))
                       dpls =  two*(s(i,j+1,k) - s(i,j  ,k))
                       dmin =  two*(s(i,j  ,k) - s(i,j-1,k))
                       slim = min(abs(dpls),abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       sly(i,j,k)= sflg*min(slim,abs(del))
                    end do
                    if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                       sly(i,jmin-1,k) = zero
                       del  = (s(i,jmin+1,k)+three*s(i,jmin,k)-four*s(i,jmin-1,k))*third
                       dpls = two*(s(i,jmin+1,k) - s(i,jmin  ,k))
                       dmin = two*(s(i,jmin  ,k) - s(i,jmin-1,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       sly(i,jmin,k)= sflg*min(slim,abs(del))
                    end if
                    if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                       sly(i,jmax+1,k) = zero
                       del  = -(s(i,jmax-1,k)+three*s(i,jmax,k)-four*s(i,jmax+1,k))*third
                       dpls = two*(s(i,jmax+1,k) - s(i,jmax  ,k))
                       dmin = two*(s(i,jmax  ,k) - s(i,jmax-1,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       sly(i,jmax,k)= sflg*min(slim,abs(del))
                    end if
                 end do
              end do
           end if
        end if
!c
!c     ------------------------ z slopes
!c
        if ( (dir.eq.ZVEL) .or. (dir.eq.ALL) ) then
           if (use_unlimited_slopes) then
              do j = jmin-1,jmax+1
                 do i = imin-1,imax+1
                    do k = kmin-1,kmax+1
                       slz(i,j,k) = half*(s(i,j,k+1)-s(i,j,k-1))
                    end do
                    if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
                       slz(i,j,kmin-1) = zero
                       slz(i,j,kmin  ) = (s(i,j,kmin+1)+three*s(i,j,kmin)-four*s(i,j,kmin-1))*third
                    end if
                    if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
                       slz(i,j,kmax+1) = zero
                       slz(i,j,kmax  ) = -(s(i,j,kmax-1)+three*s(i,j,kmax)-four*s(i,j,kmax+1))*third
                    end if
                 end do
              end do
           else
              do j = jmin-1,jmax+1
                 do i = imin-1,imax+1
                    do k = kmin-1,kmax+1
                       del  = half*(s(i,j,k+1) - s(i,j,k-1))
                       dpls =  two*(s(i,j,k+1) - s(i,j,k  ))
                       dmin =  two*(s(i,j,k  ) - s(i,j,k-1))
                       slim = min(abs(dpls),abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       slz(i,j,k)= sflg*min(slim,abs(del))
                    end do
                    if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
                       slz(i,j,kmin-1) = zero
                       del  = (s(i,j,kmin+1)+three*s(i,j,kmin)-four*s(i,j,kmin-1))*third
                       dpls = two*(s(i,j,kmin+1) - s(i,j,kmin  ))
                       dmin = two*(s(i,j,kmin  ) - s(i,j,kmin-1))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       slz(i,j,kmin)= sflg*min(slim,abs(del))
                    end if
                    if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
                       slz(i,j,kmax+1) = zero
                       del  = -(s(i,j,kmax-1)+three*s(i,j,kmax)-four*s(i,j,kmax+1))*third
                       dpls = two*(s(i,j,kmax+1) - s(i,j,kmax  ))
                       dmin = two*(s(i,j,kmax  ) - s(i,j,kmax-1))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       slz(i,j,kmax)= sflg*min(slim,abs(del))
                    end if
                 end do
              end do
           end if
        end if
!c
!c ... end, if slope_order .eq. 2
!c
      end if
!c
!c     COMPUTE 4TH order slopes
!c
      if (slope_order.eq.4)then
!c
!c     ------------------------ x slopes
!c
        if ( (dir.eq.XVEL) .or. (dir.eq.ALL) ) then
           if (use_unlimited_slopes) then
              do k = kmin-1,kmax+1
                 do j = jmin-1,jmax+1
                    do i = imin-2,imax+2
                       slxscr(i,cen)  = half*(s(i+1,j,k)-s(i-1,j,k))
                    end do
                    do i = imin-1,imax+1
                       slx(i,j,k) = two * two3rd * slxscr(i,cen) - &
                           sixth * (slxscr(i+1,cen) + slxscr(i-1,cen))
                    end do
                    if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                       slx(imin,j,k) = -sixteen15ths*s(imin-1,j,k) + half*s(imin,j,k) +  &
                           two3rd*s(imin+1,j,k) - tenth*s(imin+2,j,k)
                       slx(imin-1,j,k) = zero
                    end if
                    if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                       slx(imax,j,k) = -( -sixteen15ths*s(imax+1,j,k) + half*s(imax,j,k) +  &
                           two3rd*s(imax-1,j,k) - tenth*s(imax-2,j,k) )
                       slx(imax+1,j,k) = zero
                    end if
                 end do
              end do
           else
              do k = kmin-1,kmax+1
                 do j = jmin-1,jmax+1 
                    do i = imin-2,imax+2
                       dmin           =  two*(s(i,  j,k)-s(i-1,j,k))
                       dpls           =  two*(s(i+1,j,k)-s(i,  j,k))
                       slxscr(i,cen)  = half*(s(i+1,j,k)-s(i-1,j,k))
                       slxscr(i,lim)  = min(abs(dmin),abs(dpls))
                       slxscr(i,lim)  = merge(slxscr(i,lim),zero,(dpls*dmin) .ge. 0.0d0)
                       slxscr(i,flag) = sign(one,slxscr(i,cen))
                       slxscr(i,fromm)= slxscr(i,flag)* &
                           min(slxscr(i,lim),abs(slxscr(i,cen)))
                    end do
                    do i = imin-1,imax+1
                       ds = two * two3rd * slxscr(i,cen) - &
                           sixth * (slxscr(i+1,fromm) + slxscr(i-1,fromm))
                       slx(i,j,k) = slxscr(i,flag)*min(abs(ds),slxscr(i,lim))
                    end do

                    if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                       del  = -sixteen15ths*s(imin-1,j,k) + half*s(imin,j,k) +  &
                           two3rd*s(imin+1,j,k) -  tenth*s(imin+2,j,k)
                       dmin = two*(s(imin  ,j,k)-s(imin-1,j,k))
                       dpls = two*(s(imin+1,j,k)-s(imin  ,j,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       slx(imin-1,j,k) = zero
                       slx(imin,  j,k) = sflg*min(slim,abs(del))

!c                      Recalculate the slope at imin+1 using the revised slxscr(imin,fromm)
                       slxscr(imin,fromm) = slx(imin,j,k)
                       ds = two * two3rd * slxscr(imin+1,cen) - &
                         sixth * (slxscr(imin+2,fromm) + slxscr(imin,fromm))
                       slx(imin+1,j,k) = slxscr(imin+1,flag)*min(abs(ds),slxscr(imin+1,lim))
                    end if

                    if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                       del  = -( -sixteen15ths*s(imax+1,j,k) + half*s(imax,j,k) +  &
                           two3rd*s(imax-1,j,k) - tenth*s(imax-2,j,k) )
                       dmin = two*(s(imax  ,j,k)-s(imax-1,j,k))
                       dpls = two*(s(imax+1,j,k)-s(imax  ,j,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       slx(imax,  j,k) = sflg*min(slim,abs(del))
                       slx(imax+1,j,k) = zero

!c                      Recalculate the slope at imax-1 using the revised slxscr(imax,fromm)
                       slxscr(imax,fromm) = slx(imax,j,k)
                       ds = two * two3rd * slxscr(imax-1,cen) - &
                         sixth * (slxscr(imax-2,fromm) + slxscr(imax,fromm))
                       slx(imax-1,j,k) = slxscr(imax-1,flag)*min(abs(ds),slxscr(imax-1,lim))
                    end if
                 end do
              end do
           end if
        end if
!c
!c     ------------------------ y slopes
!c
        if ( (dir.eq.YVEL) .or. (dir.eq.ALL) ) then
           if (use_unlimited_slopes) then
              do k = kmin-1,kmax+1
                 do i = imin-1,imax+1
                    do j = jmin-2,jmax+2
                       slyscr(j,cen)  = half*(s(i,j+1,k)-s(i,j-1,k))
                    end do
                    do j = jmin-1,jmax+1
                       sly(i,j,k) = two * two3rd * slyscr(j,cen) - &
                           sixth * (slyscr(j+1,cen) + slyscr(j-1,cen))
                    end do
                    if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                       sly(i,jmin-1,k) = zero
                       sly(i,jmin,k) = -sixteen15ths*s(i,jmin-1,k) + half*s(i,jmin,k) + &
                           two3rd*s(i,jmin+1,k) - tenth*s(i,jmin+2,k)
                    end if
                    if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                       sly(i,jmax,k) = -( -sixteen15ths*s(i,jmax+1,k) + half*s(i,jmax,k) + &
                           two3rd*s(i,jmax-1,k) - tenth*s(i,jmax-2,k) )
                       sly(i,jmax+1,k) = zero
                    end if
                 end do
              end do
           else
              do k = kmin-1,kmax+1
                 do i = imin-1,imax+1 
                    do j = jmin-2,jmax+2
                       dmin           =  two*(s(i,j,  k)-s(i,j-1,k))
                       dpls           =  two*(s(i,j+1,k)-s(i,j,  k))
                       slyscr(j,cen)  = half*(s(i,j+1,k)-s(i,j-1,k))
                       slyscr(j,lim)  = min(abs(dmin),abs(dpls))
                       slyscr(j,lim)  = merge(slyscr(j,lim),zero,(dpls*dmin) .ge. 0.0d0)
                       slyscr(j,flag) = sign(one,slyscr(j,cen))
                       slyscr(j,fromm)= slyscr(j,flag)* &
                           min(slyscr(j,lim),abs(slyscr(j,cen)))
                    end do
                    do j = jmin-1,jmax+1
                       ds = two * two3rd * slyscr(j,cen) - &
                           sixth * (slyscr(j+1,fromm) + slyscr(j-1,fromm))
                       sly(i,j,k) = slyscr(j,flag)*min(abs(ds),slyscr(j,lim))
                    end do
!c
                    if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                       del  = -sixteen15ths*s(i,jmin-1,k) + half*s(i,jmin,k) + &
                           two3rd*s(i,jmin+1,k) - tenth*s(i,jmin+2,k)
                       dmin = two*(s(i,jmin  ,k)-s(i,jmin-1,k))
                       dpls = two*(s(i,jmin+1,k)-s(i,jmin  ,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       sly(i,jmin-1,k) = zero
                       sly(i,jmin,  k) = sflg*min(slim,abs(del))

!c                      Recalculate the slope at jmin+1 using the revised slyscr(jmin,fromm)
                       slyscr(jmin,fromm) = sly(i,jmin,k)
                       ds = two * two3rd * slyscr(jmin+1,cen) - &
                         sixth * (slyscr(jmin+2,fromm) + slyscr(jmin,fromm))
                       sly(i,jmin+1,k) = slyscr(jmin+1,flag)*min(abs(ds),slyscr(jmin+1,lim))
                    end if
                    if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                       del  = -( -sixteen15ths*s(i,jmax+1,k) + half*s(i,jmax,k) + &
                           two3rd*s(i,jmax-1,k) - tenth*s(i,jmax-2,k) )
                       dmin = two*(s(i,jmax  ,k)-s(i,jmax-1,k))
                       dpls = two*(s(i,jmax+1,k)-s(i,jmax  ,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       sly(i,jmax, k)  = sflg*min(slim,abs(del))
                       sly(i,jmax+1,k) = zero

!c                      Recalculate the slope at jmax-1 using the revised slyscr(jmax,fromm)
                       slyscr(jmax,fromm) = sly(i,jmax,k)
                       ds = two * two3rd * slyscr(jmax-1,cen) - &
                         sixth * (slyscr(jmax-2,fromm) + slyscr(jmax,fromm))
                       sly(i,jmax-1,k) = slyscr(jmax-1,flag)*min(abs(ds),slyscr(jmax-1,lim))
                    end if
                 end do
              end do
           end if
        end if
!c
!c     ------------------------ z slopes
!c
        if ( (dir.eq.ZVEL) .or. (dir.eq.ALL) ) then
           if (use_unlimited_slopes) then
              do j = jmin-1,jmax+1
                 do i = imin-1,imax+1
                    do k = kmin-2,kmax+2
                       slzscr(k,cen)  = half*(s(i,j,k+1)-s(i,j,k-1))
                    end do
                    do k = kmin-1,kmax+1
                       slz(i,j,k) = two * two3rd * slzscr(k,cen) - &
                           sixth * (slzscr(k+1,cen) + slzscr(k-1,cen))
                    end do
                 end do
                 if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
                    slz(i,j,kmin-1) = zero
                    slz(i,j,kmin) = -sixteen15ths*s(i,j,kmin-1) + half*s(i,j,kmin) + &
                        two3rd*s(i,j,kmin+1) - tenth*s(i,j,kmin+2)
                 end if
                 if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
                    slz(i,j,kmax) = -( -sixteen15ths*s(i,j,kmax+1) + half*s(i,j,kmax) + &
                        two3rd*s(i,j,kmax-1) - tenth*s(i,j,kmax-2) )
                    slz(i,j,kmax+1) = zero
                 end if
              end do
           else
              do j = jmin-1,jmax+1
                 do i = imin-1,imax+1
                    do k = kmin-2,kmax+2
                       dmin           =  two*(s(i,j,k  )-s(i,j,k-1))
                       dpls           =  two*(s(i,j,k+1)-s(i,j,k  ))
                       slzscr(k,cen)  = half*(s(i,j,k+1)-s(i,j,k-1))
                       slzscr(k,lim)  = min(abs(dmin),abs(dpls))
                       slzscr(k,lim)  = merge(slzscr(k,lim),zero,(dpls*dmin) .ge. 0.0d0)
                       slzscr(k,flag) = sign(one,slzscr(k,cen))
                       slzscr(k,fromm)= slzscr(k,flag)* &
                           min(slzscr(k,lim),abs(slzscr(k,cen)))
                    end do
                    do k = kmin-1,kmax+1
                       ds = two * two3rd * slzscr(k,cen) - &
                           sixth * (slzscr(k+1,fromm) + slzscr(k-1,fromm))
                       slz(i,j,k) = slzscr(k,flag)*min(abs(ds),slzscr(k,lim))
                    end do
!c
                    if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
                       del  = -sixteen15ths*s(i,j,kmin-1) + half*s(i,j,kmin) + &
                           two3rd*s(i,j,kmin+1) - tenth*s(i,j,kmin+2)
                       dmin = two*(s(i,j,kmin  )-s(i,j,kmin-1))
                       dpls = two*(s(i,j,kmin+1)-s(i,j,kmin  ))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       slz(i,j,kmin-1) = zero
                       slz(i,j,kmin  ) = sflg*min(slim,abs(del))

!c                      Recalculate the slope at jmin+1 using the revised slzscr(kmin,fromm)
                       slzscr(kmin,fromm) = slz(i,j,kmin)
                       ds = two * two3rd * slzscr(kmin+1,cen) - &
                         sixth * (slzscr(kmin+2,fromm) + slzscr(kmin,fromm))
                       slz(i,j,kmin+1) = slzscr(kmin+1,flag)*min(abs(ds),slzscr(kmin+1,lim))
                    end if
                    if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
                       del  = sixteen15ths*s(i,j,kmax+1) - half*s(i,j,kmax) &
                           - two3rd*s(i,j,kmax-1) + tenth*s(i,j,kmax-2)
                       dmin = two*(s(i,j,kmax  )-s(i,j,kmax-1))
                       dpls = two*(s(i,j,kmax+1)-s(i,j,kmax  ))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                       sflg = sign(one,del)
                       slz(i,j,kmax  ) = sflg*min(slim,abs(del))
                       slz(i,j,kmax+1) = zero

!c                      Recalculate the slope at jmax-1 using the revised slzscr(kmax,fromm)
                       slzscr(kmax,fromm) = slz(i,j,kmax)
                       ds = two * two3rd * slzscr(kmax-1,cen) - &
                         sixth * (slzscr(kmax-2,fromm) + slzscr(kmax,fromm))
                       slz(i,j,kmax-1) = slzscr(kmax-1,flag)*min(abs(ds),slzscr(kmax-1,lim))
                    end if
                 end do
              end do
           end if
        end if
!c
!c ... end, if slope_order .eq. 4
!c
      end if

      end subroutine slopes

      subroutine ppm(lo,hi,&
         s,s_lo,s_hi,&
         u,u_lo,u_hi,&
         v,v_lo,v_hi,&
         w,w_lo,w_hi,&
         Ipx,Ipx_lo,Ipx_hi,&
         Imx,Imx_lo,Imx_hi,&
         Ipy,Ipy_lo,Ipy_hi,&
         Imy,Imy_lo,Imy_hi,&
         Ipz,Ipz_lo,Ipz_hi,&
         Imz,Imz_lo,Imz_hi,&
         sm,sm_lo,sm_hi,&
         sp,sp_lo,sp_hi,&
         dsvl,dsvl_lo,dsvl_hi,&
         sedgex,sedgex_lo,sedgex_hi,&
         sedgey,sedgey_lo,sedgey_hi,&
         sedgez,sedgez_lo,sedgez_hi,&
         dx, dt, bc, eps, ppm_type)

      implicit none

      integer, intent(in) :: bc(SDIM,2)
      integer, dimension(3), intent(in) :: &
           s_lo,s_hi,u_lo,u_hi,v_lo,v_hi,w_lo,w_hi, &
           Ipx_lo,Ipx_hi,Imx_lo,Imx_hi,Ipy_lo,Ipy_hi,Imy_lo,Imy_hi,Ipz_lo,Ipz_hi,Imz_lo,Imz_hi,&
           sm_lo,sm_hi,sp_lo,sp_hi,dsvl_lo,dsvl_hi,sedgex_lo,sedgex_hi,sedgey_lo,sedgey_hi,sedgez_lo,sedgez_hi,&
           lo,hi

      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
      real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
      real(rt), intent(in) :: v(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
      real(rt), intent(in) :: w(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
      real(rt), intent(inout) :: Ipx(Ipx_lo(1):Ipx_hi(1),Ipx_lo(2):Ipx_hi(2),Ipx_lo(3):Ipx_hi(3))
      real(rt), intent(inout) :: Imx(Imx_lo(1):Imx_hi(1),Imx_lo(2):Imx_hi(2),Imx_lo(3):Imx_hi(3))
      real(rt), intent(inout) :: Ipy(Ipy_lo(1):Ipy_hi(1),Ipy_lo(2):Ipy_hi(2),Ipy_lo(3):Ipy_hi(3))
      real(rt), intent(inout) :: Imy(Imy_lo(1):Imy_hi(1),Imy_lo(2):Imy_hi(2),Imy_lo(3):Imy_hi(3))
      real(rt), intent(inout) :: Ipz(Ipz_lo(1):Ipz_hi(1),Ipz_lo(2):Ipz_hi(2),Ipz_lo(3):Ipz_hi(3))
      real(rt), intent(inout) :: Imz(Imz_lo(1):Imz_hi(1),Imz_lo(2):Imz_hi(2),Imz_lo(3):Imz_hi(3))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2),sp_lo(3):sp_hi(3))
      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2),dsvl_lo(3):dsvl_hi(3))
      real(rt), intent(inout) :: sedgex(sedgex_lo(1):sedgex_hi(1),sedgex_lo(2):sedgex_hi(2),sedgex_lo(3):sedgex_hi(3))
      real(rt), intent(inout) :: sedgey(sedgey_lo(1):sedgey_hi(1),sedgey_lo(2):sedgey_hi(2),sedgey_lo(3):sedgey_hi(3))
      real(rt), intent(inout) :: sedgez(sedgez_lo(1):sedgez_hi(1),sedgez_lo(2):sedgez_hi(2),sedgez_lo(3):sedgez_hi(3))
      real(rt), intent(in) :: eps, dx(SDIM), dt
      integer ppm_type

      integer i, j, k

      real(rt) sigma, s6, idtx, idty, idtz

      call ppm_xdir(s,s_lo,s_hi, &
                    sm,sm_lo,sm_hi, &
                    sp,sp_lo,sp_hi, &
                    dsvl,dsvl_lo,dsvl_hi, &
                    sedgex,sedgex_lo,sedgex_hi,&
                    lo,hi,bc,ppm_type)

      idtx = dt / dx(1)

      !
      ! Compute x-component of Ip and Im.
      !
      do k=lo(3)-1,hi(3)+1
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               s6    = 6.0d0*s(i,j,k) - 3.0d0*(sm(i,j,k)+sp(i,j,k))
               sigma = abs(u(i,j,k))*idtx
               if (u(i,j,k) .gt. eps) then
                  Ipx(i,j,k) = sp(i,j,k) - (sigma*half)* &
                      (sp(i,j,k)-sm(i,j,k)-(1.0d0-two3rd*sigma)*s6)
                  Imx(i,j,k) = s(i,j,k)
               else if (u(i,j,k) .lt. -eps) then
                  Ipx(i,j,k) = s(i,j,k)
                  Imx(i,j,k) = sm(i,j,k) + (sigma*half)* &
                      (sp(i,j,k)-sm(i,j,k)+(1.0d0-two3rd*sigma)*s6)
               else
                  Ipx(i,j,k) = s(i,j,k)
                  Imx(i,j,k) = s(i,j,k)
               end if
            end do
         end do
      end do

      call ppm_ydir(s,s_lo,s_hi, &
                    sm,sm_lo,sm_hi, &
                    sp,sp_lo,sp_hi, &
                    dsvl,dsvl_lo,dsvl_hi, &
                    sedgey,sedgey_lo,sedgey_hi,&
                    lo,hi,bc,ppm_type)

      idty = dt / dx(2)

      !
      ! Compute y-component of Ip and Im.
      !
      do k=lo(3)-1,hi(3)+1
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               s6    = 6.0d0*s(i,j,k) - 3.0d0*(sm(i,j,k)+sp(i,j,k))
               sigma = abs(v(i,j,k))*idty
               if (v(i,j,k) .gt. eps) then
                  Ipy(i,j,k) = sp(i,j,k) - (sigma*half)* &
                      (sp(i,j,k)-sm(i,j,k)-(1.0d0-two3rd*sigma)*s6)
                  Imy(i,j,k) = s(i,j,k)
               else if (v(i,j,k) .lt. -eps) then
                  Ipy(i,j,k) = s(i,j,k)
                  Imy(i,j,k) = sm(i,j,k) + (sigma*half)* &
                      (sp(i,j,k)-sm(i,j,k)+(1.0d0-two3rd*sigma)*s6)
               else
                  Ipy(i,j,k) = s(i,j,k)
                  Imy(i,j,k) = s(i,j,k)
               end if
            end do
         end do
      end do
      
      call ppm_zdir(s,s_lo,s_hi, &
                    sm,sm_lo,sm_hi, &
                    sp,sp_lo,sp_hi, &
                    dsvl,dsvl_lo,dsvl_hi, &
                    sedgez,sedgez_lo,sedgez_hi,&
                    lo,hi,bc,ppm_type)

      idtz = dt / dx(3)

      !
      ! Compute z-component of Ip and Im.
      !
      do k=lo(3)-1,hi(3)+1
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               s6    = 6.0d0*s(i,j,k) - 3.0d0*(sm(i,j,k)+sp(i,j,k))
               sigma = abs(w(i,j,k))*idtz
               if (w(i,j,k) .gt. eps) then
                  Ipz(i,j,k) = sp(i,j,k) - (sigma*half)* &
                      (sp(i,j,k)-sm(i,j,k)-(1.0d0-two3rd*sigma)*s6)
                  Imz(i,j,k) = s(i,j,k)
               else if (w(i,j,k) .lt. -eps) then
                  Ipz(i,j,k) = s(i,j,k)
                  Imz(i,j,k) = sm(i,j,k) + (sigma*half)* &
                      (sp(i,j,k)-sm(i,j,k)+(1.0d0-two3rd*sigma)*s6)
               else
                  Ipz(i,j,k) = s(i,j,k)
                  Imz(i,j,k) = s(i,j,k)
               end if
            end do
         end do
      end do

      end subroutine ppm
            
      subroutine ppm_fpu(lo,hi,&
         s,s_lo,s_hi,&
         uedge,uedge_lo,uedge_hi,&
         vedge,vedge_lo,vedge_hi,&
         wedge,wedge_lo,wedge_hi,&
         Ipx,Ipx_lo,Ipx_hi,&
         Imx,Imx_lo,Imx_hi,&
         Ipy,Ipy_lo,Ipy_hi,&
         Imy,Imy_lo,Imy_hi,&
         Ipz,Ipz_lo,Ipz_hi,&
         Imz,Imz_lo,Imz_hi,&
         sm,sm_lo,sm_hi,&
         sp,sp_lo,sp_hi,&
         dsvl,dsvl_lo,dsvl_hi,&
         sedgex,sedgex_lo,sedgex_hi,&
         sedgey,sedgey_lo,sedgey_hi,&
         sedgez,sedgez_lo,sedgez_hi,&
         dx, dt, bc, eps,ppm_type)

      implicit none

      integer, intent(in) :: bc(SDIM,2)
      integer, dimension(3), intent(in) :: &
           s_lo,s_hi,uedge_lo,uedge_hi,vedge_lo,vedge_hi,wedge_lo,wedge_hi,&
           Ipx_lo,Ipx_hi,Imx_lo,Imx_hi,Ipy_lo,Ipy_hi,Imy_lo,Imy_hi,Ipz_lo,Ipz_hi,Imz_lo,Imz_hi,&
           sm_lo,sm_hi,sp_lo,sp_hi,dsvl_lo,dsvl_hi,sedgex_lo,sedgex_hi,sedgey_lo,sedgey_hi,sedgez_lo,sedgez_hi,&
           lo,hi

      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
      real(rt), intent(in) :: uedge(uedge_lo(1):uedge_hi(1),uedge_lo(2):uedge_hi(2),uedge_lo(3):uedge_hi(3))
      real(rt), intent(in) :: vedge(vedge_lo(1):vedge_hi(1),vedge_lo(2):vedge_hi(2),vedge_lo(3):vedge_hi(3))
      real(rt), intent(in) :: wedge(wedge_lo(1):wedge_hi(1),wedge_lo(2):wedge_hi(2),wedge_lo(3):wedge_hi(3))
      real(rt), intent(inout) :: Ipx(Ipx_lo(1):Ipx_hi(1),Ipx_lo(2):Ipx_hi(2),Ipx_lo(3):Ipx_hi(3))
      real(rt), intent(inout) :: Imx(Imx_lo(1):Imx_hi(1),Imx_lo(2):Imx_hi(2),Imx_lo(3):Imx_hi(3))
      real(rt), intent(inout) :: Ipy(Ipy_lo(1):Ipy_hi(1),Ipy_lo(2):Ipy_hi(2),Ipy_lo(3):Ipy_hi(3))
      real(rt), intent(inout) :: Imy(Imy_lo(1):Imy_hi(1),Imy_lo(2):Imy_hi(2),Imy_lo(3):Imy_hi(3))
      real(rt), intent(inout) :: Ipz(Ipz_lo(1):Ipz_hi(1),Ipz_lo(2):Ipz_hi(2),Ipz_lo(3):Ipz_hi(3))
      real(rt), intent(inout) :: Imz(Imz_lo(1):Imz_hi(1),Imz_lo(2):Imz_hi(2),Imz_lo(3):Imz_hi(3))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2),sp_lo(3):sp_hi(3))
      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2),dsvl_lo(3):dsvl_hi(3))
      real(rt), intent(inout) :: sedgex(sedgex_lo(1):sedgex_hi(1),sedgex_lo(2):sedgex_hi(2),sedgex_lo(3):sedgex_hi(3))
      real(rt), intent(inout) :: sedgey(sedgey_lo(1):sedgey_hi(1),sedgey_lo(2):sedgey_hi(2),sedgey_lo(3):sedgey_hi(3))
      real(rt), intent(inout) :: sedgez(sedgez_lo(1):sedgez_hi(1),sedgez_lo(2):sedgez_hi(2),sedgez_lo(3):sedgez_hi(3))
      real(rt), intent(in) :: eps, dx(SDIM), dt
      integer ppm_type

      integer i, j, k

      real(rt) sigmam, sigmap, s6, idtx, idty, idtz

      call ppm_xdir(s,s_lo,s_hi, &
                    sm,sm_lo,sm_hi, &
                    sp,sp_lo,sp_hi, &
                    dsvl,dsvl_lo,dsvl_hi, &
                    sedgex,sedgex_lo,sedgex_hi,&
                    lo,hi,bc,ppm_type)


      idtx = dt / dx(1)

      !
      ! Compute x-component of Ip and Im.
      !
      do k=lo(3)-1,hi(3)+1
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
             s6     = 6.0d0*s(i,j,k) - 3.0d0*(sm(i,j,k)+sp(i,j,k))
             sigmap = abs(uedge(i+1,j,k))*idtx
             sigmam = abs(uedge(i,j,k)  )*idtx
             if (uedge(i+1,j,k) .gt. eps) then
                Ipx(i,j,k) = sp(i,j,k) - (sigmap*half)* &
                    (sp(i,j,k)-sm(i,j,k)-(1.0d0-two3rd*sigmap)*s6)
             else
                Ipx(i,j,k) = s(i,j,k)
             end if
             if (uedge(i,j,k) .lt. -eps) then
                Imx(i,j,k) = sm(i,j,k) + (sigmam*half)* &
                    (sp(i,j,k)-sm(i,j,k)+(1.0d0-two3rd*sigmam)*s6)
             else
                Imx(i,j,k) = s(i,j,k)
             end if
            end do
         end do
      end do

      call ppm_ydir(s,s_lo,s_hi, &
                    sm,sm_lo,sm_hi, &
                    sp,sp_lo,sp_hi, &
                    dsvl,dsvl_lo,dsvl_hi, &
                    sedgey,sedgey_lo,sedgey_hi,&
                    lo,hi,bc,ppm_type)


      idty = dt / dx(2)

      !
      ! Compute y-component of Ip and Im.
      !
      do k=lo(3)-1,hi(3)+1
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
             s6     = 6.0d0*s(i,j,k) - 3.0d0*(sm(i,j,k)+sp(i,j,k))
             sigmap = abs(vedge(i,j+1,k))*idty
             sigmam = abs(vedge(i,j,k)  )*idty
             if (vedge(i,j+1,k) .gt. eps) then
                Ipy(i,j,k) = sp(i,j,k) - (sigmap*half)* &
                    (sp(i,j,k)-sm(i,j,k)-(1.0d0-two3rd*sigmap)*s6)
             else
                Ipy(i,j,k) = s(i,j,k)
             end if
             if (vedge(i,j,k) .lt. -eps) then
                Imy(i,j,k) = sm(i,j,k) + (sigmam*half)* &
                    (sp(i,j,k)-sm(i,j,k)+(1.0d0-two3rd*sigmam)*s6)
             else
                Imy(i,j,k) = s(i,j,k)
             end if
            end do
         end do
      end do

      call ppm_zdir(s,s_lo,s_hi, &
                    sm,sm_lo,sm_hi, &
                    sp,sp_lo,sp_hi, &
                    dsvl,dsvl_lo,dsvl_hi, &
                    sedgez,sedgez_lo,sedgez_hi,&
                    lo,hi,bc,ppm_type)


      idtz = dt / dx(3)

      !
      ! Compute z-component of Ip and Im.
      !
      do k=lo(3)-1,hi(3)+1
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
             s6     = 6.0d0*s(i,j,k) - 3.0d0*(sm(i,j,k)+sp(i,j,k))
             sigmap = abs(wedge(i,j,k+1))*idtz
             sigmam = abs(wedge(i,j,k)  )*idtz
             if (wedge(i,j,k+1) .gt. eps) then
                Ipz(i,j,k) = sp(i,j,k) - (sigmap*half)* &
                    (sp(i,j,k)-sm(i,j,k)-(1.0d0-two3rd*sigmap)*s6)
             else
                Ipz(i,j,k) = s(i,j,k)
             end if
             if (wedge(i,j,k) .lt. -eps) then
                Imz(i,j,k) = sm(i,j,k) + (sigmam*half)* &
                    (sp(i,j,k)-sm(i,j,k)+(1.0d0-two3rd*sigmam)*s6)
             else
                Imz(i,j,k) = s(i,j,k)
             end if
            end do
         end do
      end do

      end subroutine ppm_fpu      
  
      subroutine adv_forcing( &
          aofs,DIMS(aofs), &
          xflux,DIMS(xflux), &
          uedge,DIMS(uedge), &
          areax,DIMS(ax), &
          yflux,DIMS(yflux), &
          vedge,DIMS(vedge), &
          areay,DIMS(ay), &
          zflux,DIMS(zflux), &
          wedge,DIMS(wedge), &
          areaz,DIMS(az), &
          vol,DIMS(vol), &
          lo,hi,iconserv ) bind(C,name="adv_forcing")
!c
!c     This subroutine uses scalar edge states to compute
!c     an advective tendency
!c
      implicit none
      integer i,j,k
      integer iconserv
      real(rt) divux,divuy,divuz
      integer imin,jmin,kmin,imax,jmax,kmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(aofs)
      integer DIMDEC(vol)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
      integer DIMDEC(wedge)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer DIMDEC(az)
      real(rt) aofs(DIMV(aofs))
      real(rt) vol(DIMV(vol))
      real(rt) uedge(DIMV(uedge))
      real(rt) vedge(DIMV(vedge))
      real(rt) wedge(DIMV(wedge))
      real(rt) xflux(DIMV(xflux))
      real(rt) yflux(DIMV(yflux))
      real(rt) zflux(DIMV(zflux))
      real(rt) areax(DIMV(ax))
      real(rt) areay(DIMV(ay))
      real(rt) areaz(DIMV(az))

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
!c
!c     if nonconservative initialize the advective tendency as -U*grad(S)
!c

      if ( iconserv .ne. 1 ) then
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  divux = ( &
                      areax(i+1,j,k)*uedge(i+1,j,k)- &
                      areax(i,  j,k)*uedge(i,  j,k))
                  divuy = ( &
                      areay(i,j+1,k)*vedge(i,j+1,k)- &
                      areay(i,j,  k)*vedge(i,j,  k))
                  divuz = ( &
                      areaz(i,j,k+1)*wedge(i,j,k+1)- &
                      areaz(i,j,k  )*wedge(i,j,k  ))
                  aofs(i,j,k) = &
                     ( - divux*half*(xflux(i+1,j,k)+xflux(i,j,k)) &
                       - divuy*half*(yflux(i,j+1,k)+yflux(i,j,k)) &
                       - divuz*half*(zflux(i,j,k+1)+zflux(i,j,k)) ) /vol(i,j,k)
              
               end do
            end do
         end do
      end if
!c
!c     convert edge states to fluxes
!c

      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax+1
               xflux(i,j,k) = xflux(i,j,k)*uedge(i,j,k)*areax(i,j,k)
            end do
         end do
      end do

      do k = kmin,kmax
         do j = jmin,jmax+1
            do i = imin,imax
               yflux(i,j,k) = yflux(i,j,k)*vedge(i,j,k)*areay(i,j,k)
            end do
         end do
      end do

      do k = kmin,kmax+1
         do j = jmin,jmax
            do i = imin,imax
               zflux(i,j,k) = zflux(i,j,k)*wedge(i,j,k)*areaz(i,j,k)
            end do
         end do
      end do

!c
!c     compute the part of the advective tendency 
!c     that depends on the flux convergence
!c
      if ( iconserv .ne. 1 ) then
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  aofs(i,j,k) = aofs(i,j,k) + ( &
                      xflux(i+1,j,k) - xflux(i,j,k) + &
                      yflux(i,j+1,k) - yflux(i,j,k) + &
                      zflux(i,j,k+1) - zflux(i,j,k))/vol(i,j,k)
               end do
            end do
         end do
      else
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  aofs(i,j,k) = ( &
                      xflux(i+1,j,k) - xflux(i,j,k) + &
                      yflux(i,j+1,k) - yflux(i,j,k) + &
                      zflux(i,j,k+1) - zflux(i,j,k))/vol(i,j,k)
               end do
            end do
         end do
      end if

      end subroutine adv_forcing

      subroutine sync_adv_forcing( &
          sync ,DIMS(sync), &
          xflux,DIMS(xflux), &
          ucor ,DIMS(ucor), &
          areax,DIMS(ax), &
          yflux,DIMS(yflux), &
          vcor ,DIMS(vcor), &
          areay,DIMS(ay), &
          zflux,DIMS(zflux), &
          wcor ,DIMS(wcor), &
          areaz,DIMS(az),  &
          vol ,DIMS(vol), &
          lo,hi) bind(C,name="sync_adv_forcing")
!c
!c     This subroutine computes the sync advective tendency
!c     for a state variable
!c
      implicit none
      integer i,j,k
      integer imin,jmin,kmin,imax,jmax,kmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(sync)
      integer DIMDEC(vol)
      integer DIMDEC(ucor)
      integer DIMDEC(vcor)
      integer DIMDEC(wcor)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer DIMDEC(az)
      real(rt) sync(DIMV(sync))
      real(rt) vol(DIMV(vol))
      real(rt) ucor(DIMV(ucor))
      real(rt) vcor(DIMV(vcor))
      real(rt) wcor(DIMV(wcor))
      real(rt) xflux(DIMV(xflux))
      real(rt) yflux(DIMV(yflux))
      real(rt) zflux(DIMV(zflux))
      real(rt) areax(DIMV(ax))
      real(rt) areay(DIMV(ay))
      real(rt) areaz(DIMV(az))

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
!c
!c     compute corrective fluxes from edge states 
!c     and perform conservative update
!c

      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax+1
               xflux(i,j,k) = xflux(i,j,k)*ucor(i,j,k)*areax(i,j,k)
            end do
         end do
      end do

      do k = kmin,kmax
         do j = jmin,jmax+1
            do i = imin,imax
               yflux(i,j,k) = yflux(i,j,k)*vcor(i,j,k)*areay(i,j,k)
            end do
         end do
      end do

      do k = kmin,kmax+1
         do j = jmin,jmax
            do i = imin,imax
               zflux(i,j,k) = zflux(i,j,k)*wcor(i,j,k)*areaz(i,j,k)
            end do
         end do
      end do

      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               sync(i,j,k) = sync(i,j,k) + ( &
                   xflux(i+1,j,k)-xflux(i,j,k) + &
                   yflux(i,j+1,k)-yflux(i,j,k) + &
                   zflux(i,j,k+1)-zflux(i,j,k) )/vol(i,j,k)
            end do
         end do
      end do

      end subroutine sync_adv_forcing

      subroutine ppm_zdir_colella (s,s_lo,s_hi,&
                                   sm,sm_lo,sm_hi, &
                                   sp,sp_lo,sp_hi, &
                                   sedgez,sedgez_lo,sedgez_hi,lo,hi,klo,khi)

      implicit none

      integer, intent(in) :: lo(SDIM), hi(SDIM), klo, khi
      integer, dimension(3), intent(in) :: s_lo,s_hi,sm_lo,sm_hi,sp_lo,sp_hi,sedgez_lo,sedgez_hi
      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2),sp_lo(3):sp_hi(3))
      real(rt), intent(in) :: sedgez(sedgez_lo(1):sedgez_hi(1),sedgez_lo(2):sedgez_hi(2),sedgez_lo(3):sedgez_hi(3))

      integer i, j, k

      logical extremum, bigp, bigm

      real(rt) D2, D2C, D2L, D2R, D2LIM, alphap, alpham
      real(rt) dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin
      real(rt) sgn, amax, delam, delap, dachkm, dachkp

      real(rt), PARAMETER :: C          = 1.25d0
      real(rt), PARAMETER :: three4ths  = 3.d0/4.d0
      real(rt), PARAMETER :: one20th    = 1.d0/20.0d0
      real(rt), PARAMETER :: one12th    = 1.d0/12.d0
      real(rt), PARAMETER :: seven12ths = 7.d0/12.d0
       !
       ! Use Colella 2008 limiters.
       ! This is a new version of the algorithm
       ! to eliminate sensitivity to roundoff.
       !

      do k=klo,khi
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1

               alphap = sedgez(i,j,k+1)-s(i,j,k)
               alpham = sedgez(i,j,k  )-s(i,j,k)
               bigp = abs(alphap).gt.2.d0*abs(alpham)
               bigm = abs(alpham).gt.2.d0*abs(alphap)
               extremum = .false.

               if (alpham*alphap .ge. 0.d0) then
                  extremum = .true.
               else if (bigp .or. bigm) then
                  !
                  ! Possible extremum. We look at cell centered
                  ! values and face centered values for a change
                  ! in sign in the differences adjacent to
                  ! the cell. We use the pair of differences whose
                  ! minimum magnitude is the largest, and thus least
                  ! susceptible to sensitivity to roundoff.
                  !
                  dafacem = sedgez(i,j,k) - sedgez(i,j,k-1)
                  dafacep = sedgez(i,j,k+2) - sedgez(i,j,k+1)
                  dabarm = s(i,j,k) - s(i,j,k-1)
                  dabarp = s(i,j,k+1) - s(i,j,k)
                  dafacemin = min(abs(dafacem),abs(dafacep))
                  dabarmin= min(abs(dabarm),abs(dabarp))
                  if (dafacemin.ge.dabarmin) then
                     dachkm = dafacem
                     dachkp = dafacep
                  else
                     dachkm = dabarm
                     dachkp = dabarp
                  endif
                  extremum = (dachkm*dachkp .le. 0.d0)
               end if

               if (extremum) then
                  D2  = 6.d0*(alpham + alphap)
                  D2L = s(i,j,k-2)-2.d0*s(i,j,k-1)+s(i,j,k)
                  D2R = s(i,j,k)-2.d0*s(i,j,k+1)+s(i,j,k+2)
                  D2C = s(i,j,k-1)-2.d0*s(i,j,k)+s(i,j,k+1)
                  sgn = sign(1.d0,D2)
                  D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                  alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                  alphap = alphap*D2LIM/max(abs(D2),1.d-10)
               else
                  if (bigp) then
                     sgn = sign(1.d0,alpham)
                     amax = -alphap**2 / (4*(alpham + alphap))
                     delam = s(i,j,k-1) - s(i,j,k)
                     if (sgn*amax .ge. sgn*delam) then
                        if (sgn*(delam - alpham).ge.1.d-10) then
                           alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                        else 
                           alphap = -2.d0*alpham
                        endif
                     endif
                  end if
                  if (bigm) then
                     sgn = sign(1.d0,alphap)
                     amax = -alpham**2 / (4*(alpham + alphap))
                     delap = s(i,j,k+1) - s(i,j,k)
                     if (sgn*amax .ge. sgn*delap) then
                        if (sgn*(delap - alphap).ge.1.d-10) then
                           alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                        else
                           alpham = -2.d0*alphap
                        endif
                     endif
                  end if
               end if
               
               sm(i,j,k) = s(i,j,k) + alpham
               sp(i,j,k) = s(i,j,k) + alphap

            end do
         end do
      end do

      end subroutine ppm_zdir_colella

      subroutine ppm_zdir (s,s_lo,s_hi, &
                           sm,sm_lo,sm_hi, &
                           sp,sp_lo,sp_hi, &
                           dsvl,dsvl_lo,dsvl_hi, &
                           sedgez,sedgez_lo,sedgez_hi,lo,hi,bc,ppm_type)

      implicit none
      
      integer, intent(in) :: lo(SDIM), hi(SDIM), bc(SDIM,2), ppm_type
      integer, dimension(3), intent(in) :: s_lo,s_hi,sm_lo,sm_hi,sp_lo,sp_hi,dsvl_lo,dsvl_hi,sedgez_lo,sedgez_hi
      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2),sp_lo(3):sp_hi(3))
      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2),dsvl_lo(3):dsvl_hi(3))
      real(rt), intent(inout) :: sedgez(sedgez_lo(1):sedgez_hi(1),sedgez_lo(2):sedgez_hi(2),sedgez_lo(3):sedgez_hi(3))

      integer i, j, k

      real(rt) dsl, dsr, dsc, D2, D2L, D2R, D2LIM, sgn

      real(rt), PARAMETER :: C          = 1.25d0
      real(rt), PARAMETER :: three4ths  = 3.d0/4.d0
      real(rt), PARAMETER :: one20th    = 1.d0/20.0d0
      real(rt), PARAMETER :: one12th    = 1.d0/12.d0
      real(rt), PARAMETER :: seven12ths = 7.d0/12.d0

      if (ppm_type .eq. 1) then

         dsvl = 0.d0

         !
         ! Compute van Leer slopes.
         !
         do k=lo(3)-2,hi(3)+2
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-1,hi(1)+1
                  dsc = 0.5d0 * (s(i,j,k+1) - s(i,j,k-1))
                  dsl = 2.d0  * (s(i,j,k  ) - s(i,j,k-1))
                  dsr = 2.d0  * (s(i,j,k+1) - s(i,j,k  ))
                  if (dsl*dsr .gt. 0.d0)  &
                      dsvl(i,j,k) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
               end do
            end do
         end do

         !
         ! Interpolate s to edges.
         !
         do k=lo(3)-1,hi(3)+2
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-1,hi(1)+1
                  sedgez(i,j,k) = 0.5d0*(s(i,j,k)+s(i,j,k-1)) - (sixth)*(dsvl(i,j,k)-dsvl(i,j,k-1))
                  !
                  ! Make sure edge lies between adjacent
                  ! cell-centered values.
                  !
                  sedgez(i,j,k) = max(sedgez(i,j,k),min(s(i,j,k),s(i,j,k-1)))
                  sedgez(i,j,k) = min(sedgez(i,j,k),max(s(i,j,k),s(i,j,k-1)))
               end do
            end do
         end do

         do k=lo(3)-1,hi(3)+1
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-1,hi(1)+1
                  !
                  ! Copy sedge into sp and sm.
                  !
                  sp(i,j,k) = sedgez(i,j,k+1)
                  sm(i,j,k) = sedgez(i,j,k  )
                  !
                  ! Modify using quadratic limiters.
                  !
                  if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. 0.d0) then
                     sp(i,j,k) = s(i,j,k)
                     sm(i,j,k) = s(i,j,k)
                  else if (abs(sp(i,j,k)-s(i,j,k)) .ge. 2.d0*abs(sm(i,j,k)-s(i,j,k))) then
                     sp(i,j,k) = 3.d0*s(i,j,k) - 2.d0*sm(i,j,k)
                  else if (abs(sm(i,j,k)-s(i,j,k)) .ge. 2.d0*abs(sp(i,j,k)-s(i,j,k))) then
                     sm(i,j,k) = 3.d0*s(i,j,k) - 2.d0*sp(i,j,k)
                  end if
               end do
            end do
         end do

         !
         ! Different stencil needed for z-component of
         ! EXT_DIR and HOEXTRAP bc's.
         !
         if (bc(3,1) .eq. EXT_DIR  .or. bc(3,1) .eq. HOEXTRAP) then
            !
            ! The value in the first cc ghost cell
            ! represents the edge value.
            !
            sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)

            k = lo(3)+1

            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-1,hi(1)+1
                  !
                  ! Use a modified stencil to get sedgez
                  ! on the first interior edge.
                  !
                  sedgez(i,j,lo(3)+1) = - fifth     *s(i,j,lo(3)-1) &
                                       + three4ths *s(i,j,lo(3)  ) &
                                       + half      *s(i,j,lo(3)+1) &
                                       - one20th   *s(i,j,lo(3)+2)
                  !
                  ! Make sure sedgez lies in between adjacent
                  ! cell-centered values.
                  !
                  sedgez(i,j,lo(3)+1) = max(sedgez(i,j,lo(3)+1),min(s(i,j,lo(3)+1),s(i,j,lo(3))))
                  sedgez(i,j,lo(3)+1) = min(sedgez(i,j,lo(3)+1),max(s(i,j,lo(3)+1),s(i,j,lo(3))))
                  !
                  ! Copy sedgez into sp and sm.
                  !
                  sp(i,j,lo(3)  ) = sedgez(i,j,lo(3)+1)
                  sm(i,j,lo(3)+1) = sedgez(i,j,lo(3)+1)
                  !
                  ! Reset sp on second interior edge.
                  !
                  sp(i,j,lo(3)+1) = sedgez(i,j,lo(3)+2)
                  !
                  ! Modify using quadratic limiters.
                  !
                  if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. 0.0d0) then
                     sp(i,j,k) = s(i,j,k)
                     sm(i,j,k) = s(i,j,k)
                  else if (abs(sp(i,j,k)-s(i,j,k)) .ge. 2.0d0*abs(sm(i,j,k)-s(i,j,k))) then
                     sp(i,j,k) = 3.0d0*s(i,j,k) - 2.0d0*sm(i,j,k)
                  else if (abs(sm(i,j,k)-s(i,j,k)) .ge. 2.0d0*abs(sp(i,j,k)-s(i,j,k))) then
                     sm(i,j,k) = 3.0d0*s(i,j,k) - 2.0d0*sp(i,j,k)
                  end if
               end do
            end do

         end if

         if (bc(3,2) .eq. EXT_DIR  .or. bc(3,2) .eq. HOEXTRAP) then
            !
            ! The value in the first cc ghost cell
            ! represents the edge value.
            !
            sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)) = s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)

            k = hi(3)-1
            
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-1,hi(1)+1
                  !
                  ! Use a modified stencil to get sedgez on
                  ! the first interior edge.
                  !
                  sedgez(i,j,hi(3)) = - fifth     *s(i,j,hi(3)+1) &
                                     + three4ths *s(i,j,hi(3)  ) &
                                     + half      *s(i,j,hi(3)-1) &
                                     - one20th   *s(i,j,hi(3)-2)
                  !
                  ! Make sure sedgez lies in between adjacent
                  ! cell-centered values.
                  !
                  sedgez(i,j,hi(3)) = max(sedgez(i,j,hi(3)),min(s(i,j,hi(3)-1),s(i,j,hi(3))))
                  sedgez(i,j,hi(3)) = min(sedgez(i,j,hi(3)),max(s(i,j,hi(3)-1),s(i,j,hi(3))))
                  !
                  ! Copy sedgez into sp and sm.
                  !
                  sp(i,j,hi(3)-1) = sedgez(i,j,hi(3))
                  sm(i,j,hi(3)  ) = sedgez(i,j,hi(3))
                  !
                  ! Reset sm on second interior edge.
                  !
                  sm(i,j,hi(3)-1) = sedgez(i,j,hi(3)-1)
                  !
                  ! Modify using quadratic limiters.
                  !
                  if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. 0.0d0) then
                     sp(i,j,k) = s(i,j,k)
                     sm(i,j,k) = s(i,j,k)
                  else if (abs(sp(i,j,k)-s(i,j,k)) .ge. 2.0d0*abs(sm(i,j,k)-s(i,j,k))) then
                     sp(i,j,k) = 3.0d0*s(i,j,k) - 2.0d0*sm(i,j,k)
                  else if (abs(sm(i,j,k)-s(i,j,k)) .ge. 2.0d0*abs(sp(i,j,k)-s(i,j,k))) then
                     sm(i,j,k) = 3.0d0*s(i,j,k) - 2.0d0*sp(i,j,k)
                  end if
               end do
            end do

         end if

      else if (ppm_type .eq. 2) then

         !
         ! Interpolate s to z-edges.
         !
         do k=lo(3)-2,hi(3)+3
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-1,hi(1)+1
                  sedgez(i,j,k) = (seven12ths)*(s(i,j,k-1)+s(i,j,k)) &
                      - (one12th)*(s(i,j,k-2)+s(i,j,k+1))
                  !
                  ! Limit sedgez.
                  !
                  if ((sedgez(i,j,k)-s(i,j,k-1))*(s(i,j,k)-sedgez(i,j,k)) .lt. 0.d0) then
                     D2  = 3.d0*(s(i,j,k-1)-2.d0*sedgez(i,j,k)+s(i,j,k))
                     D2L = s(i,j,k-2)-2.d0*s(i,j,k-1)+s(i,j,k)
                     D2R = s(i,j,k-1)-2.d0*s(i,j,k)+s(i,j,k+1)
                     sgn = sign(1.d0,D2)
                     D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                     sedgez(i,j,k) = 0.5d0*(s(i,j,k-1)+s(i,j,k)) - (sixth)*D2LIM
                  end if
               end do
            end do
         end do

         call ppm_zdir_colella(s,s_lo,s_hi,sm,sm_lo,sm_hi,sp,sp_lo,sp_hi, &
             sedgez,sedgez_lo,sedgez_hi,lo,hi,lo(3)-1,hi(3)+1)
         !
         ! Different stencil needed for z-component of
         ! EXT_DIR and HOEXTRAP bc's.
         !
         if (bc(3,1) .eq. EXT_DIR  .or. bc(3,1) .eq. HOEXTRAP) then
            !
            ! The value in the first cc ghost cell
            ! represents the edge value.
            !
            sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3))     = s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)
            sedgez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)

            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-1,hi(1)+1
                  !
                  ! Use a modified stencil to get sedgez on the
                  ! first interior edge.
                  !
                  sedgez(i,j,lo(3)+1) = - fifth     *s(i,j,lo(3)-1) &
                                       + three4ths *s(i,j,lo(3)  ) &
                                       + 0.5d0     *s(i,j,lo(3)+1) &
                                       - one20th   *s(i,j,lo(3)+2)
                  !
                  ! Make sure sedgez lies in between adjacent
                  ! cell-centered values.
                  !
                  sedgez(i,j,lo(3)+1) = max(sedgez(i,j,lo(3)+1),min(s(i,j,lo(3)+1),s(i,j,lo(3))))
                  sedgez(i,j,lo(3)+1) = min(sedgez(i,j,lo(3)+1),max(s(i,j,lo(3)+1),s(i,j,lo(3))))
                  !
                  ! Copy sedgez into sp.
                  !
                  sp(i,j,lo(3)  ) = sedgez(i,j,lo(3)+1)
               end do
            end do
            !
            ! Apply Colella 2008 limiters to compute sm and sp
            ! in the 2nd and 3rd inner cells.
            !
            call ppm_zdir_colella(s,s_lo,s_hi,sm,sm_lo,sm_hi,sp,sp_lo,sp_hi, &
                sedgez,sedgez_lo,sedgez_hi,lo,hi,lo(3)+1,lo(3)+2)
         end if

         if (bc(3,2) .eq. EXT_DIR  .or. bc(3,2) .eq. HOEXTRAP) then
            !
            ! The value in the first cc ghost cell
            ! represents the edge value.
            !
            sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3))       = s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)
            sedgez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1) = s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)

            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-1,hi(1)+1
                  !
                  ! Use a modified stencil to get sedgez on the
                  ! first interior edge.
                  !
                  sedgez(i,j,hi(3)) = - fifth     *s(i,j,hi(3)+1) &
                                     + three4ths *s(i,j,hi(3)  ) &
                                     + 0.5d0     *s(i,j,hi(3)-1) &
                                     - one20th   *s(i,j,hi(3)-2)
                  !
                  ! Make sure sedgez lies in between adjacent
                  ! cell-centered values.
                  !
                  sedgez(i,j,hi(3)) = max(sedgez(i,j,hi(3)),min(s(i,j,hi(3)-1),s(i,j,hi(3))))
                  sedgez(i,j,hi(3)) = min(sedgez(i,j,hi(3)),max(s(i,j,hi(3)-1),s(i,j,hi(3))))
                  !
                  ! Copy sedgez into sp.
                  !
                  sm(i,j,hi(3)  ) = sedgez(i,j,hi(3))
               end do
            end do

            !
            ! Apply Colella 2008 limiters to compute sm and sp
            ! in the 2nd and 3rd inner cells.
            !
            call ppm_zdir_colella(s,s_lo,s_hi,sm,sm_lo,sm_hi,sp,sp_lo,sp_hi, &
                sedgez,sedgez_lo,sedgez_hi,lo,hi,hi(3)-2,hi(3)-1)
         end if

      end if

      end subroutine ppm_zdir


      subroutine ppm_ydir_colella (s,s_lo,s_hi,&
                                   sm,sm_lo,sm_hi, &
                                   sp,sp_lo,sp_hi, &
                                   sedgey,sedgey_lo,sedgey_hi,lo,hi,jlo,jhi)

      implicit none

      integer, intent(in) :: lo(SDIM), hi(SDIM), jlo, jhi
      integer, dimension(3), intent(in) :: s_lo,s_hi,sm_lo,sm_hi,sp_lo,sp_hi,sedgey_lo,sedgey_hi
      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2),sp_lo(3):sp_hi(3))
      real(rt), intent(in) :: sedgey(sedgey_lo(1):sedgey_hi(1),sedgey_lo(2):sedgey_hi(2),sedgey_lo(3):sedgey_hi(3))

      integer i, j, k

      logical extremum, bigp, bigm

      real(rt) D2, D2C, D2L, D2R, D2LIM, alphap, alpham
      real(rt) dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin
      real(rt) sgn, amax, delam, delap, dachkm, dachkp

      real(rt), PARAMETER :: C          = 1.25d0
      real(rt), PARAMETER :: three4ths  = 3.d0/4.d0
      real(rt), PARAMETER :: one20th    = 1.d0/20.0d0
      real(rt), PARAMETER :: one12th    = 1.d0/12.d0
      real(rt), PARAMETER :: seven12ths = 7.d0/12.d0
       !
       ! Use Colella 2008 limiters.
       ! This is a new version of the algorithm
       ! to eliminate sensitivity to roundoff.
       !

      do k=lo(3)-1,hi(3)+1
         do j=jlo,jhi
            do i=lo(1)-1,hi(1)+1

               alphap = sedgey(i,j+1,k)-s(i,j,k)
               alpham = sedgey(i,j  ,k)-s(i,j,k)
               bigp = abs(alphap).gt.2.d0*abs(alpham)
               bigm = abs(alpham).gt.2.d0*abs(alphap)
               extremum = .false.

               if (alpham*alphap .ge. 0.d0) then
                  extremum = .true.
               else if (bigp .or. bigm) then
                  !
                  ! Possible extremum. We look at cell centered
                  ! values and face centered values for a change
                  ! in sign in the differences adjacent to
                  ! the cell. We use the pair of differences whose
                  ! minimum magnitude is the largest, and thus least
                  ! susceptible to sensitivity to roundoff.
                  !
                  dafacem = sedgey(i,j,k) - sedgey(i,j-1,k)
                  dafacep = sedgey(i,j+2,k) - sedgey(i,j+1,k)
                  dabarm = s(i,j,k) - s(i,j-1,k)
                  dabarp = s(i,j+1,k) - s(i,j,k)
                  dafacemin = min(abs(dafacem),abs(dafacep))
                  dabarmin= min(abs(dabarm),abs(dabarp))
                  if (dafacemin.ge.dabarmin) then
                     dachkm = dafacem
                     dachkp = dafacep
                  else
                     dachkm = dabarm
                     dachkp = dabarp
                  endif
                  extremum = (dachkm*dachkp .le. 0.d0)
               end if

               if (extremum) then
                  D2  = 6.d0*(alpham + alphap)
                  D2L = s(i,j-2,k)-2.d0*s(i,j-1,k)+s(i,j,k)
                  D2R = s(i,j,k)-2.d0*s(i,j+1,k)+s(i,j+2,k)
                  D2C = s(i,j-1,k)-2.d0*s(i,j,k)+s(i,j+1,k)
                  sgn = sign(1.d0,D2)
                  D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                  alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                  alphap = alphap*D2LIM/max(abs(D2),1.d-10)
               else
                  if (bigp) then
                     sgn = sign(1.d0,alpham)
                     amax = -alphap**2 / (4*(alpham + alphap))
                     delam = s(i,j-1,k) - s(i,j,k)
                     if (sgn*amax .ge. sgn*delam) then
                        if (sgn*(delam - alpham).ge.1.d-10) then
                           alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                        else 
                           alphap = -2.d0*alpham
                        endif
                     endif
                  end if
                  if (bigm) then
                     sgn = sign(1.d0,alphap)
                     amax = -alpham**2 / (4*(alpham + alphap))
                     delap = s(i,j+1,k) - s(i,j,k)
                     if (sgn*amax .ge. sgn*delap) then
                        if (sgn*(delap - alphap).ge.1.d-10) then
                           alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                        else
                           alpham = -2.d0*alphap
                        endif
                     endif
                  end if
               end if
               
               sm(i,j,k) = s(i,j,k) + alpham
               sp(i,j,k) = s(i,j,k) + alphap

            end do
         end do
      end do

      end subroutine ppm_ydir_colella

      subroutine ppm_ydir (s,s_lo,s_hi, &
                           sm,sm_lo,sm_hi, &
                           sp,sp_lo,sp_hi, &
                           dsvl,dsvl_lo,dsvl_hi, &
                           sedgey,sedgey_lo,sedgey_hi,lo,hi,bc,ppm_type)

      implicit none

      integer, intent(in) :: lo(SDIM), hi(SDIM), bc(SDIM,2), ppm_type
      integer, dimension(3), intent(in) :: s_lo,s_hi,sm_lo,sm_hi,sp_lo,sp_hi,dsvl_lo,dsvl_hi,sedgey_lo,sedgey_hi
      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2),sp_lo(3):sp_hi(3))
      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2),dsvl_lo(3):dsvl_hi(3))
      real(rt), intent(inout) :: sedgey(sedgey_lo(1):sedgey_hi(1),sedgey_lo(2):sedgey_hi(2),sedgey_lo(3):sedgey_hi(3))

      integer i, j, k

      real(rt) dsl, dsr, dsc, D2, D2L, D2R, D2LIM, sgn

      real(rt), PARAMETER :: C          = 1.25d0
      real(rt), PARAMETER :: three4ths  = 3.d0/4.d0
      real(rt), PARAMETER :: one20th    = 1.d0/20.0d0
      real(rt), PARAMETER :: one12th    = 1.d0/12.d0
      real(rt), PARAMETER :: seven12ths = 7.d0/12.d0

      if (ppm_type .eq. 1) then

         dsvl = 0.d0

         !
         ! Compute van Leer slopes.
         !
         do k=lo(3)-1,hi(3)+1
            do j=lo(2)-2,hi(2)+2
               do i=lo(1)-1,hi(1)+1
                  dsc = 0.5d0 * (s(i,j+1,k) - s(i,j-1,k))
                  dsl = 2.d0  * (s(i,j  ,k) - s(i,j-1,k))
                  dsr = 2.d0  * (s(i,j+1,k) - s(i,j  ,k))
                  if (dsl*dsr .gt. 0.d0) &
                      dsvl(i,j,k) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
               end do
            end do
         end do

         !
         ! Interpolate s to edges.
         !
         do k=lo(3)-1,hi(3)+1
            do j=lo(2)-1,hi(2)+2
               do i=lo(1)-1,hi(1)+1
                  sedgey(i,j,k) = 0.5d0*(s(i,j,k)+s(i,j-1,k)) - (sixth)*(dsvl(i,j,k)-dsvl(i,j-1,k))
                  !
                  ! Make sure edge lies between adjacent
                  ! cell-centered values.
                  !
                  sedgey(i,j,k) = max(sedgey(i,j,k),min(s(i,j,k),s(i,j-1,k)))
                  sedgey(i,j,k) = min(sedgey(i,j,k),max(s(i,j,k),s(i,j-1,k)))
               end do
            end do
         end do

         do k=lo(3)-1,hi(3)+1
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-1,hi(1)+1
                  !
                  ! Copy sedge into sp and sm.
                  !
                  sp(i,j,k) = sedgey(i,j+1,k)
                  sm(i,j,k) = sedgey(i,j  ,k)
                  !
                  ! Modify using quadratic limiters.
                  !
                  if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. 0.0d0) then
                     sp(i,j,k) = s(i,j,k)
                     sm(i,j,k) = s(i,j,k)
                  else if (abs(sp(i,j,k)-s(i,j,k)) .ge. 2.d0*abs(sm(i,j,k)-s(i,j,k))) then
                     sp(i,j,k) = 3.d0*s(i,j,k) - 2.d0*sm(i,j,k)
                  else if (abs(sm(i,j,k)-s(i,j,k)) .ge. 2.d0*abs(sp(i,j,k)-s(i,j,k))) then
                     sm(i,j,k) = 3.d0*s(i,j,k) - 2.d0*sp(i,j,k)
                  end if
               end do
            end do
         end do

         !
         ! Different stencil needed for y-component of
         ! EXT_DIR and HOEXTRAP bc's.
         !
         if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
            !
            ! The value in the first cc ghost cell
            ! represents the edge value.
            !
            sm(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)

            j = lo(2)+1

            do k=lo(3)-1,hi(3)+1
               do i=lo(1)-1,hi(1)+1
                  !
                  ! Use a modified stencil to get sedgey
                  ! on the first interior edge.
                  !
                  sedgey(i,lo(2)+1,k) = - fifth     *s(i,lo(2)-1,k) &
                                       + three4ths *s(i,lo(2)  ,k) &
                                       + half      *s(i,lo(2)+1,k) &
                                       - one20th   *s(i,lo(2)+2,k)
                  !
                  ! Make sure sedgey lies in between adjacent
                  ! cell-centered values.
                  !
                  sedgey(i,lo(2)+1,k) = max(sedgey(i,lo(2)+1,k),min(s(i,lo(2)+1,k),s(i,lo(2),k)))
                  sedgey(i,lo(2)+1,k) = min(sedgey(i,lo(2)+1,k),max(s(i,lo(2)+1,k),s(i,lo(2),k)))
                  !
                  ! Copy sedgey into sp and sm.
                  !
                  sp(i,lo(2)  ,k) = sedgey(i,lo(2)+1,k)
                  sm(i,lo(2)+1,k) = sedgey(i,lo(2)+1,k)
                  !
                  ! Reset sp on second interior edge.
                  !
                  sp(i,lo(2)+1,k) = sedgey(i,lo(2)+2,k)
                  !
                  ! Modify using quadratic limiters.
                  !
                  if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. 0.0d0) then
                     sp(i,j,k) = s(i,j,k)
                     sm(i,j,k) = s(i,j,k)
                  else if (abs(sp(i,j,k)-s(i,j,k)) .ge. 2.0d0*abs(sm(i,j,k)-s(i,j,k))) then
                     sp(i,j,k) = 3.0d0*s(i,j,k) - 2.0d0*sm(i,j,k)
                  else if (abs(sm(i,j,k)-s(i,j,k)) .ge. 2.0d0*abs(sp(i,j,k)-s(i,j,k))) then
                     sm(i,j,k) = 3.0d0*s(i,j,k) - 2.0d0*sp(i,j,k)
                  end if
               end do
            end do

         end if

         if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
            !
            ! The value in the first cc ghost cell
            ! represents the edge value.
            !
            sp(lo(1)-1:hi(1)+1,hi(2),lo(3)-1:hi(3)+1) = s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)

            j = hi(2)-1

            do k=lo(3)-1,hi(3)+1
               do i=lo(1)-1,hi(1)+1
                  !
                  ! Use a modified stencil to get sedgey
                  ! on the first interior edge.
                  !
                  sedgey(i,hi(2),k) = - fifth     *s(i,hi(2)+1,k) &
                      + three4ths *s(i,hi(2)  ,k) &
                                     + half      *s(i,hi(2)-1,k) &
                                     - one20th   *s(i,hi(2)-2,k)
                  !
                  ! Make sure sedgey lies in between adjacent
                  ! cell-centered values.
                  !
                  sedgey(i,hi(2),k) = max(sedgey(i,hi(2),k),min(s(i,hi(2)-1,k),s(i,hi(2),k)))
                  sedgey(i,hi(2),k) = min(sedgey(i,hi(2),k),max(s(i,hi(2)-1,k),s(i,hi(2),k)))
                  !
                  ! Copy sedgey into sp and sm.
                  !
                  sp(i,hi(2)-1,k) = sedgey(i,hi(2),k)
                  sm(i,hi(2)  ,k) = sedgey(i,hi(2),k)
                  !
                  ! Reset sm on second interior edge.
                  !
                  sm(i,hi(2)-1,k) = sedgey(i,hi(2)-1,k)
                  !
                  ! Modify using quadratic limiters.
                  !
                  if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. 0.0d0) then
                     sp(i,j,k) = s(i,j,k)
                     sm(i,j,k) = s(i,j,k)
                  else if (abs(sp(i,j,k)-s(i,j,k)) .ge. 2.0d0*abs(sm(i,j,k)-s(i,j,k))) then
                     sp(i,j,k) = 3.0d0*s(i,j,k) - 2.0d0*sm(i,j,k)
                  else if (abs(sm(i,j,k)-s(i,j,k)) .ge. 2.0d0*abs(sp(i,j,k)-s(i,j,k))) then
                     sm(i,j,k) = 3.0d0*s(i,j,k) - 2.0d0*sp(i,j,k)
                  end if
               end do
            end do

         end if

      else if (ppm_type .eq. 2) then

         !
         ! Interpolate s to y-edges.
         !
         do k=lo(3)-1,hi(3)+1
            do j=lo(2)-2,hi(2)+3
               do i=lo(1)-1,hi(1)+1
                  sedgey(i,j,k) = (seven12ths)*(s(i,j-1,k)+s(i,j,k)) &
                      - (one12th)*(s(i,j-2,k)+s(i,j+1,k))
                  !
                  ! Limit sedgey.
                  !
                  if ((sedgey(i,j,k)-s(i,j-1,k))*(s(i,j,k)-sedgey(i,j,k)) .lt. 0.d0) then
                     D2  = 3.d0*(s(i,j-1,k)-2.d0*sedgey(i,j,k)+s(i,j,k))
                     D2L = s(i,j-2,k)-2.d0*s(i,j-1,k)+s(i,j,k)
                     D2R = s(i,j-1,k)-2.d0*s(i,j,k)+s(i,j+1,k)
                     sgn = sign(1.d0,D2)
                     D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                     sedgey(i,j,k) = 0.5d0*(s(i,j-1,k)+s(i,j,k)) - (sixth)*D2LIM
                  end if
               end do
            end do
         end do

         call ppm_ydir_colella(s,s_lo,s_hi,sm,sm_lo,sm_hi,sp,sp_lo,sp_hi, &
             sedgey,sedgey_lo,sedgey_hi,lo,hi,lo(2)-1,hi(2)+1)
         !
         ! Different stencil needed for y-component of
         ! EXT_DIR and HOEXTRAP bc's
         !
         if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
            !
            ! The value in the first cc ghost cell
            ! represents the edge value.
            !
            sm(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1)     = s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)
            sedgey(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)

            do k=lo(3)-1,hi(3)+1
               do i=lo(1)-1,hi(1)+1
                  !
                  ! Use a modified stencil to get sedgey
                  ! on the first interior edge.
                  !
                  sedgey(i,lo(2)+1,k) = - fifth     *s(i,lo(2)-1,k) &
                                       + three4ths *s(i,lo(2)  ,k) &
                                       + 0.5d0     *s(i,lo(2)+1,k) &
                                       - one20th   *s(i,lo(2)+2,k)
                  !
                  ! Make sure sedgey lies in between adjacent
                  ! cell-centered values.
                  !
                  sedgey(i,lo(2)+1,k) = max(sedgey(i,lo(2)+1,k),min(s(i,lo(2)+1,k),s(i,lo(2),k)))
                  sedgey(i,lo(2)+1,k) = min(sedgey(i,lo(2)+1,k),max(s(i,lo(2)+1,k),s(i,lo(2),k)))
                  !
                  ! Copy sedgey into sp.
                  !
                  sp(i,lo(2)  ,k) = sedgey(i,lo(2)+1,k)
               end do
            end do

            !
            ! Apply Colella 2008 limiters to compute sm and sp
            ! in the 2nd and 3rd inner cells.
            !
            call ppm_ydir_colella(s,s_lo,s_hi,sm,sm_lo,sm_hi,sp,sp_lo,sp_hi, &
                sedgey,sedgey_lo,sedgey_hi,lo,hi,lo(2)+1,lo(2)+2)
         end if

         if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
            !
            ! The value in the first cc ghost cell
            ! represents the edge value.
            !
            sp(lo(1)-1:hi(1)+1,hi(2),lo(3)-1:hi(3)+1) = s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)
            sedgey(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1) = s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)

            do k=lo(3)-1,hi(3)+1
               do i=lo(1)-1,hi(1)+1
                  !
                  ! Use a modified stencil to get sedgey
                  ! on the first interior edge.
                  !
                  sedgey(i,hi(2),k) = - fifth     *s(i,hi(2)+1,k) &
                                     + three4ths *s(i,hi(2)  ,k) &
                                     + 0.5d0     *s(i,hi(2)-1,k) &
                                     - one20th   *s(i,hi(2)-2,k)
                  !
                  ! Make sure sedgey lies in between adjacent
                  ! cell-centered values.
                  !
                  sedgey(i,hi(2),k) = max(sedgey(i,hi(2),k),min(s(i,hi(2)-1,k),s(i,hi(2),k)))
                  sedgey(i,hi(2),k) = min(sedgey(i,hi(2),k),max(s(i,hi(2)-1,k),s(i,hi(2),k)))
                  !
                  ! Copy sedgey into sp.
                  !
                  sm(i,hi(2)  ,k) = sedgey(i,hi(2),k)
               end do
            end do

            !
            ! Apply Colella 2008 limiters to compute sm and sp
            ! in the 2nd and 3rd inner cells.
            !
            call ppm_ydir_colella(s,s_lo,s_hi,sm,sm_lo,sm_hi,sp,sp_lo,sp_hi, &
                sedgey,sedgey_lo,sedgey_hi,lo,hi,hi(2)-2,hi(2)-1)
         end if

      end if

      end subroutine ppm_ydir


      subroutine ppm_xdir_colella (s,s_lo,s_hi,&
                                   sm,sm_lo,sm_hi, &
                                   sp,sp_lo,sp_hi, &
                                   sedgex,sedgex_lo,sedgex_hi,lo,hi,ilo,ihi)

      implicit none

      integer, intent(in) :: lo(SDIM), hi(SDIM), ilo, ihi
      integer, dimension(3), intent(in) :: s_lo,s_hi,sm_lo,sm_hi,sp_lo,sp_hi,sedgex_lo,sedgex_hi
      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2),sp_lo(3):sp_hi(3))
      real(rt), intent(in) :: sedgex(sedgex_lo(1):sedgex_hi(1),sedgex_lo(2):sedgex_hi(2),sedgex_lo(3):sedgex_hi(3))

      integer i, j, k

      logical extremum, bigp, bigm

      real(rt) D2, D2C, D2L, D2R, D2LIM, alphap, alpham
      real(rt) dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin
      real(rt) sgn, amax, delam, delap, dachkm, dachkp

      real(rt), PARAMETER :: C          = 1.25d0
      real(rt), PARAMETER :: three4ths  = 3.d0/4.d0
      real(rt), PARAMETER :: one20th    = 1.d0/20.0d0
      real(rt), PARAMETER :: one12th    = 1.d0/12.d0
      real(rt), PARAMETER :: seven12ths = 7.d0/12.d0
       !
       ! Use Colella 2008 limiters.
       ! This is a new version of the algorithm
       ! to eliminate sensitivity to roundoff.
       !

      do k=lo(3)-1,hi(3)+1
         do j=lo(2)-1,hi(2)+1
            do i=ilo,ihi

               alphap = sedgex(i+1,j,k)-s(i,j,k)
               alpham = sedgex(i  ,j,k)-s(i,j,k)
               bigp = abs(alphap).gt.2.d0*abs(alpham)
               bigm = abs(alpham).gt.2.d0*abs(alphap)
               extremum = .false.

               if (alpham*alphap .ge. 0.d0) then
                  extremum = .true.
               else if (bigp .or. bigm) then
                  !
                  ! Possible extremum. We look at cell centered
                  ! values and face centered values for a change
                  ! in sign in the differences adjacent to
                  ! the cell. We use the pair of differences whose
                  ! minimum magnitude is the largest, and thus least
                  ! susceptible to sensitivity to roundoff.
                  !
                  dafacem = sedgex(i,j,k) - sedgex(i-1,j,k)
                  dafacep = sedgex(i+2,j,k) - sedgex(i+1,j,k)
                  dabarm = s(i,j,k) - s(i-1,j,k)
                  dabarp = s(i+1,j,k) - s(i,j,k)
                  dafacemin = min(abs(dafacem),abs(dafacep))
                  dabarmin= min(abs(dabarm),abs(dabarp))
                  if (dafacemin.ge.dabarmin) then
                     dachkm = dafacem
                     dachkp = dafacep
                  else
                     dachkm = dabarm
                     dachkp = dabarp
                  endif
                  extremum = (dachkm*dachkp .le. 0.d0)
               end if

               if (extremum) then
                  D2  = 6.d0*(alpham + alphap)
                  D2L = s(i-2,j,k)-2.d0*s(i-1,j,k)+s(i,j,k)
                  D2R = s(i,j,k)-2.d0*s(i+1,j,k)+s(i+2,j,k)
                  D2C = s(i-1,j,k)-2.d0*s(i,j,k)+s(i+1,j,k)
                  sgn = sign(1.d0,D2)
                  D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                  alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                  alphap = alphap*D2LIM/max(abs(D2),1.d-10)
               else
                  if (bigp) then
                     sgn = sign(1.d0,alpham)
                     amax = -alphap**2 / (4*(alpham + alphap))
                     delam = s(i-1,j,k) - s(i,j,k)
                     if (sgn*amax .ge. sgn*delam) then
                        if (sgn*(delam - alpham).ge.1.d-10) then
                           alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                        else 
                           alphap = -2.d0*alpham
                        endif
                     endif
                  end if
                  if (bigm) then
                     sgn = sign(1.d0,alphap)
                     amax = -alpham**2 / (4*(alpham + alphap))
                     delap = s(i+1,j,k) - s(i,j,k)
                     if (sgn*amax .ge. sgn*delap) then
                        if (sgn*(delap - alphap).ge.1.d-10) then
                           alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                        else
                           alpham = -2.d0*alphap
                        endif
                     endif
                  end if
               end if
               
               sm(i,j,k) = s(i,j,k) + alpham
               sp(i,j,k) = s(i,j,k) + alphap

            end do
         end do
      end do

      end subroutine ppm_xdir_colella

      subroutine ppm_xdir (s,s_lo,s_hi, &
                           sm,sm_lo,sm_hi, &
                           sp,sp_lo,sp_hi, &
                           dsvl,dsvl_lo,dsvl_hi, &
                           sedgex,sedgex_lo,sedgex_hi,lo,hi,bc,ppm_type)

      implicit none
      
      integer, intent(in) :: lo(SDIM), hi(SDIM), bc(SDIM,2), ppm_type
      integer, dimension(3), intent(in) :: s_lo,s_hi,sm_lo,sm_hi,sp_lo,sp_hi,dsvl_lo,dsvl_hi,sedgex_lo,sedgex_hi
      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2),sp_lo(3):sp_hi(3))
      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2),dsvl_lo(3):dsvl_hi(3))
      real(rt), intent(inout) :: sedgex(sedgex_lo(1):sedgex_hi(1),sedgex_lo(2):sedgex_hi(2),sedgex_lo(3):sedgex_hi(3))

      integer i, j, k

      real(rt) dsl, dsr, dsc, D2, D2L, D2R, D2LIM, sgn

      real(rt), PARAMETER :: C          = 1.25d0
      real(rt), PARAMETER :: three4ths  = 3.d0/4.d0
      real(rt), PARAMETER :: one20th    = 1.d0/20.0d0
      real(rt), PARAMETER :: one12th    = 1.d0/12.d0
      real(rt), PARAMETER :: seven12ths = 7.d0/12.d0

      if (ppm_type .eq. 1) then

         dsvl = 0.d0

         !
         ! Compute van Leer slopes.
         !
         do k=lo(3)-1,hi(3)+1
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-2,hi(1)+2
                  dsc = 0.5d0 * (s(i+1,j,k) - s(i-1,j,k))
                  dsl = 2.d0  * (s(i  ,j,k) - s(i-1,j,k))
                  dsr = 2.d0  * (s(i+1,j,k) - s(i  ,j,k))
                  if (dsl*dsr .gt. 0.d0) &
                      dsvl(i,j,k) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
               end do
            end do
         end do

         !
         ! Interpolate s to edges.
         !
         do k=lo(3)-1,hi(3)+1
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-1,hi(1)+2
                  sedgex(i,j,k) = 0.5d0*(s(i,j,k)+s(i-1,j,k)) - (sixth)*(dsvl(i,j,k)-dsvl(i-1,j,k))
                  !
                  ! Make sure edge lies between adjacent
                  ! cell-centered values.
                  !
                  sedgex(i,j,k) = max(sedgex(i,j,k),min(s(i,j,k),s(i-1,j,k)))
                  sedgex(i,j,k) = min(sedgex(i,j,k),max(s(i,j,k),s(i-1,j,k)))
               end do
            end do
         end do

         do k=lo(3)-1,hi(3)+1
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-1,hi(1)+1
                  !
                  ! Copy sedgex into sp and sm.
                  !
                  sp(i,j,k) = sedgex(i+1,j,k)
                  sm(i,j,k) = sedgex(i  ,j,k)
                  !
                  ! Modify using quadratic limiters.
                  !
                  if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. 0.d0) then
                     sp(i,j,k) = s(i,j,k)
                     sm(i,j,k) = s(i,j,k)
                  else if (abs(sp(i,j,k)-s(i,j,k)) .ge. 2.d0*abs(sm(i,j,k)-s(i,j,k))) then
                     sp(i,j,k) = 3.d0*s(i,j,k) - 2.d0*sm(i,j,k)
                  else if (abs(sm(i,j,k)-s(i,j,k)) .ge. 2.d0*abs(sp(i,j,k)-s(i,j,k))) then
                     sm(i,j,k) = 3.d0*s(i,j,k) - 2.d0*sp(i,j,k)
                  end if
               end do
            end do
         end do

         !
         ! Different stencil needed for x-component of
         ! EXT_DIR and HOEXTRAP bc's.
         !
         if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
            !
            ! The value in the first cc ghost cell represents
            ! the edge value.
            !
            sm(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

            i = lo(1)+1

            do k=lo(3)-1,hi(3)+1
               do j=lo(2)-1,hi(2)+1
                  !
                  ! Use a modified stencil to get sedgex
                  ! on the first interior edge.
                  !
                  sedgex(lo(1)+1,j,k) = - fifth     *s(lo(1)-1,j,k) &
                                       + three4ths *s(lo(1)  ,j,k) &
                                       + half      *s(lo(1)+1,j,k) &
                                       - one20th   *s(lo(1)+2,j,k)
                  !
                  ! Make sure sedgex lies in between adjacent
                  ! cell-centered values.
                  !
                  sedgex(lo(1)+1,j,k) = max(sedgex(lo(1)+1,j,k),min(s(lo(1)+1,j,k),s(lo(1),j,k)))
                  sedgex(lo(1)+1,j,k) = min(sedgex(lo(1)+1,j,k),max(s(lo(1)+1,j,k),s(lo(1),j,k)))
                  !
                  ! Copy sedgex into sp and sm.
                  !
                  sp(lo(1)  ,j,k) = sedgex(lo(1)+1,j,k)
                  sm(lo(1)+1,j,k) = sedgex(lo(1)+1,j,k)
                  !
                  ! Reset sp on second interior edge.
                  !
                  sp(lo(1)+1,j,k) = sedgex(lo(1)+2,j,k)
                  !
                  ! Modify using quadratic limiters.
                  !
                  if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. 0.0d0) then
                     sp(i,j,k) = s(i,j,k)
                     sm(i,j,k) = s(i,j,k)
                  else if (abs(sp(i,j,k)-s(i,j,k)) .ge. 2.0d0*abs(sm(i,j,k)-s(i,j,k))) then
                     sp(i,j,k) = 3.0d0*s(i,j,k) - 2.0d0*sm(i,j,k)
                  else if (abs(sm(i,j,k)-s(i,j,k)) .ge. 2.0d0*abs(sp(i,j,k)-s(i,j,k))) then
                     sm(i,j,k) = 3.0d0*s(i,j,k) - 2.0d0*sp(i,j,k)
                  end if
               end do
            end do

         end if

         if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then 
!c     the value in the first cc ghost cell represents the edge value
            sp(hi(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

            i = hi(1)-1

            do k=lo(3)-1,hi(3)+1
               do j=lo(2)-1,hi(2)+1
                  !
                  ! Use a modified stencil to get sedgex
                  ! on the first interior edge.
                  !
                  sedgex(hi(1),j,k) = - fifth     *s(hi(1)+1,j,k) &
                                     + three4ths *s(hi(1)  ,j,k) &
                                     + half      *s(hi(1)-1,j,k) &
                                     - one20th   *s(hi(1)-2,j,k)
                  !
                  ! Make sure sedgex lies in between adjacent
                  ! cell-centered values.
                  !
                  sedgex(hi(1),j,k) = max(sedgex(hi(1),j,k),min(s(hi(1)-1,j,k),s(hi(1),j,k)))
                  sedgex(hi(1),j,k) = min(sedgex(hi(1),j,k),max(s(hi(1)-1,j,k),s(hi(1),j,k)))
                  !
                  ! Copy sedgex into sp and sm.
                  !
                  sp(hi(1)-1,j,k) = sedgex(hi(1),j,k)
                  sm(hi(1)  ,j,k) = sedgex(hi(1),j,k)
                  !
                  ! Reset sm on second interior edge.
                  !
                  sm(hi(1)-1,j,k) = sedgex(hi(1)-1,j,k)
                  !
                  ! Modify using quadratic limiters.
                  !
                  if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. 0.0d0) then
                     sp(i,j,k) = s(i,j,k)
                     sm(i,j,k) = s(i,j,k)
                  else if (abs(sp(i,j,k)-s(i,j,k)) .ge. 2.0d0*abs(sm(i,j,k)-s(i,j,k))) then
                     sp(i,j,k) = 3.0d0*s(i,j,k) - 2.0d0*sm(i,j,k)
                  else if (abs(sm(i,j,k)-s(i,j,k)) .ge. 2.0d0*abs(sp(i,j,k)-s(i,j,k))) then
                     sm(i,j,k) = 3.0d0*s(i,j,k) - 2.0d0*sp(i,j,k)
                  end if
               end do
            end do

         end if

      else if (ppm_type .eq. 2) then

         !
         ! Interpolate s to x-edges.
         !
         do k=lo(3)-1,hi(3)+1
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)-2,hi(1)+3
                  sedgex(i,j,k) = (seven12ths)*(s(i-1,j,k)+s(i,j,k)) &
                      - (one12th)*(s(i-2,j,k)+s(i+1,j,k))
                  !
                  ! Limit sedgex.
                  !
                  if ((sedgex(i,j,k)-s(i-1,j,k))*(s(i,j,k)-sedgex(i,j,k)) .lt. 0.d0) then
                     D2  = 3.d0*(s(i-1,j,k)-2.d0*sedgex(i,j,k)+s(i,j,k))
                     D2L = s(i-2,j,k)-2.d0*s(i-1,j,k)+s(i,j,k)
                     D2R = s(i-1,j,k)-2.d0*s(i,j,k)+s(i+1,j,k)
                     sgn = sign(1.d0,D2)
                     D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                     sedgex(i,j,k) = 0.5d0*(s(i-1,j,k)+s(i,j,k)) - (sixth)*D2LIM
                  end if
               end do
            end do
         end do

      call ppm_xdir_colella(s,s_lo,s_hi,sm,sm_lo,sm_hi,sp,sp_lo,sp_hi, &
          sedgex,sedgex_lo,sedgex_hi,lo,hi,lo(1)-1,hi(1)+1)
         !
         ! Different stencil needed for x-component of
         ! EXT_DIR and HOEXTRAP bc's.
         !
         if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
            !
            ! The value in the first cc ghost cell
            ! represents the edge value.
            !
            sm(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)     = s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
            sedgex(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                !
                ! Use a modified stencil to get sedgex
                ! on the first interior edge.
                !
                sedgex(lo(1)+1,j,k) = - fifth     *s(lo(1)-1,j,k) &
                                     + three4ths *s(lo(1)  ,j,k) &
                                     + 0.5d0     *s(lo(1)+1,j,k) &
                                     - one20th   *s(lo(1)+2,j,k)
                !
                ! Make sure sedgex lies in between adjacent
                ! cell-centered values.
                !
                sedgex(lo(1)+1,j,k) = max(sedgex(lo(1)+1,j,k),min(s(lo(1)+1,j,k),s(lo(1),j,k)))
                sedgex(lo(1)+1,j,k) = min(sedgex(lo(1)+1,j,k),max(s(lo(1)+1,j,k),s(lo(1),j,k)))
                !
                ! Copy sedgex into sp.
                !
                sp(lo(1)  ,j,k) = sedgex(lo(1)+1,j,k)
             end do
          end do

            !
            ! Apply Colella 2008 limiters to compute
            ! sm and sp in the 2nd and 3rd inner cells.
            !
            call ppm_xdir_colella(s,s_lo,s_hi,sm,sm_lo,sm_hi,sp,sp_lo,sp_hi, &
                sedgex,sedgex_lo,sedgex_hi,lo,hi,lo(1)+1,lo(1)+2)
         end if

         if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
            !
            ! The value in the first cc ghost cell
            ! represents the edge value.
            !
            sp(hi(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

            do k=lo(3)-1,hi(3)+1
               do j=lo(2)-1,hi(2)+1
                  !
                  ! Use a modified stencil to get sedgex
                  ! on the first interior edge.
                  !
                  sedgex(hi(1),j,k) = - fifth     *s(hi(1)+1,j,k) &
                                     + three4ths *s(hi(1)  ,j,k) &
                                     + 0.5d0     *s(hi(1)-1,j,k) &
                                     - one20th   *s(hi(1)-2,j,k)
                  !
                  ! Make sure sedgex lies in between adjacent
                  ! cell-centered values.
                  !
                  sedgex(hi(1),j,k) = max(sedgex(hi(1),j,k),min(s(hi(1)-1,j,k),s(hi(1),j,k)))
                  sedgex(hi(1),j,k) = min(sedgex(hi(1),j,k),max(s(hi(1)-1,j,k),s(hi(1),j,k)))
                  !
                  ! Copy sedgex into sm.
                  !
                  sm(hi(1)  ,j,k) = sedgex(hi(1),j,k)
               end do
            end do
            !
            ! Apply Colella 2008 limiters to compute sm and sp
            ! in the 2nd and 3rd inner cells.
            !
            call ppm_xdir_colella(s,s_lo,s_hi,sm,sm_lo,sm_hi,sp,sp_lo,sp_hi, &
                sedgex,sedgex_lo,sedgex_hi,lo,hi,hi(1)-2,hi(1)-1)
         end if

      end if

      end subroutine ppm_xdir


      subroutine convscalminmax (s,DIMS(s),sn,DIMS(sn), &
                                smin,smax,DIMS(smin),lo,hi,bc)bind(C,name="convscalminmax")
!c
!c     correct an advected field for under/over shoots
!c
      implicit none
      integer  i, j, k, imin, imax, jmin, jmax, kmin, kmax
      integer  DIMDEC(s)
      integer  DIMDEC(sn)
      integer  DIMDEC(smin)
      integer  lo(SDIM), hi(SDIM)
      integer  bc(SDIM,2)
      real(rt)   s(DIMV(s))
      real(rt)   sn(DIMV(sn))
      integer  km, kk, kp
      real(rt)   smn, smx
      real(rt)   smin(DIM12(smin),0:2)
      real(rt)   smax(DIM12(smin),0:2)

      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)

!c     
!c     ::::: compute min/max a slab at a time
!c     ::::: compute min and max of neighbors on kmin-1 slab
!c
      km = 0
      kk = 1
      kp = 2

      k = kmin-1
      do j = jmin, jmax         
         do i = imin, imax
            smin(i,j,km) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k), &
                s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k), &
                s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
            smax(i,j,km) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k), &
                s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k), &
                s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
         end do         
      end do
!c
!c     ::::: compute min and max of neighbors on kmin slab
!c
      k = kmin
      do j = jmin, jmax         
         do i = imin, imax
            smin(i,j,kk) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k), &
                s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k), &
                s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
            smax(i,j,kk) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k), &
                s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k), &
                s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
         end do         
      end do

      do k = kmin, kmax
!c
!c        ::::: compute min and max of neighbors on k+1 slab
!c
         do j = jmin, jmax     
            do i = imin, imax   
               smin(i,j,kp) = min(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1), &
                   s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1), &
                   s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
               smax(i,j,kp) = max(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1), &
                   s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1), &
                   s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
!c
!c        ::::: compute min/max of cell
!c
               smn = min(smin(i,j,km),smin(i,j,kk),smin(i,j,kp))
               smx = max(smax(i,j,km),smax(i,j,kk),smax(i,j,kp))
               sn(i,j,k) = max(sn(i,j,k),smn)
               sn(i,j,k) = min(sn(i,j,k),smx)
               
            end do
         end do
!c
!c        ::::: roll indices for next slab
!c
         km = mod(km+1,3)
         kk = mod(kk+1,3)
         kp = mod(kp+1,3)
      end do

      end subroutine convscalminmax

      subroutine consscalminmax(s,rho,DIMS(s),sn,rhon,DIMS(sn), &
                                smin,smax,DIMS(smin),lo,hi,bc) bind(C,name="consscalminmax")
!c
!c     correct an conservatively-advected field for under/over shoots
!c
      implicit none
      integer  i, j, k, imin, imax, jmin, jmax, kmin, kmax
      integer  DIMDEC(s),DIMDEC(sn)
      integer  DIMDEC(smin)
      integer  lo(SDIM), hi(SDIM)
      integer  bc(SDIM,2)
      real(rt)      s(DIMV(s))
      real(rt)    rho(DIMV(s))
      real(rt)     sn(DIMV(sn))
      real(rt)   rhon(DIMV(sn))
      integer  km, kk, kp
      real(rt)   smn, smx
      real(rt)   smin(DIM12(smin),0:2)
      real(rt)   smax(DIM12(smin),0:2)

      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)

      do k = kmin-1, kmax+1
      do j = jmin-1, jmax+1
         do i = imin-1, imax+1
            s(i,j,k) = s(i,j,k) / rho(i,j,k)
         end do
      end do
      end do
!c
!c     ::::: compute min/max a slab at a time
!c     ::::: compute min and max of neighbors on kmin-1 slab
!c
      km = 0
      kk = 1
      kp = 2

      k = kmin-1
      do j = jmin, jmax         
         do i = imin, imax
            smin(i,j,km) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k), &
                s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k), &
                s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
            smax(i,j,km) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k), &
                s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k), &
                s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k)) 
         end do         
      end do
!c
!c     ::::: compute min and max of neighbors on kmin slab
!c
      k = kmin
      do j = jmin, jmax         
         do i = imin, imax
            smin(i,j,kk) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k), &
                s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k), &
                s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
            smax(i,j,kk) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k), &
                s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k), &
                s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
         end do         
      end do

      do k = kmin, kmax
!c
!c        ::::: compute min and max of neighbors on k+1 slab
!c
         do j = jmin, jmax     
            do i = imin, imax   
               smin(i,j,kp) = min(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1), &
                   s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1), &
                   s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
               smax(i,j,kp) = max(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1), &
                   s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1), &
                   s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
!c
!c        ::::: compute min/max of cell
!c
               smn = min(smin(i,j,km),smin(i,j,kk),smin(i,j,kp))
               smx = max(smax(i,j,km),smax(i,j,kk),smax(i,j,kp))
               sn(i,j,k) = max(sn(i,j,k)/rhon(i,j,k),smn) * rhon(i,j,k)
               sn(i,j,k) = min(sn(i,j,k)/rhon(i,j,k),smx) * rhon(i,j,k)
               
            end do
         end do
!c
!c        ::::: roll indices for next slab
!c
         km = mod(km+1,3)
         kk = mod(kk+1,3)
         kp = mod(kp+1,3)
      end do

      do k = kmin-1, kmax+1
      do j = jmin-1, jmax+1
         do i = imin-1, imax+1
            s(i,j,k) = s(i,j,k) * rho(i,j,k)
         end do
      end do
      end do

      end subroutine consscalminmax

      subroutine fort_sum_tf_gp( &
          tforces,DIMS(tf), &
          gp,DIMS(gp), &
          rho,DIMS(rho), &
          lo,hi ) bind(C,name="fort_sum_tf_gp")

!c
!c     sum pressure forcing into tforces
!c
      implicit none
      integer i, j, k, n
      integer DIMDEC(tf)
      integer DIMDEC(gp)
      integer DIMDEC(rho)
      integer lo(SDIM), hi(SDIM)
      real(rt) tforces(DIMV(tf),SDIM)
      real(rt) gp(DIMV(gp),SDIM)
      real(rt) rho(DIMV(rho))
      real(rt), allocatable :: irho(:,:,:)

      allocate(irho(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               irho(i,j,k) = 1.0d0/rho(i,j,k)
            end do
         end do
      end do

      do n = 1, SDIM
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,k,n) = (tforces(i,j,k,n) - gp(i,j,k,n))*irho(i,j,k)
               end do
            end do
         end do
      end do
      end subroutine fort_sum_tf_gp

      subroutine fort_sum_tf_gp_visc( &
          tforces,DIMS(tf), &
          visc,DIMS(visc), &
          gp,DIMS(gp), &
          rho,DIMS(rho), &
          lo,hi ) bind(C,name="fort_sum_tf_gp_visc")
!c
!c     sum pressure forcing and viscous forcing into
!c     tforces
!c
      implicit none
      integer i, j, k, n
      integer DIMDEC(tf)
      integer DIMDEC(visc)
      integer DIMDEC(gp)
      integer DIMDEC(rho)
      integer lo(SDIM), hi(SDIM)
      real(rt) tforces(DIMV(tf),SDIM)
      real(rt) visc(DIMV(visc),SDIM)
      real(rt) gp(DIMV(gp),SDIM)
      real(rt) rho(DIMV(rho))
      real(rt), allocatable :: irho(:,:,:)

      allocate(irho(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               irho(i,j,k) = 1.0d0/rho(i,j,k)
            end do
         end do
      end do

      do n = 1, SDIM
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,k,n) = (tforces(i,j,k,n) + visc(i,j,k,n) &
                      -    gp(i,j,k,n) )*irho(i,j,k)
               end do
            end do
         end do
      end do
      end subroutine fort_sum_tf_gp_visc

      subroutine fort_sum_tf_divu( &
          S,DIMS(S), &
          tforces,DIMS(tf), &
          divu,DIMS(divu), &
          rho,DIMS(rho), &
          lo,hi,nvar,iconserv ) bind(C,name="fort_sum_tf_divu")
!c
!c     sum tforces, viscous forcing and divU*S into tforces
!c     depending on the value of iconserv
!c
      implicit none
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, k, n

      integer DIMDEC(S)
      integer DIMDEC(tf)
      integer DIMDEC(divu)
      integer DIMDEC(rho)

      real(rt) S(DIMV(S),nvar)
      real(rt) tforces(DIMV(tf),nvar)
      real(rt) divu(DIMV(divu))
      real(rt) rho(DIMV(rho))

      if ( iconserv .eq. 1 ) then
         do n = 1, nvar
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) =  &
                    tforces(i,j,k,n) - S(i,j,k,n)*divu(i,j,k)
                  end do
               end do
            end do
         end do
      else
         do n = 1, nvar
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = tforces(i,j,k,n)/rho(i,j,k)
                  end do
               end do
            end do
         end do
      end if

      end subroutine fort_sum_tf_divu

      subroutine fort_sum_tf_divu_visc( &
          S,DIMS(S), &
          tforces,DIMS(tf), &
          divu,DIMS(divu), &
          visc,DIMS(visc), &
          rho,DIMS(rho), &
          lo,hi,nvar,iconserv ) bind(C,name="fort_sum_tf_divu_visc")
!c
!c     sum tforces, viscous forcing and divU*S into tforces
!c     depending on the value of iconserv
!c
      implicit none
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, k, n

      integer DIMDEC(S)
      integer DIMDEC(tf)
      integer DIMDEC(divu)
      integer DIMDEC(visc)
      integer DIMDEC(rho)

      real(rt) S(DIMV(S),nvar)
      real(rt) tforces(DIMV(tf),nvar)
      real(rt) divu(DIMV(divu))
      real(rt) visc(DIMV(visc),nvar)
      real(rt) rho(DIMV(rho))


      if ( iconserv .eq. 1 ) then
         do n = 1, nvar
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = tforces(i,j,k,n) +  visc(i,j,k,n) &
                         - S(i,j,k,n)*divu(i,j,k)
                  end do
               end do
            end do
         end do
      else
         do n = 1, nvar
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = (tforces(i,j,k,n) + visc(i,j,k,n))/rho(i,j,k)
                  end do
               end do
            end do
         end do
      end if
      
      end subroutine fort_sum_tf_divu_visc

      subroutine update_tf ( &
          s,       DIMS(s), &
          sn,      DIMS(sn), &
          tforces, DIMS(tf), &
          lo,hi,dt,nvar) bind(C,name="update_tf")
!c
!c     update a field with a forcing term
!c
      implicit none
      integer i, j, k, n, nvar
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(tf)
      integer lo(SDIM), hi(SDIM)
      real(rt) dt
      real(rt) s(DIMV(s),nvar)
      real(rt) sn(DIMV(sn),nvar)
      real(rt) tforces(DIMV(tf),nvar)

      do n = 1,nvar
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  sn(i,j,k,n) = s(i,j,k,n) + dt*tforces(i,j,k,n)
               end do
            end do
         end do
      end do

      end subroutine update_tf 

      subroutine update_aofs_tf ( &
          s,       DIMS(s), &
          sn,      DIMS(sn), &
          aofs,    DIMS(aofs), &
          tforces, DIMS(tf), &
          lo,hi,dt,nvar) bind(C,name="update_aofs_tf")
!c
!c     update a field with an advective tendency
!c     and a forcing term
!c
      implicit none
      integer i, j, k, n, nvar
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(aofs)
      integer DIMDEC(tf)
      integer lo(SDIM), hi(SDIM)
      real(rt) dt
      real(rt) s(DIMV(s),nvar)
      real(rt) sn(DIMV(sn),nvar)
      real(rt) aofs(DIMV(aofs),nvar)
      real(rt) tforces(DIMV(tf),nvar)

      do n = 1,nvar
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  sn(i,j,k,n) = s(i,j,k,n) + dt*(tforces(i,j,k,n) - aofs(i,j,k,n))
               end do
            end do
         end do
      end do
      end subroutine update_aofs_tf 

      subroutine update_aofs_tf_gp ( &
          u,       DIMS(u), &
          un,      DIMS(un), &
          aofs,    DIMS(aofs), &
          tforces, DIMS(tf), &
          gp,      DIMS(gp), &
          rho,     DIMS(rho), &
          lo, hi, dt) bind(C,name="update_aofs_tf_gp")
      !
      ! update the velocities
      !
      implicit none
      integer i, j, k, n
      integer DIMDEC(u)
      integer DIMDEC(un)
      integer DIMDEC(aofs)
      integer DIMDEC(rho)
      integer DIMDEC(gp)
      integer DIMDEC(tf)
      integer lo(SDIM), hi(SDIM)
      real(rt) u(DIMV(u),SDIM)
      real(rt) un(DIMV(un),SDIM)
      real(rt) aofs(DIMV(aofs),SDIM)
      real(rt) rho(DIMV(rho))
      real(rt) gp(DIMV(gp),SDIM)
      real(rt) tforces(DIMV(tf),SDIM)
      real(rt) dt

      do n = 1, SDIM
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  un(i,j,k,n) = u(i,j,k,n) + dt * &
                      ( (tforces(i,j,k,n) - gp(i,j,k,n)) / rho(i,j,k) - aofs(i,j,k,n) )
                  
               end do
            end do
         end do
      end do

      end subroutine update_aofs_tf_gp
      
 end module godunov_3d_module
