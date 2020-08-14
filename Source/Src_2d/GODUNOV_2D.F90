
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <GODUNOV_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2
#define XVEL 1
#define YVEL 2


module godunov_2d_module

  use amrex_fort_module, only : rt=>amrex_real

  implicit none

  private

  public :: extrap_vel_to_faces, fort_maxchng_velmag, &
       fort_test_umac_rho, adv_forcing, &
       sync_adv_forcing, convscalminmax,consscalminmax,&
       fort_sum_tf_gp,fort_sum_tf_gp_visc,fort_sum_tf_divu,&
       fort_sum_tf_divu_visc, update_tf,&
       update_aofs_tf, update_aofs_tf_gp

contains

  subroutine extrap_vel_to_faces(lo,hi,&
       u,u_lo,u_hi,&
       ubc, tfx,tfx_lo,tfx_hi, umac,umac_lo,umac_hi,&
       vbc, tfy,tfy_lo,tfy_hi, vmac,vmac_lo,vmac_hi,&
       dt, dx, use_forces_in_trans, ppm_type)  bind(C,name="extrap_vel_to_faces")

    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate

    implicit none
    real(rt), intent(in) :: dt, dx(SDIM)
    integer,  intent(in) ::  ubc(SDIM,2),vbc(SDIM,2), use_forces_in_trans, ppm_type
    integer,  intent(in), dimension(2) :: lo,hi,u_lo,u_hi,&
         tfx_lo,tfx_hi,tfy_lo,tfy_hi,umac_lo,umac_hi,vmac_lo,vmac_hi

    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),2)
    real(rt), intent(in) :: tfx(tfx_lo(1):tfx_hi(1),tfx_lo(2):tfx_hi(2))
    real(rt), intent(inout) :: umac(umac_lo(1):umac_hi(1),umac_lo(2):umac_hi(2)) ! result
    real(rt), intent(in) :: tfy(tfy_lo(1):tfy_hi(1),tfy_lo(2):tfy_hi(2))
    real(rt), intent(inout) :: vmac(vmac_lo(1):vmac_hi(1),vmac_lo(2):vmac_hi(2)) ! result

    integer, dimension(2) :: wklo,wkhi,uwlo,uwhi,vwlo,vwhi,eblo,ebhi,ebxhi,ebyhi,g2lo,g2hi
    real(rt), dimension(:,:), pointer, contiguous :: xlo,xhi,sx,uad
    real(rt), dimension(:,:), pointer, contiguous :: ylo,yhi,sy,vad
    real(rt), dimension(:,:), pointer, contiguous :: Imx,Ipx,sedgex
    real(rt), dimension(:,:), pointer, contiguous :: Imy,Ipy,sedgey
    real(rt), dimension(:,:), pointer, contiguous :: sm,sp,dsvl
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
    call amrex_allocate(xlo,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(xhi,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(sx,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(ylo,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(yhi,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(sy,wklo(1),wkhi(1),wklo(2),wkhi(2))

    uwlo = wklo
    uwhi = wkhi
    uwhi(1) = wkhi(1) + 1
    vwlo = wklo
    vwhi = wkhi
    vwhi(2) = vwhi(2) + 1

    call amrex_allocate(uad,uwlo(1),uwhi(1),uwlo(2),uwhi(2))
    call amrex_allocate(vad,vwlo(1),vwhi(1),vwlo(2),vwhi(2))

    ! Note: Intel unhappy to pass pointers through subroutines if not allocated
    !     We just allocate something small (and still promise not to use it)
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
    call amrex_allocate(Imx,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(Ipx,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(Imy,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(Ipy,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(sm,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(sp,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(sedgex,eblo(1),ebxhi(1),eblo(2),ebxhi(2))
    call amrex_allocate(sedgey,eblo(1),ebyhi(1),eblo(2),ebyhi(2))
    g2lo = lo - 2
    g2hi = hi + 2
    call amrex_allocate(dsvl,g2lo(1),g2hi(1),g2lo(2),g2hi(2))

    ! get velocities that resolve upwind directions on faces used to compute transverse derivatives (uad,vad)
    ! Note that this is done only in this "pre-mac" situation, to get velocities on faces that will be projected.
    ! When advecting the other states, we will use the projected velocities, not these approximate versions.
    ! These face-centered arrays for each direction, dir, are needed on surroundingNodes(gBox(lo,hi),dir), where
    ! gBox is Box(lo,hi) grown in all directions but dir.
    call transvel(lo, hi,&
         u(u_lo(1),u_lo(2),1),u_lo,u_hi,&
         uad,uwlo,uwhi,&
         xhi,wklo,wkhi,&
         sx,wklo,wkhi,&
         ubc,&
         Imx,wklo,wkhi,&
         Ipx,wklo,wkhi,&
         sedgex,eblo,ebxhi,&
         u(u_lo(1),u_lo(2),2),u_lo,u_hi,&
         vad,vwlo,vwhi,&
         yhi,wklo,wkhi,&
         sy,wklo,wkhi,&
         vbc,&
         Imy,wklo,wkhi,&
         Ipy,wklo,wkhi,&
         sedgey,eblo,ebyhi,&
         dsvl,g2lo,g2hi,&
         sm,wklo,wkhi,&
         sp,wklo,wkhi,&
         tfx,tfx_lo,tfx_hi,&
         tfy,tfy_lo,tfy_hi,&
         dt, dx, use_forces_in_trans, ppm_type)

    ! Note that the final edge velocites are are resolved using the average of the velocities predicted from both
    ! sides (including the transverse terms resolved using uad,vad from above).

    ! get velocity on x-face, predict from cc u (because velpred=1 && n=XVEL, only predicts to xfaces)
    call estate_premac(lo,hi,&
         u(u_lo(1),u_lo(2),1),u_lo,u_hi,&
         tfx,tfx_lo,tfx_hi,&
         u(u_lo(1),u_lo(2),1),u_lo,u_hi,&
         xlo,wklo,wkhi,&
         xhi,wklo,wkhi,&
         sx,wklo,wkhi,&
         uad,uwlo,uwhi,&
         uad,uwlo,uwhi,&  ! Unused since velpred = 1
         umac,umac_lo,umac_hi,&
         Imx,wklo,wkhi,&
         Ipx,wklo,wkhi,&
         sedgex,eblo,ebxhi,&
         u(u_lo(1),u_lo(2),2),u_lo,u_hi,&
         ylo,wklo,wkhi,&
         yhi,wklo,wkhi,&
         sy,wklo,wkhi,&
         vad,vwlo,vwhi,&
         vad,vwlo,vwhi,&  ! Unused since velpred = 1
         vmac,vmac_lo,vmac_hi,&
         Imy,wklo,wkhi,&
         Ipy,wklo,wkhi,&
         sedgey,eblo,ebyhi,&
         dsvl,g2lo,g2hi,&
         sm,wklo,wkhi,&
         sp,wklo,wkhi,&
         ubc, dt, dx, XVEL, 1, velpred, use_forces_in_trans, ppm_type)

    ! get velocity on y-face, predict from cc v (because velpred=1 && n=YVEL, only predicts to yfaces)
    call estate_premac(lo,hi,&
         u(u_lo(1),u_lo(2),2),u_lo,u_hi,&
         tfy,tfy_lo,tfy_hi,&
         u(u_lo(1),u_lo(2),1),u_lo,u_hi,&
         xlo,wklo,wkhi,&
         xhi,wklo,wkhi,&
         sx,wklo,wkhi,&
         uad,uwlo,uwhi,&
         uad,uwlo,uwhi,& ! Unused since velpred = 1
         umac,umac_lo,umac_hi,&
         Imx,wklo,wkhi,&
         Ipx,wklo,wkhi,&
         sedgex,eblo,ebxhi,&
         u(u_lo(1),u_lo(2),2),u_lo,u_hi,&
         ylo,wklo,wkhi,&
         yhi,wklo,wkhi,&
         sy,wklo,wkhi,&
         vad,vwlo,vwhi,&
         vad,vwlo,vwhi,& ! Unused since velpred = 1
         vmac,vmac_lo,vmac_hi,&
         Imy,wklo,wkhi,&
         Ipy,wklo,wkhi,&
         sedgey,eblo,ebyhi,&
         dsvl,g2lo,g2hi,&
         sm,wklo,wkhi,&
         sp,wklo,wkhi,&
         vbc, dt, dx, YVEL, 1, velpred, use_forces_in_trans, ppm_type)

    call amrex_deallocate(xlo)
    call amrex_deallocate(xhi)
    call amrex_deallocate(sx)
    call amrex_deallocate(ylo)
    call amrex_deallocate(yhi)
    call amrex_deallocate(sy)
    call amrex_deallocate(uad)
    call amrex_deallocate(vad)
    call amrex_deallocate(Imx)
    call amrex_deallocate(Ipx)
    call amrex_deallocate(Imy)
    call amrex_deallocate(Ipy)
    call amrex_deallocate(sm)
    call amrex_deallocate(sp)
    call amrex_deallocate(sedgex)
    call amrex_deallocate(sedgey)
    call amrex_deallocate(dsvl)

  end subroutine extrap_vel_to_faces

  subroutine extrap_state_to_faces(lo,hi,&
       s,s_lo,s_hi,nc,              tf, tf_lo,tf_hi,              divu,divu_lo,divu_hi,&
       umac,umac_lo,umac_hi,        xstate,xstate_lo,xstate_hi,&
       vmac,vmac_lo,vmac_hi,        ystate,ystate_lo,ystate_hi,&
       dt, dx, bc, state_ind, use_forces_in_trans, ppm_type, iconserv)  bind(C,name="extrap_state_to_faces")

    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate

    implicit none
    integer, intent(in) ::  nc, bc(SDIM,2,nc), state_ind, use_forces_in_trans, ppm_type, iconserv(nc)
    integer, dimension(2), intent(in) :: lo,hi,s_lo,s_hi,tf_lo,tf_hi,&
         divu_lo,divu_hi,xstate_lo,xstate_hi,ystate_lo,ystate_hi,umac_lo,umac_hi,vmac_lo,vmac_hi

    real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),nc)
    real(rt), intent(in) :: tf(tf_lo(1):tf_hi(1),tf_lo(2):tf_hi(2),nc)
    real(rt), intent(in) :: divu(divu_lo(1):divu_hi(1),divu_lo(2):divu_hi(2))
    real(rt), intent(in) :: umac(umac_lo(1):umac_hi(1),umac_lo(2):umac_hi(2))
    real(rt), intent(inout) :: xstate(xstate_lo(1):xstate_hi(1),xstate_lo(2):xstate_hi(2),nc) ! result
    real(rt), intent(in) :: vmac(vmac_lo(1):vmac_hi(1),vmac_lo(2):vmac_hi(2))
    real(rt), intent(inout) :: ystate(ystate_lo(1):ystate_hi(1),ystate_lo(2):ystate_hi(2),nc) ! result

    real(rt), intent(in) :: dt, dx(SDIM)

    integer, dimension(2) :: wklo,wkhi,eblo,ebhi,ebxhi,ebyhi,g2lo,g2hi
    real(rt), dimension(:,:), pointer, contiguous :: xlo,xhi,sx
    real(rt), dimension(:,:), pointer, contiguous :: ylo,yhi,sy
    real(rt), dimension(:,:), pointer, contiguous :: Imx,Ipx,sedgex
    real(rt), dimension(:,:), pointer, contiguous :: Imy,Ipy,sedgey
    real(rt), dimension(:,:), pointer, contiguous :: sm,sp,dsvl

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
    call amrex_allocate(xlo,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(xhi,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(sx,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(ylo,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(yhi,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(sy,wklo(1),wkhi(1),wklo(2),wkhi(2))

    ! Note: Intel unhappy to pass pointers through subroutines if not allocated
    !     We just allocate something small (and still promise not to use it)
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
    call amrex_allocate(Imx,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(Ipx,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(Imy,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(Ipy,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(sm,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(sp,wklo(1),wkhi(1),wklo(2),wkhi(2))
    call amrex_allocate(sedgex,eblo(1),ebxhi(1),eblo(2),ebxhi(2))
    call amrex_allocate(sedgey,eblo(1),ebyhi(1),eblo(2),ebyhi(2))
    g2lo = lo - 2
    g2hi = hi + 2
    call amrex_allocate(dsvl,g2lo(1),g2hi(1),g2lo(2),g2hi(2))

    call estate_fpu(lo,hi,&
         s,s_lo,s_hi,&
         tf,tf_lo,tf_hi,&
         divu,divu_lo,divu_hi,&
         xlo,wklo,wkhi,&
         xhi,wklo,wkhi,&
         sx,wklo,wkhi,&
         umac,umac_lo,umac_hi,&
         xstate,xstate_lo,xstate_hi,&
         Imx,wklo,wkhi,&
         Ipx,wklo,wkhi,&
         sedgex,eblo,ebxhi,&
         ylo,wklo,wkhi,&
         yhi,wklo,wkhi,&
         sy,wklo,wkhi,&
         vmac,vmac_lo,vmac_hi,&
         ystate,ystate_lo,ystate_hi,&
         Imy,wklo,wkhi,&
         Ipy,wklo,wkhi,&
         sedgey,eblo,ebyhi,&
         dsvl,g2lo,g2hi,&
         sm,wklo,wkhi,&
         sp,wklo,wkhi,&
         bc, dt, dx, state_ind, nc, use_forces_in_trans, iconserv, ppm_type)

    call amrex_deallocate(xlo)
    call amrex_deallocate(xhi)
    call amrex_deallocate(sx)
    call amrex_deallocate(ylo)
    call amrex_deallocate(yhi)
    call amrex_deallocate(sy)
    call amrex_deallocate(Imx)
    call amrex_deallocate(Ipx)
    call amrex_deallocate(Imy)
    call amrex_deallocate(Ipy)
    call amrex_deallocate(sm)
    call amrex_deallocate(sp)
    call amrex_deallocate(sedgex)
    call amrex_deallocate(sedgey)
    call amrex_deallocate(dsvl)

  end subroutine extrap_state_to_faces


    subroutine fort_maxchng_velmag (&
          old_vel,DIMS(old_vel),&
          new_vel,DIMS(new_vel),&
          lo,hi,max_change) bind(C,name="fort_maxchng_velmag")
! c
! c     ----------------------------------------------------------
! c     Given the velocity field at the previous and current time steps
! c     (old_vel and new_vel, respectively), find the largest change in
! c     velocity magnitude between the two.
! c     ----------------------------------------------------------
! c
      implicit none
      REAL_T   old_velmag, new_velmag
      integer  i, j
      integer  lo(SDIM), hi(SDIM)
      REAL_T   max_change

      integer DIMDEC(old_vel)
      integer DIMDEC(new_vel)

      REAL_T  old_vel(DIMV(old_vel),SDIM)
      REAL_T  new_vel(DIMV(new_vel),SDIM)

      max_change = zero

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            old_velmag = sqrt(old_vel(i,j,1)**2 + old_vel(i,j,2)**2)
            new_velmag = sqrt(new_vel(i,j,1)**2 + new_vel(i,j,2)**2)
            max_change = max(max_change, abs(new_velmag - old_velmag))
         end do
      end do

    end subroutine fort_maxchng_velmag

      subroutine fort_test_umac_rho(&
          umac,DIMS(umac),&
          vmac,DIMS(vmac),&
          rho,DIMS(rho),&
          lo,hi,dt,dx,cflmax,u_max) bind(C,name="fort_test_umac_rho")
! c
! c     This subroutine computes the extrema of the density
! c     and mac edge velocities
! c
      implicit none
      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(rho)
      integer imin, imax, jmin, jmax
      integer i, j
      integer lo(SDIM),hi(SDIM)
      REAL_T  dt, dx(SDIM), u_max(SDIM), cflmax
      REAL_T  hx, hy
      REAL_T  umax, vmax, rhomax
      REAL_T  umin, vmin, rhomin
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      REAL_T  rho(DIMV(rho))

      hx   = dx(1)
      hy   = dx(2)
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      umax = -1.d200
      vmax = -1.d200
      umin =  1.d200
      vmin =  1.d200
      rhomax = -1.d200
      rhomin =  1.d200

      do j = jmin, jmax
         do i = imin, imax+1
            umax = max(umax,umac(i,j))
            umin = min(umin,umac(i,j))
         end do
      end do
      do j = jmin, jmax+1
         do i = imin, imax
            vmax = max(vmax,vmac(i,j))
            vmin = min(vmin,vmac(i,j))
         end do
      end do
      do j = jmin, jmax
         do i = imin, imax
            rhomax = max(rhomax,rho(i,j))
            rhomin = min(rhomin,rho(i,j))
         end do
      end do

      u_max(1) = max(abs(umax), abs(umin))
      u_max(2) = max(abs(vmax), abs(vmin))
      cflmax   = dt*max(u_max(1)/hx,u_max(2)/hy)

      write(6,1000) umax,umin,u_max(1)
      write(6,1001) vmax,vmin,u_max(2)
      write(6,1002) rhomax,rhomin

 1000 format('UMAC MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1001 format('VMAC MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1002 format('RHO  MAX/MIN      ',e21.14,2x,e21.14)

#ifndef	BL_NO_FLUSH
!c      call flush(6)
#endif

    end subroutine fort_test_umac_rho

    subroutine transvel(lo,hi,&
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
           dsvl,dsvl_lo,dsvl_hi,&
           sm,sm_lo,sm_hi,&
           sp,sp_lo,sp_hi,&
           tfx,tfx_lo,tfx_hi,&
           tfy,tfy_lo,tfy_hi,&
           dt, dx, use_minion, ppm_type)
! c
! c     This subroutine computes the advective velocities used in
! c     the transverse derivatives of the Godunov box
! c
      implicit none
      integer, intent(in) ::  ubc(SDIM,2),vbc(SDIM,2), use_minion, ppm_type
      integer, dimension(2), intent(in) :: lo,hi,&
           u_lo,u_hi,ulo_lo,ulo_hi,uhi_lo,uhi_hi,sx_lo,sx_hi,&
           Imx_lo,Imx_hi,Ipx_lo,Ipx_hi,sedgex_lo,sedgex_hi,&
           v_lo,v_hi,vlo_lo,vlo_hi,vhi_lo,vhi_hi,sy_lo,sy_hi,&
           Imy_lo,Imy_hi,Ipy_lo,Ipy_hi,sedgey_lo,sedgey_hi,&
           dsvl_lo,dsvl_hi,sm_lo,sm_hi,sp_lo,sp_hi,tfx_lo,tfx_hi,tfy_lo,tfy_hi

      real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2))
      real(rt), intent(inout) :: ulo(ulo_lo(1):ulo_hi(1),ulo_lo(2):ulo_hi(2))
      real(rt), intent(inout) :: uhi(uhi_lo(1):uhi_hi(1),uhi_lo(2):uhi_hi(2))
      real(rt), intent(inout) :: sx(sx_lo(1):sx_hi(1),sx_lo(2):sx_hi(2))
      real(rt), intent(inout) :: Imx(Imx_lo(1):Imx_hi(1),Imx_lo(2):Imx_hi(2))
      real(rt), intent(inout) :: Ipx(Ipx_lo(1):Ipx_hi(1),Ipx_lo(2):Ipx_hi(2))
      real(rt), intent(inout) :: sedgex(sedgex_lo(1):sedgex_hi(1),sedgex_lo(2):sedgex_hi(2))
      real(rt), intent(in) :: tfx(tfx_lo(1):tfx_hi(1),tfx_lo(2):tfx_hi(2))

      real(rt), intent(in) :: v(v_lo(1):v_hi(1),v_lo(2):v_hi(2))
      real(rt), intent(inout) :: vlo(vlo_lo(1):vlo_hi(1),vlo_lo(2):vlo_hi(2))
      real(rt), intent(inout) :: vhi(vhi_lo(1):vhi_hi(1),vhi_lo(2):vhi_hi(2))
      real(rt), intent(inout) :: sy(sy_lo(1):sy_hi(1),sy_lo(2):sy_hi(2))
      real(rt), intent(inout) :: Imy(Imy_lo(1):Imy_hi(1),Imy_lo(2):Imy_hi(2))
      real(rt), intent(inout) :: Ipy(Ipy_lo(1):Ipy_hi(1),Ipy_lo(2):Ipy_hi(2))
      real(rt), intent(inout) :: sedgey(sedgey_lo(1):sedgey_hi(1),sedgey_lo(2):sedgey_hi(2))
      real(rt), intent(in) :: tfy(tfy_lo(1):tfy_hi(1),tfy_lo(2):tfy_hi(2))

      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2))

      integer :: i,j, imin,jmin,imax,jmax
      real(rt) :: hx, hy, dt, dth, dthx, dthy, dx(SDIM), uad, vad
      real(rt) :: eps,eps_for_bc
      logical :: ltm
      parameter( eps        = 1.0D-6 )
      parameter( eps_for_bc = 1.0D-10 )

      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      hx   = dx(1)
      hy   = dx(2)
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

! c
! c     --------------------------------------------------------------
! c     compute the x transverse velocities
! c     --------------------------------------------------------------
! c
      if (ppm_type .gt. 0) then
         call ppm(lo,hi,&
              u,u_lo,u_hi,&
              u,u_lo,u_hi,&
              v,v_lo,v_hi,&
              Ipx,Ipx_lo,Ipx_hi,&
              Imx,Imx_lo,Imx_hi,&
              Ipy,Ipy_lo,Ipy_hi,&
              Imy,Imy_lo,Imy_hi,&
              sm,sm_lo,sm_hi,&
              sp,sp_lo,sp_hi,&
              dsvl,dsvl_lo,dsvl_hi,&
              sedgex,sedgex_lo,sedgex_hi,&
              sedgey,sedgey_lo,sedgey_hi,&
              dx, dt, ubc, eps_for_bc, ppm_type)
      else
         call slopes(lo,hi,&
              u,u_lo,u_hi,&
              sx,sx_lo,sx_hi,&
              sy,sy_lo,sy_hi,&
              ubc)
      end if

      if (ppm_type .gt. 0) then
         do i = imin, imax+1
            do j = jmin-1,jmax+1
               ulo(i,j) = Ipx(i-1,j)
               uhi(i,j) = Imx(i  ,j)
            end do
         end do
      else
         do i = imin, imax+1
            do j = jmin-1,jmax+1
               ulo(i,j) = u(i-1,j) + (half  - dthx*u(i-1,j))*sx(i-1,j)
               uhi(i,j) = u(i,j)   + (-half - dthx*u(i,  j))*sx(i,j)
            end do
         end do
      end if

      if(use_minion.eq.1)then
        do i = imin, imax+1
          do j = jmin-1,jmax+1
            ulo(i,j) = ulo(i,j) + dth*tfx(i-1,j)
            uhi(i,j) = uhi(i,j) + dth*tfx(i,  j)
          end do
        end do
      end if

      call trans_xbc(lo,hi,&
           u,u_lo,u_hi,&
           ulo,ulo_lo,ulo_hi,&
           uhi,uhi_lo,uhi_hi,&
           ulo,ulo_lo,ulo_hi,&
           XVEL, ubc, eps_for_bc)

      do i = imin,imax+1
         do j = jmin-1,jmax+1
            uad = merge(ulo(i,j),uhi(i,j),(ulo(i,j)+uhi(i,j)) .ge. 0.0d0)
            ltm = ulo(i,j) .le. 0.d0  .and.  uhi(i,j) .ge. 0.d0
            ltm = ltm .or. (abs(ulo(i,j)+uhi(i,j)) .lt. eps)
            ulo(i,j) = merge(0.d0,uad,ltm)
         end do
      end do
! c
! c     compute the y transverse velocities
! c
      if (ppm_type .gt. 0) then
         call ppm(lo,hi,&
              v,v_lo,v_hi,&
              u,u_lo,u_hi,&
              v,v_lo,v_hi,&
              Ipx,Ipx_lo,Ipx_hi,&
              Imx,Imx_lo,Imx_hi,&
              Ipy,Ipy_lo,Ipy_hi,&
              Imy,Imy_lo,Imy_hi,&
              sm,sm_lo,sm_hi,&
              sp,sp_lo,sp_hi,&
              dsvl,dsvl_lo,dsvl_hi,&
              sedgex,sedgex_lo,sedgex_hi,&
              sedgey,sedgey_lo,sedgey_hi,&
              dx, dt, vbc, eps_for_bc,ppm_type)
      else
         call slopes(lo,hi,&
              v,v_lo,v_hi,&
              sx,sx_lo,sx_hi,&
              sy,sy_lo,sy_hi,&
              vbc)
      end if

      if (ppm_type .gt. 0) then
         do i = imin-1,imax+1
            do j = jmin,jmax+1
               vlo(i,j) = Ipy(i,j-1)
               vhi(i,j) = Imy(i,j  )
            end do
         end do
      else
         do i = imin-1,imax+1
            do j = jmin,jmax+1
               vlo(i,j) = v(i,j-1) + (half  - dthy*v(i,j-1))*sy(i,j-1)
               vhi(i,j) = v(i,j)   + (-half - dthy*v(i,j  ))*sy(i,j)
            end do
         end do
      end if

      if(use_minion.eq.1)then
         do i = imin-1, imax+1
            do j = jmin,jmax+1
               vlo(i,j) = vlo(i,j) + dth*tfy(i,j-1)
               vhi(i,j) = vhi(i,j) + dth*tfy(i,j)
            end do
         end do
      end if

      call trans_ybc(lo,hi,&
           v,v_lo,v_hi,&
           vlo,vlo_lo,vlo_hi,&
           vhi,vhi_lo,vhi_hi,&
           vlo,vlo_lo,vlo_hi,&
           YVEL, vbc, eps_for_bc)

      do i = imin-1,imax+1
         do j = jmin,jmax+1
            vad = merge(vlo(i,j),vhi(i,j),(vlo(i,j)+vhi(i,j)) .ge. 0.0d0)
            ltm = vlo(i,j) .le. 0.d0  .and.  vhi(i,j) .ge. 0.d0
            ltm = ltm .or. (abs(vlo(i,j)+vhi(i,j)) .lt. eps)
            vlo(i,j) = merge(0.d0,vad,ltm)
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
         vedge,vedge_lo,vedge_hi,&
         ystate,ystate_lo,ystate_hi,&
         Imy,Imy_lo,Imy_hi,&
         Ipy,Ipy_lo,Ipy_hi,&
         sedgey,sedgey_lo,sedgey_hi,&
         dsvl,dsvl_lo,dsvl_hi,&
         sm,sm_lo,sm_hi,&
         sp,sp_lo,sp_hi,&
         bc, dt, dx, n, nc, velpred, use_minion, ppm_type)

      ! Here is what this routine does
      ! 1. Trace values of state from cell-centers to faces using cell-centered velocities
      !     and excluding transverse corrections.  This produces data in xlo,xho.  Enforce
      !     boundary conditions on faces, then resolve the upwind value from xlo,xhi using
      !     uad and put this into xlo.  If uad is small, xlo=(xlo+xhi)/2
      ! 1a. Repeat for ylo
      ! 2. The transverse correction is - vbar * grady, where vbar=(vad_{j+1}+vad_{j})/2
      !      and grady is computed with the ylo values.  Results into stxlo,stxhi, and
      !      upwinded into xstate using UFACE, where
      !      if velpred!=1:
      !          UFACE = uedge
      !      else
      !          UFACE = stxlo + stxhi
      ! 2a. Repeat for transverse corrections to get stylo,styhi, and resolve into ystate
      !
      implicit none

      integer, intent(in) :: velpred, use_minion, ppm_type, bc(SDIM,2), n, nc
      real(rt), intent(in) :: dt, dx(SDIM)

      integer, dimension(2), intent(in) :: s_lo,s_hi,tf_lo,tf_hi,&
           u_lo,u_hi,xlo_lo,xlo_hi,xhi_lo,xhi_hi,sx_lo,sx_hi,uad_lo,uad_hi,&
           uedge_lo,uedge_hi,xstate_lo,xstate_hi,Imx_lo,Imx_hi,Ipx_lo,Ipx_hi,sedgex_lo,sedgex_hi,&
           v_lo,v_hi,ylo_lo,ylo_hi,yhi_lo,yhi_hi,sy_lo,sy_hi,vad_lo,vad_hi,&
           vedge_lo,vedge_hi,ystate_lo,ystate_hi,Imy_lo,Imy_hi,Ipy_lo,Ipy_hi,sedgey_lo,sedgey_hi,&
           dsvl_lo,dsvl_hi,sm_lo,sm_hi,sp_lo,sp_hi,lo,hi

      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),nc)
      real(rt), intent(in) :: tf(tf_lo(1):tf_hi(1),tf_lo(2):tf_hi(2),nc)
      real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2))
      real(rt), intent(inout) :: xlo(xlo_lo(1):xlo_hi(1),xlo_lo(2):xlo_hi(2))
      real(rt), intent(inout) :: xhi(xhi_lo(1):xhi_hi(1),xhi_lo(2):xhi_hi(2))
      real(rt), intent(inout) :: sx(sx_lo(1):sx_hi(1),sx_lo(2):sx_hi(2))
      real(rt), intent(in) :: uad(uad_lo(1):uad_hi(1),uad_lo(2):uad_hi(2))
      real(rt), intent(in) :: uedge(uedge_lo(1):uedge_hi(1),uedge_lo(2):uedge_hi(2))
      real(rt), intent(inout) :: xstate(xstate_lo(1):xstate_hi(1),xstate_lo(2):xstate_hi(2),nc) ! result
      real(rt), intent(inout) :: Imx(Imx_lo(1):Imx_hi(1),Imx_lo(2):Imx_hi(2))
      real(rt), intent(inout) :: Ipx(Ipx_lo(1):Ipx_hi(1),Ipx_lo(2):Ipx_hi(2))
      real(rt), intent(inout) :: sedgex(sedgex_lo(1):sedgex_hi(1),sedgex_lo(2):sedgex_hi(2))
      real(rt), intent(in) :: v(v_lo(1):v_hi(1),v_lo(2):v_hi(2))
      real(rt), intent(inout) :: ylo(ylo_lo(1):ylo_hi(1),ylo_lo(2):ylo_hi(2))
      real(rt), intent(inout) :: yhi(yhi_lo(1):yhi_hi(1),yhi_lo(2):yhi_hi(2))
      real(rt), intent(inout) :: sy(sy_lo(1):sy_hi(1),sy_lo(2):sy_hi(2))
      real(rt), intent(in) :: vad(vad_lo(1):vad_hi(1),vad_lo(2):vad_hi(2))
      real(rt), intent(in) :: vedge(vedge_lo(1):vedge_hi(1),vedge_lo(2):vedge_hi(2))
      real(rt), intent(inout) :: ystate(ystate_lo(1):ystate_hi(1),ystate_lo(2):ystate_hi(2),nc) ! result
      real(rt), intent(inout) :: Imy(Imy_lo(1):Imy_hi(1),Imy_lo(2):Imy_hi(2))
      real(rt), intent(inout) :: Ipy(Ipy_lo(1):Ipy_hi(1),Ipy_lo(2):Ipy_hi(2))
      real(rt), intent(inout) :: sedgey(sedgey_lo(1):sedgey_hi(1),sedgey_lo(2):sedgey_hi(2))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2))
      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2))

      real(rt) :: stxlo(lo(1)-2:hi(1)+2)
      real(rt) :: stxhi(lo(1)-2:hi(1)+2)
      real(rt) :: stylo(lo(2)-2:hi(2)+2)
      real(rt) :: styhi(lo(2)-2:hi(2)+2)
      real(rt) :: hx, hy, dth, dthx, dthy
      real(rt) :: tr,stx,sty,fu,fv,eps,eps_for_bc
      integer  :: i,j,L,imin,jmin,imax,jmax, place_to_break
      logical  :: ltx,lty
      parameter( eps        = 1.0D-6 )
      parameter( eps_for_bc = 1.0D-10 )

      dth  = half*dt
      dthx = half*dt/dx(1)
      dthy = half*dt/dx(2)
      hx = dx(1)
      hy = dx(2)
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)


      do L=1,nc
         ! c
         ! c     compute the slopes
         ! c
         ! c     trace the state to the cell edges
         ! c
         if (ppm_type .gt. 0) then
            call ppm(lo,hi,&
                 s(s_lo(1),s_lo(2),L),s_lo,s_hi,&
                 u,u_lo,u_hi,&
                 v,v_lo,v_hi,&
                 Ipx,Ipx_lo,Ipx_hi,&
                 Imx,Imx_lo,Imx_hi,&
                 Ipy,Ipy_lo,Ipy_hi,&
                 Imy,Imy_lo,Imy_hi,&
                 sm,sm_lo,sm_hi,&
                 sp,sp_lo,sp_hi,&
                 dsvl,dsvl_lo,dsvl_hi,&
                 sedgex,sedgex_lo,sedgex_hi,&
                 sedgey,sedgey_lo,sedgey_hi,&
                 dx, dt, bc, eps_for_bc, ppm_type)
         else
            call slopes(lo,hi,&
                 s(s_lo(1),s_lo(2),L),s_lo,s_hi,&
                 sx,sx_lo,sx_hi,&
                 sy,sy_lo,sy_hi,&
                 bc)
         end if
         !c
         !c     trace the state to the cell edges
         !c
         if (ppm_type .gt. 0) then
            do i = imin, imax+1
               do j = jmin-1,jmax+1
                  xlo(i,j) = Ipx(i-1,j)
                  xhi(i,j) = Imx(i  ,j)
               end do
            end do
         else
            do i = imin, imax+1
               do j = jmin-1,jmax+1
                  xlo(i,j) = s(i-1,j,L) + (half - dthx*u(i-1,j))*sx(i-1,j)
                  xhi(i,j) = s(i  ,j,L) - (half + dthx*u(i  ,j))*sx(i  ,j)
               end do
            end do
         end if

         if(use_minion.eq.1)then
            do i = imin, imax+1
               do j = jmin-1,jmax+1
                  xlo(i,j) = xlo(i,j) + dth*tf(i-1,j,L)
                  xhi(i,j) = xhi(i,j) + dth*tf(i,  j,L)
               end do
            end do
         end if

         call trans_xbc(lo,hi,&
              s(s_lo(1),s_lo(2),L),s_lo,s_hi,&
              xlo,xlo_lo,xlo_hi,&
              xhi,xhi_lo,xhi_hi,&
              uad,uad_lo,uad_hi,&
              n+L-1, bc, eps_for_bc)

         do j = jmin-1,jmax+1
            do i = imin, imax+1
               fu  = merge(0.d0,one,abs(uad(i,j)).lt.eps)
               stx = merge(xlo(i,j),xhi(i,j),uad(i,j) .ge. 0.0d0)
               xlo(i,j) = fu*stx + (one - fu)*half*(xhi(i,j)+xlo(i,j))
            end do
         end do

         if (ppm_type .gt. 0) then
            do j = jmin,jmax+1
               do i = imin-1,imax+1
                  ylo(i,j) = Ipy(i,j-1)
                  yhi(i,j) = Ipy(i,j  )
               end do
            end do
         else
            do j = jmin,jmax+1
               do i = imin-1,imax+1
                  ylo(i,j) = s(i,j-1,L) + (half - dthy*v(i,j-1))*sy(i,j-1)
                  yhi(i,j) = s(i,j  ,L) - (half + dthy*v(i,j  ))*sy(i,j  )
               end do
            end do
         end if

         if(use_minion.eq.1)then
            do i = imin-1, imax+1
               do j = jmin,jmax+1
                  ylo(i,j) = ylo(i,j) + dth*tf(i,j-1,L)
                  yhi(i,j) = yhi(i,j) + dth*tf(i,j  ,L)
               end do
            end do
         end if

         call trans_ybc(lo,hi,&
              s(s_lo(1),s_lo(2),L),s_lo,s_hi,&
              ylo,ylo_lo,ylo_hi,&
              yhi,yhi_lo,yhi_hi,&
              vad,vad_lo,vad_hi,&
              n+L-1, bc, eps_for_bc)

         do j = jmin,jmax+1
            do i = imin-1,imax+1
               fv  = merge(0.d0,one,abs(vad(i,j)).lt.eps)
               sty = merge(ylo(i,j),yhi(i,j),vad(i,j) .ge. 0.0d0)
               ylo(i,j) = fv*sty + (one - fv)*half*(yhi(i,j)+ylo(i,j))
            end do
         end do
         !c
         !c     compute the xedge state
         !c
         if ((velpred.ne.1) .or. (n+L-1.eq.XVEL)) then
            do j = jmin,jmax
               do i = imin-1,imax+1
                  tr = half*&
                       (vad(i,j+1)+vad(i,j))*&
                       (ylo(i,j+1)-ylo(i,j))/hy

                  if (ppm_type .gt. 0) then
                     stxlo(i+1) = Ipx(i,j)&
                          - dth*tr&
                          + dth*tf(i,j,L)
                     stxhi(i  ) = Imx(i,j)&
                          - dth*tr&
                          + dth*tf(i,j,L)
                  else
                     stxlo(i+1) = s(i,j,L) + (half-dthx*u(i,j))*sx(i,j)&
                          - dth*tr&
                          + dth*tf(i,j,L)
                     stxhi(i  ) = s(i,j,L) - (half+dthx*u(i,j))*sx(i,j)&
                          - dth*tr&
                          + dth*tf(i,j,L)
                  end if
               end do

               if (bc(1,1).eq.EXT_DIR .and. velpred.eq.1) then
                  stxhi(imin) = s(imin-1,j,L)
                  stxlo(imin) = s(imin-1,j,L)
               else if (bc(1,1).eq.EXT_DIR .and. uad(imin,j).ge.0.d0) then
                  stxhi(imin) = s(imin-1,j,L)
                  stxlo(imin) = s(imin-1,j,L)
               else if (bc(1,1).eq.EXT_DIR .and. uad(imin,j).lt.0.d0) then
                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP) then
                  if (n+L-1.eq.XVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWXINFLOW
                        !c     prevent backflow
                        stxhi(imin) = MIN(stxhi(imin),0.d0)
#endif
                        stxlo(imin) = stxhi(imin)
                     else
                        if (uad(imin,j).ge.0.d0) then
#ifndef ALLOWXINFLOW
                           !c     prevent backflow
                           stxhi(imin) = MIN(stxhi(imin),0.d0)
#endif
                           stxlo(imin) = stxhi(imin)
                        endif
                     endif
                  else
                     stxlo(imin) = stxhi(imin)
                  end if
               else if (bc(1,1).eq.REFLECT_EVEN) then
                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1).eq.REFLECT_ODD) then
                  stxhi(imin) = 0.d0
                  stxlo(imin) = stxhi(imin)
               end if

               if (bc(1,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stxlo(imax+1) = s(imax+1,j,L)
                  stxhi(imax+1) = s(imax+1,j,L)
               else if (bc(1,2).eq.EXT_DIR .and. uad(imax+1,j).le.0.d0) then
                  stxlo(imax+1) = s(imax+1,j,L)
                  stxhi(imax+1) = s(imax+1,j,L)
               else if (bc(1,2).eq.EXT_DIR .and. uad(imax+1,j).gt.0.d0) then
                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP) then
                  if (n+L-1.eq.XVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWXINFLOW
                        !c     prevent backflow
                        stxlo(imax+1) = MAX(stxlo(imax+1),0.d0)
#endif
                        stxhi(imax+1) = stxlo(imax+1)
                     else
                        if (uad(imax+1,j).le.0.d0) then
#ifndef ALLOWXINFLOW
                           !c     prevent backflow
                           stxlo(imax+1) = MAX(stxlo(imax+1),0.d0)
#endif
                           stxhi(imax+1) = stxlo(imax+1)
                        endif
                     endif
                  else
                     stxhi(imax+1) = stxlo(imax+1)
                  endif
               else if (bc(1,2).eq.REFLECT_EVEN) then
                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2).eq.REFLECT_ODD) then
                  stxlo(imax+1) = 0.d0
                  stxhi(imax+1) = 0.d0
               end if

               if ( velpred .eq. 1 ) then
                  do i = imin, imax+1
                     ltx = stxlo(i) .le. 0.d0  .and.  stxhi(i) .ge. 0.d0
                     ltx = ltx .or. (abs(stxlo(i)+stxhi(i)) .lt. eps)
                     stx = merge(stxlo(i),stxhi(i),(stxlo(i)+stxhi(i)) .ge. 0.0d0)
                     xstate(i,j,L) = merge(0.d0,stx,ltx)
                  end do
               else
                  do i = imin, imax+1
                     xstate(i,j,L) = merge(stxlo(i),stxhi(i),uedge(i,j) .ge. 0.0d0)
                     xstate(i,j,L) = merge(half*(stxlo(i)+stxhi(i)),xstate(i,j,L),&
                          abs(uedge(i,j)).lt.eps)
                  end do
               end if
               place_to_break = 1
            end do
         end if
         !c
         !c     compute the yedge states
         !c
         if ((velpred.ne.1) .or. (n+L-1.eq.YVEL)) then
            do i = imin, imax
               do j = jmin-1,jmax+1
                  tr = half*&
                       (uad(i+1,j)+uad(i,j))*&
                       (xlo(i+1,j)-xlo(i,j))/hx

                  if (ppm_type .gt. 0) then
                     stylo(j+1)= Ipy(i,j)&
                          - dth*tr&
                          + dth*tf(i,j,L)
                     styhi(j  )= Imy(i,j)&
                          - dth*tr&
                          + dth*tf(i,j,L)
                  else
                     stylo(j+1)= s(i,j,L) + (half-dthy*v(i,j))*sy(i,j)&
                          - dth*tr&
                          + dth*tf(i,j,L)
                     styhi(j  )= s(i,j,L) - (half+dthy*v(i,j))*sy(i,j)&
                          - dth*tr&
                          + dth*tf(i,j,L)
                  end if
               end do

               if (bc(2,1).eq.EXT_DIR .and. velpred .eq. 1) then
                  styhi(jmin) = s(i,jmin-1,L)
                  stylo(jmin) = s(i,jmin-1,L)
               else if (bc(2,1).eq.EXT_DIR .and. vad(i,jmin).ge.0.d0) then
                  styhi(jmin) = s(i,jmin-1,L)
                  stylo(jmin) = s(i,jmin-1,L)
               else if (bc(2,1).eq.EXT_DIR .and. vad(i,jmin).lt.0.d0) then
                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP) then
                  if (n+L-1.eq.YVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWYINFLOW
                        !c     prevent backflow
                        styhi(jmin) = MIN(styhi(jmin),0.d0)
#endif
                        stylo(jmin) = styhi(jmin)
                     else
                        if (vad(i,jmin).ge.0.d0) then
#ifndef ALLOWYINFLOW
                           !c     prevent backflow
                           styhi(jmin) = MIN(styhi(jmin),0.d0)
#endif
                           stylo(jmin) = styhi(jmin)
                        endif
                     endif
                  else
                     stylo(jmin) = styhi(jmin)
                  endif
               else if (bc(2,1).eq.REFLECT_EVEN) then
                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1).eq.REFLECT_ODD) then
                  styhi(jmin) = 0.d0
                  stylo(jmin) = 0.d0
               end if

               if (bc(2,2).eq.EXT_DIR .and. velpred .eq. 1) then
                  stylo(jmax+1) = s(i,jmax+1,L)
                  styhi(jmax+1) = s(i,jmax+1,L)
               else if (bc(2,2).eq.EXT_DIR .and. vad(i,jmax+1).le.0.d0) then
                  stylo(jmax+1) = s(i,jmax+1,L)
                  styhi(jmax+1) = s(i,jmax+1,L)
               else if (bc(2,2).eq.EXT_DIR .and. vad(i,jmax+1).gt.0.d0) then
                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP) then
                  if (n+L-1.eq.YVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWYINFLOW
                        !c     prevent backflow
                        stylo(jmax+1) = MAX(stylo(jmax+1),0.d0)
#endif
                        styhi(jmax+1) = stylo(jmax+1)
                     else
                        if (vad(i,jmax+1).le.0.d0) then
#ifndef ALLOWYINFLOW
                           !c     prevent backflow
                           stylo(jmax+1) = MAX(stylo(jmax+1),0.d0)
#endif
                           styhi(jmax+1) = stylo(jmax+1)
                        endif
                     endif
                  else
                     styhi(jmax+1) = stylo(jmax+1)
                  endif
               else if (bc(2,2).eq.REFLECT_EVEN) then
                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2).eq.REFLECT_ODD) then
                  stylo(jmax+1) = 0.d0
                  styhi(jmax+1) = 0.d0
               end if

               if ( velpred .eq. 1 ) then
                  do j = jmin, jmax+1
                     lty = stylo(j) .le. 0.d0  .and.  styhi(j) .ge. 0.d0
                     lty = lty .or. (abs(stylo(j)+styhi(j)) .lt. eps)
                     sty = merge(stylo(j),styhi(j),(stylo(j)+styhi(j)) .ge. 0.0d0)
                     ystate(i,j,L)=merge(0.d0,sty,lty)

! EM_DEBUG
!if ((i < 4) .and.((j > 14).and.(j < 18))) then
!write(*,*) 'DEBUG IN ESTATE_PREMAC',i,j,n,L,ystate(i,j,L)
!endif
                  end do
               else
                  do j=jmin,jmax+1
                     ystate(i,j,L) = merge(stylo(j),styhi(j),vedge(i,j) .ge. 0.0d0)
                     ystate(i,j,L) = merge(half*(stylo(j)+styhi(j)),ystate(i,j,L),&
                          abs(vedge(i,j)).lt.eps)
                  end do
               end if
               place_to_break = 1
            end do
         end if
      enddo
    end subroutine estate_premac

    subroutine estate_fpu(lo,hi,&
         s,s_lo,s_hi,&
         tf,tf_lo,tf_hi,&
         divu,divu_lo,divu_hi,&
         xlo,xlo_lo,xlo_hi,&
         xhi,xhi_lo,xhi_hi,&
         sx,sx_lo,sx_hi,&
         uedge,uedge_lo,uedge_hi,&
         xstate,xstate_lo,xstate_hi,&
         Imx,Imx_lo,Imx_hi,&
         Ipx,Ipx_lo,Ipx_hi,&
         sedgex,sedgex_lo,sedgex_hi,&
         ylo,ylo_lo,ylo_hi,&
         yhi,yhi_lo,yhi_hi,&
         sy,sy_lo,sy_hi,&
         vedge,vedge_lo,vedge_hi,&
         ystate,ystate_lo,ystate_hi,&
         Imy,Imy_lo,Imy_hi,&
         Ipy,Ipy_lo,Ipy_hi,&
         sedgey,sedgey_lo,sedgey_hi,&
         dsvl,dsvl_lo,dsvl_hi,&
         sm,sm_lo,sm_hi,&
         sp,sp_lo,sp_hi,&
         bc, dt, dx, n, nc, use_minion, iconserv, ppm_type)

      implicit none

      integer, intent(in) :: nc, use_minion, iconserv(nc), ppm_type, bc(SDIM,2,nc), n
      real(rt), intent(in) :: dt, dx(SDIM)

      integer, dimension(2), intent(in) :: s_lo,s_hi,tf_lo,tf_hi,divu_lo,divu_hi,&
           xlo_lo,xlo_hi,xhi_lo,xhi_hi,sx_lo,sx_hi,&
           uedge_lo,uedge_hi,xstate_lo,xstate_hi,Imx_lo,Imx_hi,Ipx_lo,Ipx_hi,sedgex_lo,sedgex_hi,&
           ylo_lo,ylo_hi,yhi_lo,yhi_hi,sy_lo,sy_hi,&
           vedge_lo,vedge_hi,ystate_lo,ystate_hi,Imy_lo,Imy_hi,Ipy_lo,Ipy_hi,sedgey_lo,sedgey_hi,&
           dsvl_lo,dsvl_hi,sm_lo,sm_hi,sp_lo,sp_hi,lo,hi

      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),nc)
      real(rt), intent(in) :: tf(tf_lo(1):tf_hi(1),tf_lo(2):tf_hi(2),nc)
      real(rt), intent(in) :: divu(divu_lo(1):divu_hi(1),divu_lo(2):divu_hi(2))
      real(rt), intent(inout) :: xlo(xlo_lo(1):xlo_hi(1),xlo_lo(2):xlo_hi(2))
      real(rt), intent(inout) :: xhi(xhi_lo(1):xhi_hi(1),xhi_lo(2):xhi_hi(2))
      real(rt), intent(inout) :: sx(sx_lo(1):sx_hi(1),sx_lo(2):sx_hi(2))
      real(rt), intent(in) :: uedge(uedge_lo(1):uedge_hi(1),uedge_lo(2):uedge_hi(2))
      real(rt), intent(inout) :: xstate(xstate_lo(1):xstate_hi(1),xstate_lo(2):xstate_hi(2),nc) ! result
      real(rt), intent(inout) :: Imx(Imx_lo(1):Imx_hi(1),Imx_lo(2):Imx_hi(2))
      real(rt), intent(inout) :: Ipx(Ipx_lo(1):Ipx_hi(1),Ipx_lo(2):Ipx_hi(2))
      real(rt), intent(inout) :: sedgex(sedgex_lo(1):sedgex_hi(1),sedgex_lo(2):sedgex_hi(2))
      real(rt), intent(inout) :: ylo(ylo_lo(1):ylo_hi(1),ylo_lo(2):ylo_hi(2))
      real(rt), intent(inout) :: yhi(yhi_lo(1):yhi_hi(1),yhi_lo(2):yhi_hi(2))
      real(rt), intent(inout) :: sy(sy_lo(1):sy_hi(1),sy_lo(2):sy_hi(2))
      real(rt), intent(in) :: vedge(vedge_lo(1):vedge_hi(1),vedge_lo(2):vedge_hi(2))
      real(rt), intent(inout) :: ystate(ystate_lo(1):ystate_hi(1),ystate_lo(2):ystate_hi(2),nc) ! result
      real(rt), intent(inout) :: Imy(Imy_lo(1):Imy_hi(1),Imy_lo(2):Imy_hi(2))
      real(rt), intent(inout) :: Ipy(Ipy_lo(1):Ipy_hi(1),Ipy_lo(2):Ipy_hi(2))
      real(rt), intent(inout) :: sedgey(sedgey_lo(1):sedgey_hi(1),sedgey_lo(2):sedgey_hi(2))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2))
      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2))

      real(rt) :: stxlo(lo(1)-2:hi(1)+2)
      real(rt) :: stxhi(lo(1)-2:hi(1)+2)
      real(rt) :: stylo(lo(2)-2:hi(2)+2)
      real(rt) :: styhi(lo(2)-2:hi(2)+2)
      real(rt) :: hx, hy, dth, dthx, dthy
      real(rt) :: tr,ubar,vbar,stx,sty,fu,fv,eps,eps_for_bc,st
      integer  :: i,j,L,imin,jmin,imax,jmax, inc,place_to_break
      parameter( eps        = 1.0D-6 )
      parameter( eps_for_bc = 1.0D-10 )

      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      hx   = dx(1)
      hy   = dx(2)
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      do L=1,nc
!c
!c     compute the slopes
!c
      if (ppm_type .gt. 0) then
         call ppm_fpu(lo, hi,&
              s(s_lo(1),s_lo(2),L),s_lo,s_hi,&
              uedge,uedge_lo,uedge_hi,&
              vedge,vedge_lo,vedge_hi,&
              Ipx,Ipx_lo,Ipx_hi,&
              Imx,Imx_lo,Imx_hi,&
              Ipy,Ipy_lo,Ipy_hi,&
              Imy,Imy_lo,Imy_hi,&
              sm,sm_lo,sm_hi,&
              sp,sp_lo,sp_hi,&
              dsvl,dsvl_lo,dsvl_hi,&
              sedgex,sedgex_lo,sedgex_hi,&
              sedgey,sedgey_lo,sedgey_hi,&
              dx, dt, bc(1,1,L), eps_for_bc, ppm_type)
      else
         call slopes(lo,hi,&
              s(s_lo(1),s_lo(2),L),s_lo,s_hi,&
              sx,sx_lo,sx_hi,&
              sy,sy_lo,sy_hi,&
              bc(1,1,L))
      end if

!c
!c     trace the state to the cell edges
!c
      if (ppm_type .gt. 0) then
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j) = Ipx(i-1,j)
               xhi(i,j) = Imx(i  ,j)
            end do
         end do
      else
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j) = s(i-1,j,L) + (half - dthx*uedge(i,j))*sx(i-1,j)
               xhi(i,j) = s(i,  j,L) - (half + dthx*uedge(i,j))*sx(i,  j)

! EM_DEBUG
!if ((i < 4) .and.((j > 14).and.(j < 18))) then
!write(*,*) 'DEBUG COMPUTE XLO ',i,j,xlo(i,j),s(i-1,j,L),uedge(i,j),sx(i-1,j)
!endif

            end do
         end do
      end if


      if(use_minion.eq.1)then
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j) = xlo(i,j) + dth*tf(i-1,j,L)
               xhi(i,j) = xhi(i,j) + dth*tf(i,  j,L)
            end do
         end do
         if (iconserv(L) .eq. 1) then
           do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j) = xlo(i,j) - dth*s(i-1,j,L)*divu(i-1,j)
               xhi(i,j) = xhi(i,j) - dth*s(i  ,j,L)*divu(i,  j)
            end do
           end do
         end if
      end if

      call trans_xbc(lo,hi,&
           s(s_lo(1),s_lo(2),L),s_lo,s_hi,&
           xlo,xlo_lo,xlo_hi,&
           xhi,xhi_lo,xhi_hi,&
           uedge,uedge_lo,uedge_hi,&
           n+L-1, bc(1,1,L), eps_for_bc)

      do j = jmin-1,jmax+1
         do i = imin,  imax+1
            fu  = merge(0.d0,one,abs(uedge(i,j)).lt.eps)
            stx = merge(xlo(i,j),xhi(i,j),uedge(i,j) .ge. 0.0d0)
            xlo(i,j) = fu*stx + (one - fu)*half*(xhi(i,j)+xlo(i,j))
         end do
      end do

      if (ppm_type .gt. 0) then
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               ylo(i,j) = Ipy(i,j-1)
               yhi(i,j) = Imy(i,j  )
            end do
         end do
      else
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               ylo(i,j) = s(i,j-1,L) + (half - dthy*vedge(i,j))*sy(i,j-1)
               yhi(i,j) = s(i,j  ,L) - (half + dthy*vedge(i,j))*sy(i,j)

! EM_DEBUG
!if ((i < 4) .and.((j > 14).and.(j < 18))) then
!write(*,*) 'DEBUG COMPUTE YLO ',i,j,ylo(i,j),s(i,j-1,L),vedge(i,j),sy(i,j-1)
!endif

            end do
         end do
      end if

      if (use_minion.eq.1)then
         do j = jmin, jmax+1
            do i = imin-1,  imax+1
               ylo(i,j) = ylo(i,j) + dth*tf(i,j-1,L)
               yhi(i,j) = yhi(i,j) + dth*tf(i,j  ,L)
            end do
         end do
         if (iconserv(L) .eq. 1) then
           do j = jmin-1,jmax+1
            do i = imin,  imax+1
               ylo(i,j) = ylo(i,j) - dth*s(i,j-1,L)*divu(i,j-1)
               yhi(i,j) = yhi(i,j) - dth*s(i,j  ,L)*divu(i,j  )
            end do
           end do
         end if
      end if

      call trans_ybc(lo,hi,&
           s(s_lo(1),s_lo(2),L),s_lo,s_hi,&
           ylo,ylo_lo,ylo_hi,&
           yhi,yhi_lo,yhi_hi,&
           vedge,vedge_lo,vedge_hi,&
           n+L-1, bc(1,1,L), eps_for_bc)

      do j = jmin,  jmax+1
         do i = imin-1,imax+1
            fv  = merge(0.d0,one,abs(vedge(i,j)).lt.eps)
            sty = merge(ylo(i,j),yhi(i,j),vedge(i,j) .ge. 0.0d0)
            ylo(i,j) = fv*sty + (one - fv)*half*(yhi(i,j)+ylo(i,j))
         end do
      end do

!c
!c     compute the xedge states
!c
      do j = jmin,jmax

            do i = imin-1,imax+1

               if (iconserv(L).eq.1) then

                  st = -(vedge(i,j+1)*ylo(i,j+1) - vedge(i,j)*ylo(i,j))/hy&
                      + s(i,j,L)*(vedge(i,j+1)-vedge(i,j))/hy&
                      - s(i,j,L)*divu(i,j)
               else
                  if (vedge(i,j)*vedge(i,j+1).le.0.d0) then
                     vbar = 0.5d0*(vedge(i,j)+vedge(i,j+1))
                     if (vbar.lt.0.d0) then
                        inc = 1
                     else
                        inc = 0
                     endif
                     tr = vbar*(s(i,j+inc,L)-s(i,j+inc-1,L))/hy
                  else
                     tr = half*(vedge(i,j+1) + vedge(i,j)) *&
                         (  ylo(i,j+1) - ylo(i,j)  ) / hy
                  endif
                  st =  -tr
               endif

               if (ppm_type .gt. 0) then
                  stxlo(i+1)= Ipx(i,j)&
                      + dth*(st + tf(i,j,L))
                  stxhi(i  )= Imx(i,j)&
                      + dth*(st + tf(i,j,L))
               else
                  stxlo(i+1)= s(i,j,L) + (half-dthx*uedge(i+1,j))*sx(i,j)&
                      + dth*(st + tf(i,j,L))
                  stxhi(i  )= s(i,j,L) - (half+dthx*uedge(i  ,j))*sx(i,j)&
                      + dth*(st + tf(i,j,L))
               end if
            end do

            if (bc(1,1,L).eq.EXT_DIR .and. uedge(imin,j).ge.0.d0) then
               stxhi(imin) = s(imin-1,j,L)
               stxlo(imin) = s(imin-1,j,L)
            else if (bc(1,1,L).eq.EXT_DIR .and. uedge(imin,j).lt.0.d0) then
               stxlo(imin) = stxhi(imin)
            else if (bc(1,1,L).eq.FOEXTRAP.or.bc(1,1,L).eq.HOEXTRAP) then
               if (n+L-1.eq.XVEL) then
                  if (uedge(imin,j).ge.0.d0) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                     stxhi(imin) = MIN(stxhi(imin),0.d0)
#endif
                     stxlo(imin) = stxhi(imin)
                  endif
               else
                  stxlo(imin) = stxhi(imin)
               endif
            else if (bc(1,1,L).eq.REFLECT_EVEN) then
               stxlo(imin) = stxhi(imin)
            else if (bc(1,1,L).eq.REFLECT_ODD) then
               stxhi(imin) = 0.d0
               stxlo(imin) = 0.d0
            end if
            if (bc(1,2,L).eq.EXT_DIR .and. uedge(imax+1,j).le.0.d0) then
               stxlo(imax+1) = s(imax+1,j,L)
               stxhi(imax+1) = s(imax+1,j,L)
            else if (bc(1,2,L).eq.EXT_DIR .and. uedge(imax+1,j).gt.0.d0) then
               stxhi(imax+1) = stxlo(imax+1)
            else if (bc(1,2,L).eq.FOEXTRAP.or.bc(1,2,L).eq.HOEXTRAP) then
               if (n+L-1.eq.XVEL) then
                  if (uedge(imax+1,j).le.0.d0) then
#ifndef ALLOWXINFLOW
!c     prevent backflow
                     stxlo(imax+1) = MAX(stxlo(imax+1),0.d0)
#endif
                     stxhi(imax+1) = stxlo(imax+1)
                  endif
               else
                  stxhi(imax+1) = stxlo(imax+1)
               endif
            else if (bc(1,2,L).eq.REFLECT_EVEN) then
               stxhi(imax+1) = stxlo(imax+1)
            else if (bc(1,2,L).eq.REFLECT_ODD) then
               stxlo(imax+1) = 0.d0
               stxhi(imax+1) = 0.d0
            end if

            do i = imin, imax+1
! EM_DEBUG
!if ((i < 4) .and.((j > 14).and.(j < 18))) then
!write(*,*) "DEBUG IN ESTATE_FPU ",i,j,L,stxlo(i),stxhi(i),uedge(i,j)
!endif
               xstate(i,j,L) = merge(stxlo(i),stxhi(i),uedge(i,j) .ge. 0.0d0)
               xstate(i,j,L) = merge(half*(stxlo(i)+stxhi(i)),xstate(i,j,L)&
                   ,abs(uedge(i,j)).lt.eps)
            end do
            place_to_break = 1
      end do

!c
!c     compute the yedge states
!c
      do i = imin,imax

            do j = jmin-1,jmax+1

               if (iconserv(L).eq.1) then

                  st = -(uedge(i+1,j)*xlo(i+1,j) - uedge(i,j)*xlo(i,j))/hx&
                      + s(i,j,L)*(uedge(i+1,j)-uedge(i,j))/hx&
                      - s(i,j,L)*divu(i,j)

               else

                  if (uedge(i,j)*uedge(i+1,j).le.0.d0) then
                     ubar = 0.5d0*(uedge(i,j)+uedge(i+1,j))
                     if (ubar.lt.0.d0) then
                        inc = 1
                     else
                        inc = 0
                     endif
                     tr = ubar*(s(i+inc,j,L)-s(i+inc-1,j,L))/hx
                  else
                     tr = half*(uedge(i+1,j) + uedge(i,j)) *&
                                (xlo(i+1,j) -   xlo(i,j)  ) / hx
                  endif
                  st = -tr

               endif

               if (ppm_type .gt. 0) then
                  stylo(j+1)= Ipy(i,j)&
                      + dth*(st + tf(i,j,L))
                  styhi(j  )= Imy(i,j)&
                      + dth*(st + tf(i,j,L))
               else
                  stylo(j+1)= s(i,j,L) + (half-dthy*vedge(i,j+1))*sy(i,j)&
                      + dth*(st + tf(i,j,L))
                  styhi(j  )= s(i,j,L) - (half+dthy*vedge(i,j  ))*sy(i,j)&
                      + dth*(st + tf(i,j,L))
               end if

            end do

            if (bc(2,1,L).eq.EXT_DIR .and. vedge(i,jmin).ge.0.d0) then
               styhi(jmin) = s(i,jmin-1,L)
               stylo(jmin) = s(i,jmin-1,L)
            else if (bc(2,1,L).eq.EXT_DIR .and. vedge(i,jmin).lt.0.d0) then
               stylo(jmin) = styhi(jmin)
            else if (bc(2,1,L).eq.FOEXTRAP.or.bc(2,1,L).eq.HOEXTRAP) then
               if (L+n-1.eq.YVEL) then
                  if (vedge(i,jmin).ge.0.d0) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                     styhi(jmin) = MIN(styhi(jmin),0.d0)
#endif
                     stylo(jmin) = styhi(jmin)
                  endif
               else
                  stylo(jmin) = styhi(jmin)
               endif
            else if (bc(2,1,L).eq.REFLECT_EVEN) then
               stylo(jmin) = styhi(jmin)
            else if (bc(2,1,L).eq.REFLECT_ODD) then
               styhi(jmin) = 0.d0
               stylo(jmin) = 0.d0
            end if

            if (bc(2,2,L).eq.EXT_DIR .and. vedge(i,jmax+1).le.0.d0) then
               stylo(jmax+1) = s(i,jmax+1,L)
               styhi(jmax+1) = s(i,jmax+1,L)
            else if (bc(2,2,L).eq.EXT_DIR .and. vedge(i,jmax+1).gt.0.d0) then
               styhi(jmax+1) = stylo(jmax+1)
            else if (bc(2,2,L).eq.FOEXTRAP.or.bc(2,2,L).eq.HOEXTRAP) then
               if (n+L-1.eq.YVEL) then
                  if (vedge(i,jmax+1).le.0.d0) then
#ifndef ALLOWYINFLOW
!c     prevent backflow
                     stylo(jmax+1) = MAX(stylo(jmax+1),0.d0)
#endif
                     styhi(jmax+1) = stylo(jmax+1)
                  endif
               else
                  styhi(jmax+1) = stylo(jmax+1)
               endif
            else if (bc(2,2,L).eq.REFLECT_EVEN) then
               styhi(jmax+1) = stylo(jmax+1)
            else if (bc(2,2,L).eq.REFLECT_ODD) then
               stylo(jmax+1) = 0.d0
               styhi(jmax+1) = 0.d0
            end if

            do j=jmin,jmax+1
               ystate(i,j,L) = merge(stylo(j),styhi(j),vedge(i,j) .ge. 0.0d0)
               ystate(i,j,L) = merge(half*(stylo(j)+styhi(j)),ystate(i,j,L),&
                   abs(vedge(i,j)).lt.eps)
            end do
            place_to_break = 1
         end do
      end do

    end subroutine estate_fpu

    subroutine trans_xbc(lo,hi,&
         s,s_lo,s_hi,&
         xlo,xlo_lo,xlo_hi,&
         xhi,xhi_lo,xhi_hi,&
         uad,uad_lo,uad_hi,&
         n, xbc, eps)
! c
! c     This subroutine processes boundary conditions on information
! c     traced to cell faces in the x direction.  This is used for
! c     computing velocities and edge states used in calculating
! c     transverse derivatives
! c
      implicit none

      integer, intent(in) :: n,xbc(SDIM,2)
      integer, dimension(2), intent(in) :: &
           s_lo,s_hi,xlo_lo,xlo_hi,xhi_lo,xhi_hi,uad_lo,uad_hi,lo,hi

      real(rt), intent(in)    :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2))
      real(rt), intent(inout) :: xlo(xlo_lo(1):xlo_hi(1),xlo_lo(2):xlo_hi(2))
      real(rt), intent(inout) :: xhi(xhi_lo(1):xhi_hi(1),xhi_lo(2):xhi_hi(2))
      real(rt), intent(in)    :: uad(uad_lo(1):uad_hi(1),uad_lo(2):uad_hi(2))
      real(rt), intent(in)    ::  eps

      real(rt) ::  stx
      logical ltest
      integer j
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      if (xbc(1,1).eq.EXT_DIR) then
         if (n .eq. XVEL) then
            do j = jmin-1,jmax+1
              if (uad(imin,j) .ge. 0.d0) then
                 xhi(imin,j) = s(imin-1,j)
                 xlo(imin,j) = s(imin-1,j)
              else
                 xlo(imin,j) = xhi(imin,j)
              end if
            end do
         else
            do j = jmin-1,jmax+1
               ltest = uad(imin,j).le.eps
               stx   = merge(xhi(imin,j),s(imin-1,j),ltest)
               xhi(imin,j) = stx
               xlo(imin,j) = stx
            end do
         end if
      else if (xbc(1,1).eq.FOEXTRAP.or.xbc(1,1).eq.HOEXTRAP&
             .or.xbc(1,1).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            xlo(imin,j) = xhi(imin,j)
         end do
      else if (xbc(1,1).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            xhi(imin,j) = 0.d0
            xlo(imin,j) = 0.d0
         end do
      end if

      if (xbc(1,2).eq.EXT_DIR) then
         if (n .eq. XVEL) then
            do j = jmin-1,jmax+1
              if (uad(imax+1,j) .le. 0.d0) then
                 xhi(imax+1,j) = s(imax+1,j)
                 xlo(imax+1,j) = s(imax+1,j)
               else
                 xhi(imax+1,j) = xlo(imax+1,j)
               end if
             end do
         else
            do j = jmin-1,jmax+1
               ltest = uad(imax+1,j).ge.-eps
               stx   = merge(xlo(imax+1,j),s(imax+1,j),ltest)
               xhi(imax+1,j) = stx
               xlo(imax+1,j) = stx
            end do
         end if
      else if (xbc(1,2).eq.FOEXTRAP.or.xbc(1,2).eq.HOEXTRAP&
             .or.xbc(1,2).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            xhi(imax+1,j) = xlo(imax+1,j)
         end do
      else if (xbc(1,2).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            xhi(imax+1,j) = 0.d0
            xlo(imax+1,j) = 0.d0
         end do
      end if

    end subroutine trans_xbc

    subroutine trans_ybc(lo,hi,&
         s,s_lo,s_hi,&
         ylo,ylo_lo,ylo_hi,&
         yhi,yhi_lo,yhi_hi,&
         vad,vad_lo,vad_hi,&
         n, ybc, eps)
! c
! c     This subroutine processes boundary conditions on information
! c     traced to cell faces in the y direction.  This is used for
! c     computing velocities and edge states used in calculating
! c     transverse derivatives
! c
      implicit none

      integer, intent(in) :: n,ybc(SDIM,2)
      integer, dimension(2), intent(in) :: &
           s_lo,s_hi,ylo_lo,ylo_hi,yhi_lo,yhi_hi,vad_lo,vad_hi,lo,hi

      real(rt), intent(in)    :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2))
      real(rt), intent(inout) :: ylo(ylo_lo(1):ylo_hi(1),ylo_lo(2):ylo_hi(2))
      real(rt), intent(inout) :: yhi(yhi_lo(1):yhi_hi(1),yhi_lo(2):yhi_hi(2))
      real(rt), intent(in)    :: vad(vad_lo(1):vad_hi(1),vad_lo(2):vad_hi(2))
      real(rt), intent(in)    ::  eps

      real(rt) :: sty
      logical ltest
      integer i
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      if (ybc(2,1).eq.EXT_DIR) then
         if (n .eq. YVEL) then
            do i = imin-1,imax+1
              if (vad(i,jmin).ge.0.d0) then
                 yhi(i,jmin) = s(i,jmin-1)
                 ylo(i,jmin) = s(i,jmin-1)
              else
                 ylo(i,jmin) = yhi(i,jmin)
              end if
            end do
         else
            do i = imin-1,imax+1
               ltest = vad(i,jmin).le.eps
               sty   = merge(yhi(i,jmin),s(i,jmin-1),ltest)
               yhi(i,jmin) = sty
               ylo(i,jmin) = sty
            end do
         end if
      else if (ybc(2,1).eq.FOEXTRAP.or.ybc(2,1).eq.HOEXTRAP&
             .or.ybc(2,1).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            ylo(i,jmin) = yhi(i,jmin)
         end do
      else if (ybc(2,1).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            yhi(i,jmin) = 0.d0
            ylo(i,jmin) = 0.d0
         end do
      end if

      if (ybc(2,2).eq.EXT_DIR) then
         if (n .eq. YVEL) then
            do i = imin-1,imax+1
              if (vad(i,jmax+1).le.0.d0) then
                 ylo(i,jmax+1) = s(i,jmax+1)
                 yhi(i,jmax+1) = s(i,jmax+1)
              else
                 yhi(i,jmax+1) = ylo(i,jmax+1)
              end if
            end do
         else
            do i = imin-1,imax+1
               ltest = vad(i,jmax+1).ge.-eps
               sty   = merge(ylo(i,jmax+1),s(i,jmax+1),ltest)
               yhi(i,jmax+1) = sty
               ylo(i,jmax+1) = sty
            end do
         end if
      else if (ybc(2,2).eq.FOEXTRAP.or.ybc(2,2).eq.HOEXTRAP&
             .or.ybc(2,2).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            yhi(i,jmax+1) = ylo(i,jmax+1)
         end do
      else if (ybc(2,2).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            ylo(i,jmax+1) = 0.d0
            yhi(i,jmax+1) = 0.d0
         end do
      end if

    end subroutine trans_ybc

    subroutine slopes(lo,hi,&
         s,s_lo,s_hi,&
         slx,slx_lo,slx_hi,&
         sly,sly_lo,sly_hi,&
         bc)
! c
! c     this subroutine computes first, second or forth order slopes of
! c     a 2D scalar field.
! c
! c     Boundary conditions on interior slopes are handled automatically
! c     by the ghost cells
! c
! c     Boundary conditions on EXT_DIR and HOEXTRAP slopes are implemented
! c     by setting them to 0.d0 outside of the domain and using a
! c     one-sided derivative from the interior
! c
      implicit none

#include <GODCOMM_F.H>

      integer, intent(in) :: bc(SDIM,2)
      integer, dimension(2), intent(in) :: lo,hi,s_lo,s_hi,slx_lo,slx_hi,sly_lo,sly_hi
      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2))
      real(rt), intent(inout) :: slx(slx_lo(1):slx_hi(1),slx_lo(2):slx_hi(2))
      real(rt), intent(inout) :: sly(sly_lo(1):sly_hi(1),sly_lo(2):sly_hi(2))
      real(rt) :: slxscr(lo(1)-2:hi(1)+2,4)
      real(rt) :: slyscr(lo(2)-2:hi(2)+2,4)

      integer cen,lim,flag,fromm
      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flag = 3 )
      parameter( fromm = 4 )

      integer imin,jmin,imax,jmax,i,j
      integer ng
      real(rt) dpls,dmin,ds
      real(rt) del,slim,sflg

! C
! C     Determine ng in a way that covers the case of tiling where
! C     (lo:hi) is only a portion of the box s is defined on.
! C
      ng = lo(1) - s_lo(1)
      if (slope_order .eq.1) then
         if (ng .lt. 1) then
            call bl_abort('slopes: too few bndry cells for first order')
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
      imax = hi(1)
      jmax = hi(2)
!c
!c ::: ::::: added to prevent underflow for small s values
!c
!      do j = lo(2)-ng, hi(2)+ng
!        do i = lo(1)-ng, hi(1)+ng
!           s(i,j) = merge(s(i,j), 0.d0, abs(s(i,j)).gt.1.0D-20)
!       end do
!      end do
!c
!c     COMPUTE 0TH order slopes
!c
      if (slope_order.eq.1) then
        do j = jmin-1,jmax+1
           do i = imin-1,imax+1
              slx(i,j) = 0.d0
              sly(i,j) = 0.d0
  	   end do
        end do
      end if
!c
!c     COMPUTE 2nd order slopes
!c
      if (slope_order.eq.2) then
!c
!c     ------------------------ x slopes
!c
         if (use_unlimited_slopes) then
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  slx(i,j) = half*(s(i+1,j)-s(i-1,j))
               end do
            end do
            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imin-1,j) = 0.d0
                  slx(imin,j)   = (s(imin+1,j)+three*s(imin,j)-four*s(imin-1,j))/three
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imax+1,j) = 0.d0
                  slx(imax,j)   = -(s(imax-1,j)+three*s(imax,j)-four*s(imax+1,j))/three
               end do
            end if
         else
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  del  = half*(s(i+1,j)-s(i-1,j))
                  dpls = two*(s(i+1,j) - s(i ,j))
                  dmin = two*(s(i ,j) - s(i-1,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, 0.d0, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  slx(i,j)= sflg*min(slim,abs(del))
               end do
            end do

            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imin-1,j) = 0.d0
                  del  = (s(imin+1,j)+three*s(imin,j)-four*s(imin-1,j))/three
                  dpls = two*(s(imin+1,j) - s(imin  ,j))
                  dmin = two*(s(imin  ,j) - s(imin-1,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, 0.d0, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  slx(imin,j)= sflg*min(slim,abs(del))
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imax+1,j) = 0.d0
                  del  = -(s(imax-1,j)+three*s(imax,j)-four*s(imax+1,j))/three
                  dpls = two*(s(imax+1,j) - s(imax  ,j))
                  dmin = two*(s(imax  ,j) - s(imax-1,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, 0.d0, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  slx(imax,j)= sflg*min(slim,abs(del))
               end do
            end if
         end if
!c
!c     ------------------------ y slopes
!c
         if (use_unlimited_slopes) then
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  sly(i,j) = half*(s(i,j+1)-s(i,j-1))
               end do
            end do
            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmin-1) = 0.d0
                  sly(i,jmin) = (s(i,jmin+1)+three*s(i,jmin)-four*s(i,jmin-1))/three
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmax+1) = 0.d0
                  sly(i,jmax) = -(s(i,jmax-1)+three*s(i,jmax)-four*s(i,jmax+1))/three
               end do
            end if
         else
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  del  = half*(s(i,j+1)-s(i,j-1))
                  dpls = two*(s(i,j+1) - s(i,j ))
                  dmin = two*(s(i,j ) - s(i,j-1))
                  slim = min(abs(dpls),abs(dmin))
                  slim = merge(slim, 0.d0, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  sly(i,j)= sflg*min(slim,abs(del))
               end do
            end do

            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmin-1) = 0.d0
                  del  = (s(i,jmin+1)+three*s(i,jmin)-four*s(i,jmin-1))/three
                  dpls = two*(s(i,jmin+1) - s(i,jmin  ))
                  dmin = two*(s(i,jmin  ) - s(i,jmin-1))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, 0.d0, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  sly(i,jmin)= sflg*min(slim,abs(del))
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmax+1) = 0.d0
                  del  = -(s(i,jmax-1)+three*s(i,jmax)-four*s(i,jmax+1))/three
                  dpls = two*(s(i,jmax+1) - s(i,jmax  ))
                  dmin = two*(s(i,jmax  ) - s(i,jmax-1))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, 0.d0, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  sly(i,jmax)= sflg*min(slim,abs(del))
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
      if (slope_order.eq.4) then
!c
!c     ------------------------ x slopes
!c
         if (use_unlimited_slopes) then
            do j = jmin-1,jmax+1
               do i = imin-2,imax+2
                  slxscr(i,cen)  = half*(s(i+1,j)-s(i-1,j))
               end do
               do i = imin-1,imax+1
                  slx(i,j) = two * two3rd * slxscr(i,cen) -&
                      sixth * (slxscr(i+1,cen) + slxscr(i-1,cen))
               end do
            end do

            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imin,j) = -sixteen/fifteen*s(imin-1,j) + half*s(imin,j) + &
                      two3rd*s(imin+1,j) - tenth*s(imin+2,j)
                  slx(imin-1,j) = 0.d0
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imax,j) = -( -sixteen/fifteen*s(imax+1,j) + half*s(imax,j) + &
                      two3rd*s(imax-1,j) - tenth*s(imax-2,j) )
                  slx(imax+1,j) = 0.d0
               end do
            end if
         else
            do j = jmin-1,jmax+1
               do i = imin-2,imax+2
                  dmin           =  two*(s(i,  j)-s(i-1,j))
                  dpls           =  two*(s(i+1,j)-s(i  ,j))
                  slxscr(i,cen)  = half*(s(i+1,j)-s(i-1,j))
                  slxscr(i,lim)  = min(abs(dmin),abs(dpls))
                  slxscr(i,lim)  = merge(slxscr(i,lim),0.d0,(dpls*dmin) .ge. 0.0d0)
                  slxscr(i,flag) = sign(one,slxscr(i,cen))
                  slxscr(i,fromm)= slxscr(i,flag)*&
                      min(slxscr(i,lim),abs(slxscr(i,cen)))
               end do

               do i = imin-1,imax+1
                  ds = two * two3rd * slxscr(i,cen) -&
                      sixth * (slxscr(i+1,fromm) + slxscr(i-1,fromm))
                  slx(i,j) = slxscr(i,flag)*min(abs(ds),slxscr(i,lim))
               end do

               if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                  del  = -sixteen/fifteen*s(imin-1,j) + half*s(imin,j) + &
                      two3rd*s(imin+1,j) - tenth*s(imin+2,j)
                  dmin = two*(s(imin  ,j)-s(imin-1,j))
                  dpls = two*(s(imin+1,j)-s(imin  ,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, 0.d0, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  slx(imin-1,j) = 0.d0
                  slx(imin,  j) = sflg*min(slim,abs(del))

!c                 Recalculate the slope at imin+1 using the revised slxscr(imin,fromm)
                  slxscr(imin,fromm) = slx(imin,j)
                  ds = two * two3rd * slxscr(imin+1,cen) -&
                    sixth * (slxscr(imin+2,fromm) + slxscr(imin,fromm))
                  slx(imin+1,j) = slxscr(imin+1,flag)*min(abs(ds),slxscr(imin+1,lim))
               end if

               if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                  del  = -( -sixteen/fifteen*s(imax+1,j) + half*s(imax,j) + &
                      two3rd*s(imax-1,j) - tenth*s(imax-2,j) )
                  dmin = two*(s(imax  ,j)-s(imax-1,j))
                  dpls = two*(s(imax+1,j)-s(imax  ,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, 0.d0, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  slx(imax,  j) = sflg*min(slim,abs(del))
                  slx(imax+1,j) = 0.d0

!c                 Recalculate the slope at imax-1 using the revised slxscr(imax,fromm)
                  slxscr(imax,fromm) = slx(imax,j)
                  ds = two * two3rd * slxscr(imax-1,cen) -&
                    sixth * (slxscr(imax-2,fromm) + slxscr(imax,fromm))
                  slx(imax-1,j) = slxscr(imax-1,flag)*min(abs(ds),slxscr(imax-1,lim))
               end if
            end do
         end if
!c
!c     ------------------------ y slopes
!c
         if (use_unlimited_slopes) then
            do i = imin-1,imax+1
               do j = jmin-2,jmax+2
                  slyscr(j,cen)  = half*(s(i,j+1)-s(i,j-1))
               end do
               do j = jmin-1,jmax+1
                  sly(i,j) = two * two3rd * slyscr(j,cen) -&
                      sixth * (slyscr(j+1,cen) + slyscr(j-1,cen))
               end do
            end do

            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmin-1) = 0.d0
                  sly(i,jmin) = -sixteen/fifteen*s(i,jmin-1) + half*s(i,jmin) + &
                      two3rd*s(i,jmin+1) - tenth*s(i,jmin+2)
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmax) = -( -sixteen/fifteen*s(i,jmax+1) + half*s(i,jmax) + &
                      two3rd*s(i,jmax-1) - tenth*s(i,jmax-2) )
                  sly(i,jmax+1) = 0.d0
               end do
            end if
         else
            do i = imin-1,imax+1
              do j = jmin-2,jmax+2
                  dmin           =  two*(s(i,j  )-s(i,j-1))
                  dpls           =  two*(s(i,j+1)-s(i,j  ))
                  slyscr(j,cen)  = half*(s(i,j+1)-s(i,j-1))
                  slyscr(j,lim)  = min(abs(dmin),abs(dpls))
                  slyscr(j,lim)  = merge(slyscr(j,lim),0.d0,(dpls*dmin) .ge. 0.0d0)
                  slyscr(j,flag) = sign(one,slyscr(j,cen))
                  slyscr(j,fromm)= slyscr(j,flag)*&
                      min(slyscr(j,lim),abs(slyscr(j,cen)))
               end do
               do j = jmin-1,jmax+1
                  ds = two * two3rd * slyscr(j,cen) -&
                      sixth * (slyscr(j+1,fromm) + slyscr(j-1,fromm))
                  sly(i,j) = slyscr(j,flag)*min(abs(ds),slyscr(j,lim))
               end do

               if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                  del  = -sixteen/fifteen*s(i,jmin-1) + half*s(i,jmin) + &
                      two3rd*s(i,jmin+1) - tenth*s(i,jmin+2)
                  dmin = two*(s(i,jmin  )-s(i,jmin-1))
                  dpls = two*(s(i,jmin+1)-s(i,jmin  ))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, 0.d0, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  sly(i,jmin-1) = 0.d0
                  sly(i,jmin  ) = sflg*min(slim,abs(del))

!c                 Recalculate the slope at jmin+1 using the revised slyscr(jmin,fromm)
                  slyscr(jmin,fromm) = sly(i,jmin)
                  ds = two * two3rd * slyscr(jmin+1,cen) -&
                    sixth * (slyscr(jmin+2,fromm) + slyscr(jmin,fromm))
                  sly(i,jmin+1) = slyscr(jmin+1,flag)*min(abs(ds),slyscr(jmin+1,lim))
               end if

               if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                  del  = -( -sixteen/fifteen*s(i,jmax+1) + half*s(i,jmax) + &
                      two3rd*s(i,jmax-1) - tenth*s(i,jmax-2) )
                  dmin = two*(s(i,jmax  )-s(i,jmax-1))
                  dpls = two*(s(i,jmax+1)-s(i,jmax  ))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, 0.d0, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  sly(i,jmax  ) = sflg*min(slim,abs(del))
                  sly(i,jmax+1) = 0.d0

!c                 Recalculate the slope at jmax-1 using the revised slyscr(jmax,fromm)
                  slyscr(jmax,fromm) = sly(i,jmax)
                  ds = two * two3rd * slyscr(jmax-1,cen) -&
                    sixth * (slyscr(jmax-2,fromm) + slyscr(jmax,fromm))
                  sly(i,jmax-1) = slyscr(jmax-1,flag)*min(abs(ds),slyscr(jmax-1,lim))
               end if
           end do
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
         Ipx,Ipx_lo,Ipx_hi,&
         Imx,Imx_lo,Imx_hi,&
         Ipy,Ipy_lo,Ipy_hi,&
         Imy,Imy_lo,Imy_hi,&
         sm,sm_lo,sm_hi,&
         sp,sp_lo,sp_hi,&
         dsvl,dsvl_lo,dsvl_hi,&
         sedgex,sedgex_lo,sedgex_hi,&
         sedgey,sedgey_lo,sedgey_hi,&
         dx, dt, bc, eps, ppm_type)

      implicit none

      integer, intent(in) :: bc(SDIM,2)
      integer, dimension(2), intent(in) :: &
           s_lo,s_hi,u_lo,u_hi,v_lo,v_hi,Ipx_lo,Ipx_hi,Imx_lo,Imx_hi,Ipy_lo,Ipy_hi,Imy_lo,Imy_hi,&
           sm_lo,sm_hi,sp_lo,sp_hi,dsvl_lo,dsvl_hi,sedgex_lo,sedgex_hi,sedgey_lo,sedgey_hi,&
           lo,hi

      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2))
      real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2))
      real(rt), intent(in) :: v(v_lo(1):v_hi(1),v_lo(2):v_hi(2))
      real(rt), intent(inout) :: Ipx(Ipx_lo(1):Ipx_hi(1),Ipx_lo(2):Ipx_hi(2))
      real(rt), intent(inout) :: Imx(Imx_lo(1):Imx_hi(1),Imx_lo(2):Imx_hi(2))
      real(rt), intent(inout) :: Ipy(Ipy_lo(1):Ipy_hi(1),Ipy_lo(2):Ipy_hi(2))
      real(rt), intent(inout) :: Imy(Imy_lo(1):Imy_hi(1),Imy_lo(2):Imy_hi(2))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2))
      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2))
      real(rt), intent(inout) :: sedgex(sedgex_lo(1):sedgex_hi(1),sedgex_lo(2):sedgex_hi(2))
      real(rt), intent(inout) :: sedgey(sedgey_lo(1):sedgey_hi(1),sedgey_lo(2):sedgey_hi(2))
      real(rt), intent(in) :: eps, dx(SDIM), dt
      integer ppm_type

!c     local
      integer i,j

      logical extremum, bigp, bigm

      real(rt) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham
      real(rt) :: sgn, sigma, s6, amax, delam, delap
      real(rt) :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp

!c     constant used in Colella 2008
      C = 1.25d0

!c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c     x-direction
!c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!c     compute s at x-edges
      if (ppm_type .eq. 1) then

!c     compute van Leer slopes in x-direction
         dsvl = 0.d0
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-2,hi(1)+2
               dsc = 0.5d0 * (s(i+1,j) - s(i-1,j))
               dsl = 2.d0  * (s(i  ,j) - s(i-1,j))
               dsr = 2.d0  * (s(i+1,j) - s(i  ,j))
               if (dsl*dsr .gt. 0.d0)&
                   dsvl(i,j) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
            end do
         end do

!c     interpolate s to x-edges
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+2
               sedgex(i,j) = 0.5d0*(s(i,j)+s(i-1,j)) - (1.d0/6.d0)*(dsvl(i,j)-dsvl(i-1,j))
!c     make sure sedgex lies in between adjacent cell-centered values
               sedgex(i,j) = max(sedgex(i,j),min(s(i,j),s(i-1,j)))
               sedgex(i,j) = min(sedgex(i,j),max(s(i,j),s(i-1,j)))
            end do
         end do

!c     copy sedgex into sp and sm
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               sp(i,j) = sedgex(i+1,j)
               sm(i,j) = sedgex(i  ,j)
            end do
         end do

!c     modify using quadratic limiters
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end do

!c     different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
         if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
!c     the value in the first cc ghost cell represents the edge value
            sm(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

!c     use a modified stencil to get sedgex on the first interior edge
            sedgex(lo(1)+1,lo(2)-1:hi(2)+1) = &
                -(1.d0/5.d0)  *s(lo(1)-1,lo(2)-1:hi(2)+1)&
                +(3.d0/4.d0)  *s(lo(1)  ,lo(2)-1:hi(2)+1)&
                +0.5d0        *s(lo(1)+1,lo(2)-1:hi(2)+1)&
                -(1.d0/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

!c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(lo(1)+1,j) = max(sedgex(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
               sedgex(lo(1)+1,j) = min(sedgex(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
            end do

!c     copy sedgex into sp and sm
            do j=lo(2)-1,hi(2)+1
               sp(lo(1)  ,j) = sedgex(lo(1)+1,j)
               sm(lo(1)+1,j) = sedgex(lo(1)+1,j)
            end do

!c     reset sp on second interior edge
            do j=lo(2)-1,hi(2)+1
               sp(lo(1)+1,j) = sedgex(lo(1)+2,j)
            end do

!c     modify using quadratic limiters
            do j=lo(2)-1,hi(2)+1
               i = lo(1)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

         if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
!c     the value in the first cc ghost cell represents the edge value
            sp(hi(1),lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

!c     use a modified stencil to get sedgex on the first interior edge
            sedgex(hi(1),lo(2)-1:hi(2)+1) =&
                -(1.d0/5.d0)  *s(hi(1)+1,lo(2)-1:hi(2)+1)&
                +(3.d0/4.d0)  *s(hi(1)  ,lo(2)-1:hi(2)+1)&
                +0.5d0        *s(hi(1)-1,lo(2)-1:hi(2)+1)&
                -(1.d0/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

!c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(hi(1),j) = max(sedgex(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
               sedgex(hi(1),j) = min(sedgex(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
            end do

!c     copy sedgex into sp and sm
            do j=lo(2)-1,hi(2)+1
               sp(hi(1)-1,j) = sedgex(hi(1),j)
               sm(hi(1)  ,j) = sedgex(hi(1),j)
            end do

!c     reset sm on second interior edge
            do j=lo(2)-1,hi(2)+1
               sm(hi(1)-1,j) = sedgex(hi(1)-1,j)
            end do

!c     modify using quadratic limiters
            do j=lo(2)-1,hi(2)+1
               i = hi(1)-1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

      else if (ppm_type .eq. 2) then

!c     interpolate s to x-edges
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-2,hi(1)+3
               sedgex(i,j) = (7.d0/12.d0)*(s(i-1,j)+s(i,j)) - (1.d0/12.d0)*(s(i-2,j)+s(i+1,j))
!c     limit sedgex
               if ((sedgex(i,j)-s(i-1,j))*(s(i,j)-sedgex(i,j)) .lt. 0.d0) then
                  D2  = 3.d0*(s(i-1,j)-2.d0*sedgex(i,j)+s(i,j))
                  D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                  D2R = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                  sgn = sign(1.d0,D2)
                  D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                  sedgex(i,j) = 0.5d0*(s(i-1,j)+s(i,j)) - (1.d0/6.d0)*D2LIM
               end if
            end do
         end do

!c     use Colella 2008 limiters
!c     This is a new version of the algorithm
!c     to eliminate sensitivity to roundoff.
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1

               alphap = sedgex(i+1,j)-s(i,j)
               alpham = sedgex(i  ,j)-s(i,j)
               bigp = abs(alphap).gt.2.d0*abs(alpham)
               bigm = abs(alpham).gt.2.d0*abs(alphap)
               extremum = .false.

               if (alpham*alphap .ge. 0.d0) then
                  extremum = .true.
               else if (bigp .or. bigm) then
!c     Possible extremum. We look at cell centered values and face
!c     centered values for a change in sign in the differences adjacent to
!c     the cell. We use the pair of differences whose minimum magnitude is the
!c     largest, and thus least susceptible to sensitivity to roundoff.
                  dafacem = sedgex(i,j) - sedgex(i-1,j)
                  dafacep = sedgex(i+2,j) - sedgex(i+1,j)
                  dabarm = s(i,j) - s(i-1,j)
                  dabarp = s(i+1,j) - s(i,j)
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
                  D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                  D2R = s(i,j)-2.d0*s(i+1,j)+s(i+2,j)
                  D2C = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                  sgn = sign(1.d0,D2)
                  D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                  alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                  alphap = alphap*D2LIM/max(abs(D2),1.d-10)
               else
                  if (bigp) then
                     sgn = sign(1.d0,alpham)
                     amax = -alphap**2 / (4*(alpham + alphap))
                     delam = s(i-1,j) - s(i,j)
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
                     delap = s(i+1,j) - s(i,j)
                     if (sgn*amax .ge. sgn*delap) then
                        if (sgn*(delap - alphap).ge.1.d-10) then
                           alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                        else
                           alpham = -2.d0*alphap
                        endif
                     endif
                  end if
               end if

               sm(i,j) = s(i,j) + alpham
               sp(i,j) = s(i,j) + alphap

            end do
         end do

!c     different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
         if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
!c     the value in the first cc ghost cell represents the edge value
            sm(lo(1),lo(2)-1:hi(2)+1)    = s(lo(1)-1,lo(2)-1:hi(2)+1)
            sedgex(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

!c     use a modified stencil to get sedgex on the first interior edge
            sedgex(lo(1)+1,lo(2)-1:hi(2)+1) =&
                -(1.d0/5.d0)  *s(lo(1)-1,lo(2)-1:hi(2)+1)&
                +(3.d0/4.d0)  *s(lo(1)  ,lo(2)-1:hi(2)+1)&
                +0.5d0        *s(lo(1)+1,lo(2)-1:hi(2)+1)&
                -(1.d0/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

!c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(lo(1)+1,j) = max(sedgex(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
               sedgex(lo(1)+1,j) = min(sedgex(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
            end do

!c     copy sedgex into sp
            do j=lo(2)-1,hi(2)+1
               sp(lo(1)  ,j) = sedgex(lo(1)+1,j)
            end do

!c     apply Colella 2008 limiters to compute sm and sp in the second
!c     and third inner cells
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)+1,lo(1)+2

                  alphap = sedgex(i+1,j)-s(i,j)
                  alpham = sedgex(i  ,j)-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
!c     Possible extremum. We look at cell centered values and face
!c     centered values for a change in sign in the differences adjacent to
!c     the cell. We use the pair of differences whose minimum magnitude is the
!c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgex(i,j) - sedgex(i-1,j)
                     dafacep = sedgex(i+2,j) - sedgex(i+1,j)
                     dabarm = s(i,j) - s(i-1,j)
                     dabarp = s(i+1,j) - s(i,j)
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
                     D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i+1,j)+s(i+2,j)
                     D2C = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i-1,j) - s(i,j)
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
                        delap = s(i+1,j) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

         if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
!c     the value in the first cc ghost cell represents the edge value
            sp(hi(1),lo(2)-1:hi(2)+1)      = s(hi(1)+1,lo(2)-1:hi(2)+1)
            sedgex(hi(1)+1,lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

!c     use a modified stencil to get sedgex on the first interior edge
            sedgex(hi(1),lo(2)-1:hi(2)+1) =&
                -(1.d0/5.d0)  *s(hi(1)+1,lo(2)-1:hi(2)+1)&
                +(3.d0/4.d0)  *s(hi(1)  ,lo(2)-1:hi(2)+1)&
                +0.5d0        *s(hi(1)-1,lo(2)-1:hi(2)+1)&
                -(1.d0/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

!c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(hi(1),j) = max(sedgex(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
               sedgex(hi(1),j) = min(sedgex(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
            end do

!c     copy sedgex into sm
            do j=lo(2)-1,hi(2)+1
               sm(hi(1)  ,j) = sedgex(hi(1),j)
            end do

!c     reset sm on second interior edge
            do j=lo(2)-1,hi(2)+1
               sm(hi(1)-1,j) = sedgex(hi(1)-1,j)
            end do

!c     apply Colella 2008 limiters to compute sm and sp in the second
!c     and third inner cells
            do j=lo(2)-1,hi(2)+1
               do i=hi(1)-2,hi(1)-1

                  alphap = sedgex(i+1,j)-s(i,j)
                  alpham = sedgex(i  ,j)-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
!c     Possible extremum. We look at cell centered values and face
!c     centered values for a change in sign in the differences adjacent to
!c     the cell. We use the pair of differences whose minimum magnitude is the
!c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgex(i,j) - sedgex(i-1,j)
                     dafacep = sedgex(i+2,j) - sedgex(i+1,j)
                     dabarm = s(i,j) - s(i-1,j)
                     dabarp = s(i+1,j) - s(i,j)
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
                     D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i+1,j)+s(i+2,j)
                     D2C = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i-1,j) - s(i,j)
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
                        delap = s(i+1,j) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

      end if

!c     compute x-component of Ip and Im
      do j=lo(2)-1,hi(2)+1
         do i=lo(1)-1,hi(1)+1
            sigma = abs(u(i,j))*dt/dx(1)
            s6 = 6.0d0*s(i,j) - 3.0d0*(sm(i,j)+sp(i,j))
            if (u(i,j) .gt. eps) then
               Ipx(i,j) = sp(i,j) - (sigma/2.0d0)*&
                   (sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
               Imx(i,j) = s(i,j)
            else if (u(i,j) .lt. -eps) then
               Ipx(i,j) = s(i,j)
               Imx(i,j) = sm(i,j) + (sigma/2.0d0)*&
                   (sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
            else
               Ipx(i,j) = s(i,j)
               Imx(i,j) = s(i,j)
            end if
         end do
      end do

!c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c     y-direction
!c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!c     compute s at y-edges
      if (ppm_type .eq. 1) then

!c     compute van Leer slopes in y-direction
         dsvl = 0.d0
         do j=lo(2)-2,hi(2)+2
            do i=lo(1)-1,hi(1)+1
               dsc = 0.5d0 * (s(i,j+1) - s(i,j-1))
               dsl = 2.d0  * (s(i,j  ) - s(i,j-1))
               dsr = 2.d0  * (s(i,j+1) - s(i,j  ))
               if (dsl*dsr .gt. 0.d0)&
                   dsvl(i,j) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
            end do
         end do

!c     interpolate s to y-edges
         do j=lo(2)-1,hi(2)+2
            do i=lo(1)-1,hi(1)+1
               sedgey(i,j) = 0.5d0*(s(i,j)+s(i,j-1)) - (1.d0/6.d0)*(dsvl(i,j)-dsvl(i,j-1))
!c     make sure sedgey lies in between adjacent cell-centered values
               sedgey(i,j) = max(sedgey(i,j),min(s(i,j),s(i,j-1)))
               sedgey(i,j) = min(sedgey(i,j),max(s(i,j),s(i,j-1)))
            end do
         end do

!c     copy sedgey into sp and sm
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               sp(i,j) = sedgey(i,j+1)
               sm(i,j) = sedgey(i,j  )
            end do
         end do

!c     modify using quadratic limiters
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end do

!c     different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
         if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
!c     the value in the first cc ghost cell represents the edge value
            sm(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

!c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,lo(2)+1) =&
                -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,lo(2)-1)&
                +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,lo(2)  )&
                +0.5d0        *s(lo(1)-1:hi(1)+1,lo(2)+1)&
                -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

!c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,lo(2)+1) = max(sedgey(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
               sedgey(i,lo(2)+1) = min(sedgey(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
            end do

!c     copy sedgey into sp and sm
            do i=lo(1)-1,hi(1)+1
               sp(i,lo(2)  ) = sedgey(i,lo(2)+1)
               sm(i,lo(2)+1) = sedgey(i,lo(2)+1)
            end do

!c     reset sp on second interior edge
            do i=lo(1)-1,hi(1)+1
               sp(i,lo(2)+1) = sedgey(i,lo(2)+2)
            end do

!c     modify using quadratic limiters
            do i=lo(1)-1,hi(1)+1
               j = lo(2)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

         if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
!c     the value in the first cc ghost cell represents the edge value
            sp(lo(1)-1:hi(1)+1,hi(2)) = s(lo(1)-1:hi(1)+1,hi(2)+1)

!c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,hi(2)) =&
                -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,hi(2)+1)&
                +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,hi(2)  )&
                +0.5d0        *s(lo(1)-1:hi(1)+1,hi(2)-1)&
                -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

!c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,hi(2)) = max(sedgey(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
               sedgey(i,hi(2)) = min(sedgey(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
            end do

!c     copy sedgey into sp and sm
            do i=lo(1)-1,hi(1)+1
               sp(i,hi(2)-1) = sedgey(i,hi(2))
               sm(i,hi(2)  ) = sedgey(i,hi(2))
            end do

!c     reset sm on second interior edge
            do i=lo(1)-1,hi(1)+1
               sm(i,hi(2)-1) = sedgey(i,hi(2)-1)
            end do

!c     modify using quadratic limiters
            do i=lo(1)-1,hi(1)+1
               j = hi(2)-1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

      else if (ppm_type .eq. 2) then

!c     interpolate s to y-edges
         do j=lo(2)-2,hi(2)+3
            do i=lo(1)-1,hi(1)+1
               sedgey(i,j) = (7.d0/12.d0)*(s(i,j-1)+s(i,j)) - (1.d0/12.d0)*(s(i,j-2)+s(i,j+1))
!c     limit sedgey
               if ((sedgey(i,j)-s(i,j-1))*(s(i,j)-sedgey(i,j)) .lt. 0.d0) then
                  D2  = 3.d0*(s(i,j-1)-2.d0*sedgey(i,j)+s(i,j))
                  D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                  D2R = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                  sgn = sign(1.d0,D2)
                  D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                  sedgey(i,j) = 0.5d0*(s(i,j-1)+s(i,j)) - (1.d0/6.d0)*D2LIM
               end if
            end do
         end do

!c     use Colella 2008 limiters
!c     This is a new version of the algorithm
!c     to eliminate sensitivity to roundoff.
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1

               alphap = sedgey(i,j+1)-s(i,j)
               alpham = sedgey(i,j  )-s(i,j)
               bigp = abs(alphap).gt.2.d0*abs(alpham)
               bigm = abs(alpham).gt.2.d0*abs(alphap)
               extremum = .false.

               if (alpham*alphap .ge. 0.d0) then
                  extremum = .true.
               else if (bigp .or. bigm) then
!c     Possible extremum. We look at cell centered values and face
!c     centered values for a change in sign in the differences adjacent to
!c     the cell. We use the pair of differences whose minimum magnitude is the
!c     largest, and thus least susceptible to sensitivity to roundoff.
                  dafacem = sedgey(i,j) - sedgey(i,j-1)
                  dafacep = sedgey(i,j+2) - sedgey(i,j+1)
                  dabarm = s(i,j) - s(i,j-1)
                  dabarp = s(i,j+1) - s(i,j)
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
                  D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                  D2R = s(i,j)-2.d0*s(i,j+1)+s(i,j+2)
                  D2C = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                  sgn = sign(1.d0,D2)
                  D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                  alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                  alphap = alphap*D2LIM/max(abs(D2),1.d-10)
               else
                  if (bigp) then
                     sgn = sign(1.d0,alpham)
                     amax = -alphap**2 / (4*(alpham + alphap))
                     delam = s(i,j-1) - s(i,j)
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
                     delap = s(i,j+1) - s(i,j)
                     if (sgn*amax .ge. sgn*delap) then
                        if (sgn*(delap - alphap).ge.1.d-10) then
                           alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                        else
                           alpham = -2.d0*alphap
                        endif
                     endif
                  end if
               end if

               sm(i,j) = s(i,j) + alpham
               sp(i,j) = s(i,j) + alphap

            end do
         end do

!c     different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
         if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
!c     the value in the first c!c ghost cell represents the edge value
            sm(lo(1)-1:hi(1)+1,lo(2))    = s(lo(1)-1:hi(1)+1,lo(2)-1)
            sedgey(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

!c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,lo(2)+1) =&
                -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,lo(2)-1)&
                +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,lo(2)  )&
                +0.5d0        *s(lo(1)-1:hi(1)+1,lo(2)+1)&
                -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

!c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,lo(2)+1) = max(sedgey(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
               sedgey(i,lo(2)+1) = min(sedgey(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
            end do

!c     copy sedgey into sp
            do i=lo(1)-1,hi(1)+1
               sp(i,lo(2)  ) = sedgey(i,lo(2)+1)
            end do

!c     apply Colella 2008 limiters to compute sm and sp in the second
!c     and third inner cells
            do j=lo(2)+1,lo(2)+1
               do i=lo(1)-1,hi(1)+1

                  alphap = sedgey(i,j+1)-s(i,j)
                  alpham = sedgey(i,j  )-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
!c     Possible extremum. We look at cell centered values and face
!c     centered values for a change in sign in the differences adjacent to
!c     the cell. We use the pair of differences whose minimum magnitude is the
!c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgey(i,j) - sedgey(i,j-1)
                     dafacep = sedgey(i,j+2) - sedgey(i,j+1)
                     dabarm = s(i,j) - s(i,j-1)
                     dabarp = s(i,j+1) - s(i,j)
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
                     D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i,j+1)+s(i,j+2)
                     D2C = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i,j-1) - s(i,j)
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
                        delap = s(i,j+1) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

         if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
!c     the value in the first c!c ghost cell represents the edge value
            sp(lo(1)-1:hi(1)+1,hi(2))      = s(lo(1)-1:hi(1)+1,hi(2)+1)
            sedgey(lo(1)-1:hi(1)+1,hi(2)+1) = s(lo(1)-1:hi(1)+1,hi(2)+1)

!c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,hi(2)) =&
                -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,hi(2)+1)&
                +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,hi(2)  )&
                +0.5d0        *s(lo(1)-1:hi(1)+1,hi(2)-1)&
                -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

!c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,hi(2)) = max(sedgey(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
               sedgey(i,hi(2)) = min(sedgey(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
            end do

!c     copy sedgey into sm
            do i=lo(1)-1,hi(1)+1
               sm(i,hi(2)  ) = sedgey(i,hi(2))
            end do

!c     apply Colella 2008 limiters to compute sm and sp in the second
!c     and third inner cells
            do j=hi(2)-2,hi(2)-1
               do i=lo(1)-1,hi(1)+1

                  alphap = sedgey(i,j+1)-s(i,j)
                  alpham = sedgey(i,j  )-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
!c     Possible extremum. We look at cell centered values and face
!c     centered values for a change in sign in the differences adjacent to
!c     the cell. We use the pair of differences whose minimum magnitude is the
!c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgey(i,j) - sedgey(i,j-1)
                     dafacep = sedgey(i,j+2) - sedgey(i,j+1)
                     dabarm = s(i,j) - s(i,j-1)
                     dabarp = s(i,j+1) - s(i,j)
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
                     D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i,j+1)+s(i,j+2)
                     D2C = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i,j-1) - s(i,j)
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
                        delap = s(i,j+1) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

      end if

!c     compute y-component of Ip and Im
      do j=lo(2)-1,hi(2)+1
         do i=lo(1)-1,hi(1)+1
            sigma = abs(v(i,j))*dt/dx(2)
            s6 = 6.0d0*s(i,j) - 3.0d0*(sm(i,j)+sp(i,j))
            if (v(i,j) .gt. eps) then
               Ipy(i,j) = sp(i,j) - (sigma/2.0d0)*&
                   (sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
               Imy(i,j) = s(i,j)
            else if (v(i,j) .lt. -eps) then
               Ipy(i,j) = s(i,j)
               Imy(i,j) = sm(i,j) + (sigma/2.0d0)*&
                   (sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
            else
               Ipy(i,j) = s(i,j)
               Imy(i,j) = s(i,j)
            end if
         end do
      end do

    end subroutine ppm

      subroutine adv_forcing(&
          aofs,DIMS(aofs),&
          xflux,DIMS(xflux),&
          uedge,DIMS(uedge),&
          areax,DIMS(ax),&
          yflux,DIMS(yflux),&
          vedge,DIMS(vedge),&
          areay,DIMS(ay),&
          vol,DIMS(vol),&
          lo,hi,iconserv ) bind(C,name="adv_forcing")
! c
! c     This subroutine uses scalar edge states to compute
! c     an advective tendency
! c
      implicit none
      integer i,j
      integer iconserv
      REAL_T divux,divuy
      integer imin,jmin,imax,jmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(aofs)
      integer DIMDEC(vol)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      REAL_T aofs(DIMV(aofs))
      REAL_T vol(DIMV(vol))
      REAL_T uedge(DIMV(uedge))
      REAL_T vedge(DIMV(vedge))
      REAL_T xflux(DIMV(xflux))
      REAL_T yflux(DIMV(yflux))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
!c
!c     if nonconservative initialize the advective tendency as -U*grad(S)
!c

      if ( iconserv .ne. 1 ) then
         do j = jmin,jmax
            do i = imin,imax
               divux = (&
                   areax(i+1,j)*uedge(i+1,j) -&
                   areax(i,  j)*uedge(i,  j) )/vol(i,j)
               divuy = (&
                   areay(i,j+1)*vedge(i,j+1) -&
                   areay(i,j  )*vedge(i,j  ) )/vol(i,j)
               aofs(i,j) =&
                   - divux*half*(xflux(i+1,j) + xflux(i,j))&
                   - divuy*half*(yflux(i,j+1) + yflux(i,j))

            end do
         end do
      end if
!c
!c     convert edge states to fluxes
!c
      do j = jmin,jmax
         do i = imin,imax+1
            xflux(i,j) = xflux(i,j)*uedge(i,j)*areax(i,j)
         end do
      end do
      do j = jmin,jmax+1
         do i = imin,imax
            yflux(i,j) = yflux(i,j)*vedge(i,j)*areay(i,j)
         end do
      end do
!c
!c     compute part of advective tendency that depends on the flux convergence
!c
      if ( iconserv .ne. 1 ) then
         do j = jmin,jmax
            do i = imin,imax
               aofs(i,j) = aofs(i,j) + (&
                   xflux(i+1,j) - xflux(i,j) +&
                   yflux(i,j+1) - yflux(i,j))/vol(i,j)
            end do
         end do
      else
         do j = jmin,jmax
            do i = imin,imax
               aofs(i,j) = (&
                   xflux(i+1,j) - xflux(i,j) +&
                   yflux(i,j+1) - yflux(i,j))/vol(i,j)
            end do
         end do
      end if

    end subroutine adv_forcing

      subroutine sync_adv_forcing(&
          sync,DIMS(sync),&
          xflux,DIMS(xflux),&
          ucor,DIMS(ucor),&
          areax,DIMS(ax),&
          yflux,DIMS(yflux),&
          vcor,DIMS(vcor),&
          areay,DIMS(ay),&
          vol,DIMS(vol),&
          lo,hi ) bind(C,name="sync_adv_forcing")
! c
! c     This subroutine computes the sync advective tendency
! c     for a state variable
! c
      implicit none
      integer i,j
      integer imin,jmin,imax,jmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(sync)
      integer DIMDEC(vol)
      integer DIMDEC(ucor)
      integer DIMDEC(vcor)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      REAL_T sync(DIMV(sync))
      REAL_T vol(DIMV(vol))
      REAL_T ucor(DIMV(ucor))
      REAL_T vcor(DIMV(vcor))
      REAL_T xflux(DIMV(xflux))
      REAL_T yflux(DIMV(yflux))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
! c
! c     compute corrective fluxes from edge states
! c     and perform conservative update
! c
      do j = jmin,jmax
         do i = imin,imax+1
            xflux(i,j) = xflux(i,j)*ucor(i,j)*areax(i,j)
         end do
      end do
      do j = jmin,jmax+1
         do i = imin,imax
            yflux(i,j) = yflux(i,j)*vcor(i,j)*areay(i,j)
         end do
      end do

      do j = jmin,jmax
         do i = imin,imax
            sync(i,j) = sync(i,j) + (&
                 xflux(i+1,j) - xflux(i,j) +&
                yflux(i,j+1) - yflux(i,j) )/vol(i,j)
         end do
      end do

    end subroutine sync_adv_forcing


    subroutine ppm_fpu(lo,hi,&
         s,s_lo,s_hi,&
         uedge,uedge_lo,uedge_hi,&
         vedge,vedge_lo,vedge_hi,&
         Ipx,Ipx_lo,Ipx_hi,&
         Imx,Imx_lo,Imx_hi,&
         Ipy,Ipy_lo,Ipy_hi,&
         Imy,Imy_lo,Imy_hi,&
         sm,sm_lo,sm_hi,&
         sp,sp_lo,sp_hi,&
         dsvl,dsvl_lo,dsvl_hi,&
         sedgex,sedgex_lo,sedgex_hi,&
         sedgey,sedgey_lo,sedgey_hi,&
         dx, dt, bc, eps,ppm_type)

      implicit none

      integer, intent(in) :: bc(SDIM,2)
      integer, dimension(2), intent(in) :: &
           s_lo,s_hi,uedge_lo,uedge_hi,vedge_lo,vedge_hi,Ipx_lo,Ipx_hi,Imx_lo,Imx_hi,Ipy_lo,Ipy_hi,Imy_lo,Imy_hi,&
           sm_lo,sm_hi,sp_lo,sp_hi,dsvl_lo,dsvl_hi,sedgex_lo,sedgex_hi,sedgey_lo,sedgey_hi,&
           lo,hi

      real(rt), intent(in) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2))
      real(rt), intent(in) :: uedge(uedge_lo(1):uedge_hi(1),uedge_lo(2):uedge_hi(2))
      real(rt), intent(in) :: vedge(vedge_lo(1):vedge_hi(1),vedge_lo(2):vedge_hi(2))
      real(rt), intent(inout) :: Ipx(Ipx_lo(1):Ipx_hi(1),Ipx_lo(2):Ipx_hi(2))
      real(rt), intent(inout) :: Imx(Imx_lo(1):Imx_hi(1),Imx_lo(2):Imx_hi(2))
      real(rt), intent(inout) :: Ipy(Ipy_lo(1):Ipy_hi(1),Ipy_lo(2):Ipy_hi(2))
      real(rt), intent(inout) :: Imy(Imy_lo(1):Imy_hi(1),Imy_lo(2):Imy_hi(2))
      real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2))
      real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1),sp_lo(2):sp_hi(2))
      real(rt), intent(inout) :: dsvl(dsvl_lo(1):dsvl_hi(1),dsvl_lo(2):dsvl_hi(2))
      real(rt), intent(inout) :: sedgex(sedgex_lo(1):sedgex_hi(1),sedgex_lo(2):sedgex_hi(2))
      real(rt), intent(inout) :: sedgey(sedgey_lo(1):sedgey_hi(1),sedgey_lo(2):sedgey_hi(2))
      real(rt), intent(in) :: eps, dx(SDIM), dt
      integer ppm_type

!c     local
      integer i,j

      logical extremum, bigp, bigm

      real(rt) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, alphap, alpham
      real(rt) :: sgn, s6, amax, delam, delap, sigmam, sigmap
      real(rt) :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp
      real(rt), PARAMETER :: C = 1.25d0

!c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c     x-direction
!c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!c     compute s at x-edges
      if (ppm_type .eq. 1) then

!c     compute van Leer slopes in x-direction
         dsvl = 0.d0
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-2,hi(1)+2
               dsc = 0.5d0 * (s(i+1,j) - s(i-1,j))
               dsl = 2.d0  * (s(i  ,j) - s(i-1,j))
               dsr = 2.d0  * (s(i+1,j) - s(i  ,j))
               if (dsl*dsr .gt. 0.d0)&
                   dsvl(i,j) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
            end do
         end do

!c     interpolate s to x-edges
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+2
               sedgex(i,j) = 0.5d0*(s(i,j)+s(i-1,j)) - (1.d0/6.d0)*(dsvl(i,j)-dsvl(i-1,j))
!c     make sure sedgex lies in between adjacent cell-centered values
               sedgex(i,j) = max(sedgex(i,j),min(s(i,j),s(i-1,j)))
               sedgex(i,j) = min(sedgex(i,j),max(s(i,j),s(i-1,j)))
            end do
         end do

!c     copy sedgex into sp and sm
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               sp(i,j) = sedgex(i+1,j)
               sm(i,j) = sedgex(i  ,j)
            end do
         end do

!c     modify using quadrati!c limiters
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end do

!c     different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
         if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
!c     the value in the first c!c ghost cell represents the edge value
            sm(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

!c     use a modified stencil to get sedgex on the first interior edge
            sedgex(lo(1)+1,lo(2)-1:hi(2)+1) = &
                -(1.d0/5.d0)  *s(lo(1)-1,lo(2)-1:hi(2)+1)&
                +(3.d0/4.d0)  *s(lo(1)  ,lo(2)-1:hi(2)+1)&
                +0.5d0        *s(lo(1)+1,lo(2)-1:hi(2)+1)&
                -(1.d0/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

!c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(lo(1)+1,j) = max(sedgex(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
               sedgex(lo(1)+1,j) = min(sedgex(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
            end do

!c     copy sedgex into sp and sm
            do j=lo(2)-1,hi(2)+1
               sp(lo(1)  ,j) = sedgex(lo(1)+1,j)
               sm(lo(1)+1,j) = sedgex(lo(1)+1,j)
            end do

!c     reset sp on second interior edge
            do j=lo(2)-1,hi(2)+1
               sp(lo(1)+1,j) = sedgex(lo(1)+2,j)
            end do

!c     modify using quadrati!c limiters
            do j=lo(2)-1,hi(2)+1
               i = lo(1)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

         if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
!c     the value in the first c!c ghost cell represents the edge value
            sp(hi(1),lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

!c     use a modified stencil to get sedgex on the first interior edge
            sedgex(hi(1),lo(2)-1:hi(2)+1) =&
                -(1.d0/5.d0)  *s(hi(1)+1,lo(2)-1:hi(2)+1)&
                +(3.d0/4.d0)  *s(hi(1)  ,lo(2)-1:hi(2)+1)&
                +0.5d0        *s(hi(1)-1,lo(2)-1:hi(2)+1)&
                -(1.d0/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

!c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(hi(1),j) = max(sedgex(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
               sedgex(hi(1),j) = min(sedgex(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
            end do

!c     copy sedgex into sp and sm
            do j=lo(2)-1,hi(2)+1
               sp(hi(1)-1,j) = sedgex(hi(1),j)
               sm(hi(1)  ,j) = sedgex(hi(1),j)
            end do

!c     reset sm on second interior edge
            do j=lo(2)-1,hi(2)+1
               sm(hi(1)-1,j) = sedgex(hi(1)-1,j)
            end do

!c     modify using quadrati!c limiters
            do j=lo(2)-1,hi(2)+1
               i = hi(1)-1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

      else if (ppm_type .eq. 2) then

!c     interpolate s to x-edges
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-2,hi(1)+3
               sedgex(i,j) = (7.d0/12.d0)*(s(i-1,j)+s(i,j)) - (1.d0/12.d0)*(s(i-2,j)+s(i+1,j))
!c     limit sedgex
               if ((sedgex(i,j)-s(i-1,j))*(s(i,j)-sedgex(i,j)) .lt. 0.d0) then
                  D2  = 3.d0*(s(i-1,j)-2.d0*sedgex(i,j)+s(i,j))
                  D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                  D2R = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                  sgn = sign(1.d0,D2)
                  D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                  sedgex(i,j) = 0.5d0*(s(i-1,j)+s(i,j)) - (1.d0/6.d0)*D2LIM
               end if
            end do
         end do

!c     use Colella 2008 limiters
!c     This is a new version of the algorithm
!c     to eliminate sensitivity to roundoff.
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1

               alphap = sedgex(i+1,j)-s(i,j)
               alpham = sedgex(i  ,j)-s(i,j)
               bigp = abs(alphap).gt.2.d0*abs(alpham)
               bigm = abs(alpham).gt.2.d0*abs(alphap)
               extremum = .false.

               if (alpham*alphap .ge. 0.d0) then
                  extremum = .true.
               else if (bigp .or. bigm) then
!c     Possible extremum. We look at cell centered values and face
!c     centered values for a change in sign in the differences adjacent to
!c     the cell. We use the pair of differences whose minimum magnitude is the
!c     largest, and thus least susceptible to sensitivity to roundoff.
                  dafacem = sedgex(i,j) - sedgex(i-1,j)
                  dafacep = sedgex(i+2,j) - sedgex(i+1,j)
                  dabarm = s(i,j) - s(i-1,j)
                  dabarp = s(i+1,j) - s(i,j)
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
                  D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                  D2R = s(i,j)-2.d0*s(i+1,j)+s(i+2,j)
                  D2C = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                  sgn = sign(1.d0,D2)
                  D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                  alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                  alphap = alphap*D2LIM/max(abs(D2),1.d-10)
               else
                  if (bigp) then
                     sgn = sign(1.d0,alpham)
                     amax = -alphap**2 / (4*(alpham + alphap))
                     delam = s(i-1,j) - s(i,j)
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
                     delap = s(i+1,j) - s(i,j)
                     if (sgn*amax .ge. sgn*delap) then
                        if (sgn*(delap - alphap).ge.1.d-10) then
                           alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                        else
                           alpham = -2.d0*alphap
                        endif
                     endif
                  end if
               end if

               sm(i,j) = s(i,j) + alpham
               sp(i,j) = s(i,j) + alphap

            end do
         end do

!c     different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
         if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
!c     the value in the first c!c ghost cell represents the edge value
            sm(lo(1),lo(2)-1:hi(2)+1)    = s(lo(1)-1,lo(2)-1:hi(2)+1)
            sedgex(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

!c     use a modified stencil to get sedgex on the first interior edge
            sedgex(lo(1)+1,lo(2)-1:hi(2)+1) =&
                -(1.d0/5.d0)  *s(lo(1)-1,lo(2)-1:hi(2)+1)&
                +(3.d0/4.d0)  *s(lo(1)  ,lo(2)-1:hi(2)+1)&
                +0.5d0        *s(lo(1)+1,lo(2)-1:hi(2)+1)&
                -(1.d0/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

!c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(lo(1)+1,j) = max(sedgex(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
               sedgex(lo(1)+1,j) = min(sedgex(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
            end do

!c     copy sedgex into sp
            do j=lo(2)-1,hi(2)+1
               sp(lo(1)  ,j) = sedgex(lo(1)+1,j)
            end do

!c     apply Colella 2008 limiters to compute sm and sp in the second
!c     and third inner cells
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)+1,lo(1)+2

                  alphap = sedgex(i+1,j)-s(i,j)
                  alpham = sedgex(i  ,j)-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
!c     Possible extremum. We look at cell centered values and face
!c     centered values for a change in sign in the differences adjacent to
!c     the cell. We use the pair of differences whose minimum magnitude is the
!c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgex(i,j) - sedgex(i-1,j)
                     dafacep = sedgex(i+2,j) - sedgex(i+1,j)
                     dabarm = s(i,j) - s(i-1,j)
                     dabarp = s(i+1,j) - s(i,j)
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
                     D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i+1,j)+s(i+2,j)
                     D2C = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i-1,j) - s(i,j)
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
                        delap = s(i+1,j) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

         if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
!c     the value in the first c!c ghost cell represents the edge value
            sp(hi(1),lo(2)-1:hi(2)+1)      = s(hi(1)+1,lo(2)-1:hi(2)+1)
            sedgex(hi(1)+1,lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

!c     use a modified stencil to get sedgex on the first interior edge
            sedgex(hi(1),lo(2)-1:hi(2)+1) =&
                -(1.d0/5.d0)  *s(hi(1)+1,lo(2)-1:hi(2)+1)&
                +(3.d0/4.d0)  *s(hi(1)  ,lo(2)-1:hi(2)+1)&
                +0.5d0        *s(hi(1)-1,lo(2)-1:hi(2)+1)&
                -(1.d0/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

!c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(hi(1),j) = max(sedgex(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
               sedgex(hi(1),j) = min(sedgex(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
            end do

!c     copy sedgex into sm
            do j=lo(2)-1,hi(2)+1
               sm(hi(1)  ,j) = sedgex(hi(1),j)
            end do

!c     reset sm on second interior edge
            do j=lo(2)-1,hi(2)+1
               sm(hi(1)-1,j) = sedgex(hi(1)-1,j)
            end do

!c     apply Colella 2008 limiters to compute sm and sp in the second
!c     and third inner cells
            do j=lo(2)-1,hi(2)+1
               do i=hi(1)-2,hi(1)-1

                  alphap = sedgex(i+1,j)-s(i,j)
                  alpham = sedgex(i  ,j)-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
!c     Possible extremum. We look at cell centered values and face
!c     centered values for a change in sign in the differences adjacent to
!c     the cell. We use the pair of differences whose minimum magnitude is the
!c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgex(i,j) - sedgex(i-1,j)
                     dafacep = sedgex(i+2,j) - sedgex(i+1,j)
                     dabarm = s(i,j) - s(i-1,j)
                     dabarp = s(i+1,j) - s(i,j)
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
                     D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i+1,j)+s(i+2,j)
                     D2C = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i-1,j) - s(i,j)
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
                        delap = s(i+1,j) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

      end if

!c     compute x-component of Ip and Im
      do j=lo(2)-1,hi(2)+1
         do i=lo(1)-1,hi(1)+1
            sigmap = abs(uedge(i+1,j))*dt/dx(1)
            sigmam = abs(uedge(i,j))*dt/dx(1)
            s6 = 6.0d0*s(i,j) - 3.0d0*(sm(i,j)+sp(i,j))
            if (uedge(i+1,j) .gt. eps) then
               Ipx(i,j) = sp(i,j) - (sigmap/2.0d0)*&
                   (sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigmap)*s6)
            else
               Ipx(i,j) = s(i,j)
            end if
            if (uedge(i,j) .lt. -eps) then
               Imx(i,j) = sm(i,j) + (sigmam/2.0d0)*&
                   (sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigmam)*s6)
            else
               Imx(i,j) = s(i,j)
            end if
         end do
      end do

!c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c     y-direction
!c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!c     compute s at y-edges
      if (ppm_type .eq. 1) then

!c     compute van Leer slopes in y-direction
         dsvl = 0.d0
         do j=lo(2)-2,hi(2)+2
            do i=lo(1)-1,hi(1)+1
               dsc = 0.5d0 * (s(i,j+1) - s(i,j-1))
               dsl = 2.d0  * (s(i,j  ) - s(i,j-1))
               dsr = 2.d0  * (s(i,j+1) - s(i,j  ))
               if (dsl*dsr .gt. 0.d0)&
                   dsvl(i,j) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
            end do
         end do

!c     interpolate s to y-edges
         do j=lo(2)-1,hi(2)+2
            do i=lo(1)-1,hi(1)+1
               sedgey(i,j) = 0.5d0*(s(i,j)+s(i,j-1)) - (1.d0/6.d0)*(dsvl(i,j)-dsvl(i,j-1))
!c     make sure sedgey lies in between adjacent cell-centered values
               sedgey(i,j) = max(sedgey(i,j),min(s(i,j),s(i,j-1)))
               sedgey(i,j) = min(sedgey(i,j),max(s(i,j),s(i,j-1)))
            end do
         end do

!c     copy sedgey into sp and sm
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               sp(i,j) = sedgey(i,j+1)
               sm(i,j) = sedgey(i,j  )
            end do
         end do

!c     modify using quadrati!c limiters
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end do

!c     different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
         if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
!c     the value in the first c!c ghost cell represents the edge value
            sm(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

!c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,lo(2)+1) =&
                -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,lo(2)-1)&
                +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,lo(2)  )&
                +0.5d0        *s(lo(1)-1:hi(1)+1,lo(2)+1)&
                -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

!c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,lo(2)+1) = max(sedgey(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
               sedgey(i,lo(2)+1) = min(sedgey(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
            end do

!c     copy sedgey into sp and sm
            do i=lo(1)-1,hi(1)+1
               sp(i,lo(2)  ) = sedgey(i,lo(2)+1)
               sm(i,lo(2)+1) = sedgey(i,lo(2)+1)
            end do

!c     reset sp on second interior edge
            do i=lo(1)-1,hi(1)+1
               sp(i,lo(2)+1) = sedgey(i,lo(2)+2)
            end do

!c     modify using quadrati!c limiters
            do i=lo(1)-1,hi(1)+1
               j = lo(2)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

         if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
!c     the value in the first c!c ghost cell represents the edge value
            sp(lo(1)-1:hi(1)+1,hi(2)) = s(lo(1)-1:hi(1)+1,hi(2)+1)

!c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,hi(2)) =&
                -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,hi(2)+1)&
                +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,hi(2)  )&
                +0.5d0        *s(lo(1)-1:hi(1)+1,hi(2)-1)&
                -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

!c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,hi(2)) = max(sedgey(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
               sedgey(i,hi(2)) = min(sedgey(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
            end do

!c     copy sedgey into sp and sm
            do i=lo(1)-1,hi(1)+1
               sp(i,hi(2)-1) = sedgey(i,hi(2))
               sm(i,hi(2)  ) = sedgey(i,hi(2))
            end do

!c     reset sm on second interior edge
            do i=lo(1)-1,hi(1)+1
               sm(i,hi(2)-1) = sedgey(i,hi(2)-1)
            end do

!c     modify using quadrati!c limiters
            do i=lo(1)-1,hi(1)+1
               j = hi(2)-1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

      else if (ppm_type .eq. 2) then

!c     interpolate s to y-edges
         do j=lo(2)-2,hi(2)+3
            do i=lo(1)-1,hi(1)+1
               sedgey(i,j) = (7.d0/12.d0)*(s(i,j-1)+s(i,j)) - (1.d0/12.d0)*(s(i,j-2)+s(i,j+1))
!c     limit sedgey
               if ((sedgey(i,j)-s(i,j-1))*(s(i,j)-sedgey(i,j)) .lt. 0.d0) then
                  D2  = 3.d0*(s(i,j-1)-2.d0*sedgey(i,j)+s(i,j))
                  D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                  D2R = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                  sgn = sign(1.d0,D2)
                  D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                  sedgey(i,j) = 0.5d0*(s(i,j-1)+s(i,j)) - (1.d0/6.d0)*D2LIM
               end if
            end do
         end do

!c     use Colella 2008 limiters
!c     This is a new version of the algorithm
!c     to eliminate sensitivity to roundoff.
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1

               alphap = sedgey(i,j+1)-s(i,j)
               alpham = sedgey(i,j  )-s(i,j)
               bigp = abs(alphap).gt.2.d0*abs(alpham)
               bigm = abs(alpham).gt.2.d0*abs(alphap)
               extremum = .false.

               if (alpham*alphap .ge. 0.d0) then
                  extremum = .true.
               else if (bigp .or. bigm) then
!c     Possible extremum. We look at cell centered values and face
!c     centered values for a change in sign in the differences adjacent to
!c     the cell. We use the pair of differences whose minimum magnitude is the
!c     largest, and thus least susceptible to sensitivity to roundoff.
                  dafacem = sedgey(i,j) - sedgey(i,j-1)
                  dafacep = sedgey(i,j+2) - sedgey(i,j+1)
                  dabarm = s(i,j) - s(i,j-1)
                  dabarp = s(i,j+1) - s(i,j)
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
                  D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                  D2R = s(i,j)-2.d0*s(i,j+1)+s(i,j+2)
                  D2C = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                  sgn = sign(1.d0,D2)
                  D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                  alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                  alphap = alphap*D2LIM/max(abs(D2),1.d-10)
               else
                  if (bigp) then
                     sgn = sign(1.d0,alpham)
                     amax = -alphap**2 / (4*(alpham + alphap))
                     delam = s(i,j-1) - s(i,j)
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
                     delap = s(i,j+1) - s(i,j)
                     if (sgn*amax .ge. sgn*delap) then
                        if (sgn*(delap - alphap).ge.1.d-10) then
                           alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                        else
                           alpham = -2.d0*alphap
                        endif
                     endif
                  end if
               end if

               sm(i,j) = s(i,j) + alpham
               sp(i,j) = s(i,j) + alphap

            end do
         end do

!c     different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
         if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
!c     the value in the first c!c ghost cell represents the edge value
            sm(lo(1)-1:hi(1)+1,lo(2))    = s(lo(1)-1:hi(1)+1,lo(2)-1)
            sedgey(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

!c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,lo(2)+1) =&
                -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,lo(2)-1)&
                +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,lo(2)  )&
                +0.5d0        *s(lo(1)-1:hi(1)+1,lo(2)+1)&
                -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

!c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,lo(2)+1) = max(sedgey(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
               sedgey(i,lo(2)+1) = min(sedgey(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
            end do

!c     copy sedgey into sp
            do i=lo(1)-1,hi(1)+1
               sp(i,lo(2)  ) = sedgey(i,lo(2)+1)
            end do

!c     apply Colella 2008 limiters to compute sm and sp in the second
!c     and third inner cells
            do j=lo(2)+1,lo(2)+1
               do i=lo(1)-1,hi(1)+1

                  alphap = sedgey(i,j+1)-s(i,j)
                  alpham = sedgey(i,j  )-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
!c     Possible extremum. We look at cell centered values and face
!c     centered values for a change in sign in the differences adjacent to
!c     the cell. We use the pair of differences whose minimum magnitude is the
!c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgey(i,j) - sedgey(i,j-1)
                     dafacep = sedgey(i,j+2) - sedgey(i,j+1)
                     dabarm = s(i,j) - s(i,j-1)
                     dabarp = s(i,j+1) - s(i,j)
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
                     D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i,j+1)+s(i,j+2)
                     D2C = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i,j-1) - s(i,j)
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
                        delap = s(i,j+1) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

         if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
!c     the value in the first c!c ghost cell represents the edge value
            sp(lo(1)-1:hi(1)+1,hi(2))      = s(lo(1)-1:hi(1)+1,hi(2)+1)
            sedgey(lo(1)-1:hi(1)+1,hi(2)+1) = s(lo(1)-1:hi(1)+1,hi(2)+1)

!c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,hi(2)) =&
                -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,hi(2)+1)&
                +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,hi(2)  )&
                +0.5d0        *s(lo(1)-1:hi(1)+1,hi(2)-1)&
                -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

!c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,hi(2)) = max(sedgey(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
               sedgey(i,hi(2)) = min(sedgey(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
            end do

!c     copy sedgey into sm
            do i=lo(1)-1,hi(1)+1
               sm(i,hi(2)  ) = sedgey(i,hi(2))
            end do

!c     apply Colella 2008 limiters to compute sm and sp in the second
!c     and third inner cells
            do j=hi(2)-2,hi(2)-1
               do i=lo(1)-1,hi(1)+1

                  alphap = sedgey(i,j+1)-s(i,j)
                  alpham = sedgey(i,j  )-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
!c     Possible extremum. We look at cell centered values and face
!c     centered values for a change in sign in the differences adjacent to
!c     the cell. We use the pair of differences whose minimum magnitude is the
!c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgey(i,j) - sedgey(i,j-1)
                     dafacep = sedgey(i,j+2) - sedgey(i,j+1)
                     dabarm = s(i,j) - s(i,j-1)
                     dabarp = s(i,j+1) - s(i,j)
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
                     D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i,j+1)+s(i,j+2)
                     D2C = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i,j-1) - s(i,j)
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
                        delap = s(i,j+1) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

      end if

!c     compute y-component of Ip and Im
      do j=lo(2)-1,hi(2)+1
         do i=lo(1)-1,hi(1)+1
            sigmap = abs(vedge(i,j+1))*dt/dx(2)
            sigmam = abs(vedge(i,j))*dt/dx(2)
            s6 = 6.0d0*s(i,j) - 3.0d0*(sm(i,j)+sp(i,j))
            if (vedge(i,j+1) .gt. eps) then
               Ipy(i,j) = sp(i,j) - (sigmap/2.0d0)*&
                   (sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigmap)*s6)
            else
               Ipy(i,j) = s(i,j)
            end if

            if (vedge(i,j) .lt. -eps) then
               Imy(i,j) = sm(i,j) + (sigmam/2.0d0)*&
                   (sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigmam)*s6)
            else
               Imy(i,j) = s(i,j)
            end if
         end do
      end do

    end subroutine ppm_fpu

      subroutine convscalminmax(s,DIMS(s),sn,DIMS(sn),&
                                    lo,hi,bc) bind(C,name="convscalminmax")
!c
!c     correct an convectively-advected field for under/over shoots
!c
      implicit none
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer lo(SDIM), hi(SDIM)
      integer bc(SDIM,2)
      REAL_T s(DIMV(s))
      REAL_T sn(DIMV(sn))
      integer  i, j, imin, imax, jmin, jmax
      REAL_T   smin, smax

      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
!c
!c     compute extrema of s
!c
      do j = jmin, jmax
         do i = imin, imax
            smin = min(&
                s(i-1,j-1),s(i  ,j-1),s(i+1,j-1),&
                s(i-1,j  ),s(i  ,j  ),s(i+1,j  ),&
                s(i-1,j+1),s(i  ,j+1),s(i+1,j+1))
            smax = max(&
                s(i-1,j-1),s(i  ,j-1),s(i+1,j-1),&
                s(i-1,j  ),s(i  ,j  ),s(i+1,j  ),&
                s(i-1,j+1),s(i  ,j+1),s(i+1,j+1))
            sn(i,j) = max(sn(i,j),smin)
            sn(i,j) = min(sn(i,j),smax)
         end do
      end do

    end subroutine convscalminmax

      subroutine consscalminmax(s,rho,DIMS(s),sn,rhon,DIMS(sn),&
                                    lo,hi,bc) bind(C,name="consscalminmax")
!c
!c     correct an conservatively-advected field for under/over shoots
!c
      implicit none
      integer DIMDEC(s),DIMDEC(sn)
      integer lo(SDIM), hi(SDIM)
      integer bc(SDIM,2)
      REAL_T    s(DIMV(s))
      REAL_T  rho(DIMV(s))
      REAL_T   sn(DIMV(sn))
      REAL_T rhon(DIMV(sn))
      integer  i, j, imin, imax, jmin, jmax
      REAL_T   slo_min,smd_min,shi_min
      REAL_T   slo_max,smd_max,shi_max
      REAL_T   smin, smax

      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
!c
!c     compute extrema of s
!c
      do j = jmin-1, jmax+1
         do i = imin-1, imax+1
            s(i,j) = s(i,j) / rho(i,j)
         end do
      end do

      do j = jmin, jmax
         do i = imin, imax
            slo_min = min(s(i-1,j-1),s(i  ,j-1),s(i+1,j-1))
            smd_min = min(s(i-1,j  ),s(i  ,j  ),s(i+1,j  ))
            shi_min = min(s(i-1,j+1),s(i  ,j+1),s(i+1,j+1))
            smin = min(slo_min,smd_min,shi_min)

            slo_max = max(s(i-1,j-1),s(i  ,j-1),s(i+1,j-1))
            smd_max = max(s(i-1,j  ),s(i  ,j  ),s(i+1,j  ))
            shi_max = max(s(i-1,j+1),s(i  ,j+1),s(i+1,j+1))
            smax = max(slo_max,smd_max,shi_max)

            sn(i,j) = max(sn(i,j)/rhon(i,j),smin) * rhon(i,j)
            sn(i,j) = min(sn(i,j)/rhon(i,j),smax) * rhon(i,j)
         end do
      end do

      do j = jmin-1, jmax+1
         do i = imin-1, imax+1
            s(i,j) = s(i,j) * rho(i,j)
         end do
      end do

    end subroutine consscalminmax

      subroutine fort_sum_tf_gp(&
          tforces,DIMS(tf),&
          gp,DIMS(gp),&
          rho,DIMS(rho),&
          lo,hi ) bind(C,name="fort_sum_tf_gp")
!c
!c     sum pressure forcing into tforces
!c
      implicit none
      integer i, j, n
      integer DIMDEC(tf)
      integer DIMDEC(gp)
      integer DIMDEC(rho)
      integer lo(SDIM), hi(SDIM)
      REAL_T tforces(DIMV(tf),SDIM)
      REAL_T gp(DIMV(gp),SDIM)
      REAL_T rho(DIMV(rho))

      do n = 1, SDIM
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tforces(i,j,n) = (&
                   tforces(i,j,n)&
                   -    gp(i,j,n))/rho(i,j)
            end do
         end do
      end do

    end subroutine fort_sum_tf_gp

      subroutine fort_sum_tf_gp_visc(&
          tforces,DIMS(tf),&
          visc,DIMS(visc),&
          gp,DIMS(gp),&
          rho,DIMS(rho),&
          lo,hi ) bind(C,name="fort_sum_tf_gp_visc")
!c
!c     sum pressure forcing and viscous forcing into
!c     tforces
!c
      implicit none
      integer i, j, n
      integer DIMDEC(tf)
      integer DIMDEC(visc)
      integer DIMDEC(gp)
      integer DIMDEC(rho)
      integer lo(SDIM), hi(SDIM)
      REAL_T tforces(DIMV(tf),SDIM)
      REAL_T visc(DIMV(visc),SDIM)
      REAL_T gp(DIMV(gp),SDIM)
      REAL_T rho(DIMV(rho))

      do n = 1, SDIM
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

! EM_DEBUG
!if ((i < 4) .and.((j > 14).and.(j < 18))) then
!write(*,*) 'DEBUG IN SUM_TF_GP ',i,j,n,gp(i,j,n),tforces(i,j,n),visc(i,j,n),rho(i,j)
!endif

               tforces(i,j,n) = (&
                   tforces(i,j,n)&
                   +  visc(i,j,n)&
                   -    gp(i,j,n))/rho(i,j)

! EM_DEBUG
!if ((i < 4) .and.((j > 14).and.(j < 18))) then
!write(*,*) 'DEBUG IN SUM_TF_GP ',i,j,n,gp(i,j,n),tforces(i,j,n),visc(i,j,n),rho(i,j)
!endif

            end do
         end do
      end do

    end subroutine fort_sum_tf_gp_visc

      subroutine fort_sum_tf_divu(&
          s,DIMS(S),&
          tforces,DIMS(tf),&
          divu,DIMS(divu),&
          rho,DIMS(rho),&
          lo,hi,nvar,iconserv ) bind(C,name="fort_sum_tf_divu")
!c
!c     sum divU*S into tforces or divide tforces by rho
!c     depending on the value of iconserv
!c
      implicit none
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, n

      integer DIMDEC(S)
      integer DIMDEC(tf)
      integer DIMDEC(divu)
      integer DIMDEC(rho)

      REAL_T S(DIMV(S),nvar)
      REAL_T tforces(DIMV(tf),nvar)
      REAL_T divu(DIMV(divu))
      REAL_T rho(DIMV(rho))

      if ( iconserv .eq. 1 ) then
         do n = 1, nvar
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) = &
                 tforces(i,j,n) - S(i,j,n)*divu(i,j)
               end do
            end do
         end do
      else
         do n = 1, nvar
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) =&
                 tforces(i,j,n)/rho(i,j)
               end do
            end do
         end do
      end if

    end subroutine fort_sum_tf_divu

      subroutine fort_sum_tf_divu_visc(&
          S,DIMS(S),&
          tforces,DIMS(tf),&
          divu,DIMS(divu),&
          visc,DIMS(visc),&
          rho,DIMS(rho),&
          lo,hi,nvar,iconserv ) bind(C,name="fort_sum_tf_divu_visc")
!c
!c     sum tforces, viscous foricing and divU*S into tforces
!c     depending on the value of iconserv
!c
      implicit none
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, n

      integer DIMDEC(S)
      integer DIMDEC(tf)
      integer DIMDEC(divu)
      integer DIMDEC(visc)
      integer DIMDEC(rho)

      REAL_T S(DIMV(S),nvar)
      REAL_T tforces(DIMV(tf),nvar)
      REAL_T divu(DIMV(divu))
      REAL_T visc(DIMV(visc),nvar)
      REAL_T rho(DIMV(rho))

      if ( iconserv .eq. 1 ) then
         do n = 1, nvar
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) = &
                      tforces(i,j,n)&
                      +  visc(i,j,n)&
                      -     S(i,j,n)*divu(i,j)
               end do
            end do
         end do
      else
         do n = 1, nvar
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) = (&
                      tforces(i,j,n)&
                      +  visc(i,j,n) )/rho(i,j)
               end do
            end do
         end do
      end if

    end subroutine fort_sum_tf_divu_visc

      subroutine update_tf(&
          s,       DIMS(s),&
          sn,      DIMS(sn),&
          tforces, DIMS(tf),&
          lo,hi,dt,nvar) bind(C,name="update_tf")
!c
!c     update a field with a forcing term
!c
      implicit none
      integer i, j, n, nvar
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(tf)
      integer lo(SDIM), hi(SDIM)
      REAL_T dt
      REAL_T s(DIMV(s),nvar)
      REAL_T sn(DIMV(sn),nvar)
      REAL_T tforces(DIMV(tf),nvar)

      do n = 1,nvar
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               sn(i,j,n) = s(i,j,n) + dt*tforces(i,j,n)
            end do
         end do
      end do

    end subroutine update_tf

      subroutine update_aofs_tf(&
          s,       DIMS(s),&
          sn,      DIMS(sn),&
          aofs,    DIMS(aofs),&
          tforces, DIMS(tf),&
          lo,hi,dt,nvar) bind(C,name="update_aofs_tf")
!c
!c     update a field with an advective tendency
!c     and a forcing term
!c
      implicit none
      integer i, j, n, nvar
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(aofs)
      integer DIMDEC(tf)
      integer lo(SDIM), hi(SDIM)
      REAL_T dt
      REAL_T s(DIMV(s),nvar)
      REAL_T sn(DIMV(sn),nvar)
      REAL_T aofs(DIMV(aofs),nvar)
      REAL_T tforces(DIMV(tf),nvar)

      do n = 1,nvar
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               sn(i,j,n) = s(i,j,n)&
                   - dt*aofs(i,j,n)&
                   + dt*tforces(i,j,n)
            end do
         end do
      end do

    end subroutine update_aofs_tf

      subroutine update_aofs_tf_gp(&
          u,       DIMS(u),&
          un,      DIMS(un),&
          aofs,    DIMS(aofs),&
          tforces, DIMS(tf),&
          gp,      DIMS(gp),&
          rho,     DIMS(rho),&
          lo, hi, dt) bind(C,name="update_aofs_tf_gp")

!c
!c     update the velocities
!c
      implicit none
      integer i, j, n
      integer DIMDEC(u)
      integer DIMDEC(un)
      integer DIMDEC(aofs)
      integer DIMDEC(rho)
      integer DIMDEC(gp)
      integer DIMDEC(tf)
      integer lo(SDIM), hi(SDIM)
      REAL_T u(DIMV(u),SDIM)
      REAL_T un(DIMV(un),SDIM)
      REAL_T aofs(DIMV(aofs),SDIM)
      REAL_T rho(DIMV(rho))
      REAL_T gp(DIMV(gp),SDIM)
      REAL_T tforces(DIMV(tf),SDIM)
      REAL_T dt

      do n = 1, SDIM
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

! EM_DEBUG
!  if ((i < 2) .and.((j > 14).and.(j < 18))) then
!write(*,*) 'DEBUG IN update_aofs_tf_gp ',i,j,n,gp(i,j,n),rho(i,j),tforces(i,j,n),aofs(i,j,n),u(i,j,n)
!endif


               un(i,j,n) = u(i,j,n) &
                   - dt*   aofs(i,j,n)&
                   + dt*tforces(i,j,n)/rho(i,j)&
                   - dt*     gp(i,j,n)/rho(i,j)
            end do
         end do
      end do

    end subroutine update_aofs_tf_gp


      subroutine bdsslope(s,lo_1,lo_2,hi_1,hi_2,slx,sly,sc,dx)

      implicit none
      integer lo_1,lo_2,hi_1,hi_2
! C
! C     FIXME? Do the arrays s, slx, sly, sc have the same bounds
! C            here and in the calling function?
! C            Perhaps passing in loop bounds (is,js,ie,je) instead
! C            of infering from lo and hi would be better...
! C
      REAL_T      s(lo_1-3:hi_1+3,lo_2-3:hi_2+3)
      REAL_T    slx(lo_1-1:hi_1+1,lo_2-1:hi_2+1)
      REAL_T    sly(lo_1-1:hi_1+1,lo_2-1:hi_2+1)
      REAL_T   slxy(lo_1-1:hi_1+1,lo_2-1:hi_2+1)
      REAL_T   sint(lo_1-2:hi_1+2,lo_2-2:hi_2+2)
      REAL_T     sc(lo_1-1:hi_1+1,lo_2-1:hi_2+1,4)
      REAL_T dx(2)

      REAL_T   diff(lo_1-1:hi_1+1,4)
      REAL_T   smin(lo_1-1:hi_1+1,4)
      REAL_T   smax(lo_1-1:hi_1+1,4)
      REAL_T sumdif(lo_1-1:hi_1+1)
      REAL_T sgndif(lo_1-1:hi_1+1)
      integer   kdp(lo_1-1:hi_1+1)

      REAL_T hx,hy,sumloc,redfac,redmax,div
      REAL_T choice1,choice2
      REAL_T eps
      integer inc1, inc2, inc3, inc4
      integer i,j,k,ll,is,ie,js,je

      hx = dx(1)
      hy = dx(2)
      is = lo_1
      ie = hi_1
      js = lo_2
      je = hi_2

      eps = 1.d-8

      do i = is-2,ie+1
        do j = js-2,je+1
          sint(i,j) = (&
                   s(i-1,j-1) + s(i-1,j+2) + s(i+2,j-1) + s(i+2,j+2) &
          - seven*(s(i-1,j  ) + s(i-1,j+1) + s(i  ,j-1) + s(i+1,j-1) +&
                   s(i  ,j+2) + s(i+1,j+2) + s(i+2,j  ) + s(i+2,j+1)) +&
            49.d0*(s(i  ,j  ) + s(i+1,j  ) + s(i  ,j+1) + s(i+1,j+1)) ) / 144.d0
        enddo
      enddo

      do j = js-1,je+1
        do i = is-1,ie+1

          slx(i,j) = half*(sint(i  ,j) + sint(i  ,j-1) - &
                          sint(i-1,j) - sint(i-1,j-1) ) / hx
          sly(i,j) = half*(sint(i  ,j) - sint(i  ,j-1) + &
                          sint(i-1,j) - sint(i-1,j-1) ) / hy
          slxy(i,j) = (sint(i,j  ) - sint(i  ,j-1) - &
                      sint(i-1,j) + sint(i-1,j-1) ) / (hx*hy)
        enddo
      enddo

      do j = js-1,je+1

        do i = is-1,ie+1
          smin(i,4) = min(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
          smax(i,4) = max(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
          smin(i,3) = min(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
          smax(i,3) = max(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
          smin(i,2) = min(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
          smax(i,2) = max(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
          smin(i,1) = min(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
          smax(i,1) = max(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))

          sc(i,j,4) = s(i,j) + half*(hx*slx(i,j) + hy*sly(i,j)) &
                     + fourth*hx*hy*slxy(i,j)
          sc(i,j,3) = s(i,j) + half*(hx*slx(i,j) - hy*sly(i,j)) &
                     - fourth*hx*hy*slxy(i,j)
          sc(i,j,2) = s(i,j) - half*(hx*slx(i,j) - hy*sly(i,j)) &
                     - fourth*hx*hy*slxy(i,j)
          sc(i,j,1) = s(i,j) - half*(hx*slx(i,j) + hy*sly(i,j)) &
                     + fourth*hx*hy*slxy(i,j)

          sc(i,j,4) = max(min(sc(i,j,4), smax(i,4)), smin(i,4))
          sc(i,j,3) = max(min(sc(i,j,3), smax(i,3)), smin(i,3))
          sc(i,j,2) = max(min(sc(i,j,2), smax(i,2)), smin(i,2))
          sc(i,j,1) = max(min(sc(i,j,1), smax(i,1)), smin(i,1))
        enddo

        do ll = 1,3
          do i = is-1,ie+1
            sumloc = fourth*(sc(i,j,4) + sc(i,j,3) + &
                            sc(i,j,2) + sc(i,j,1))
            sumdif(i) = (sumloc - s(i,j))*4.
            sgndif(i) = sign(one,sumdif(i))

            diff(i,4) = (sc(i,j,4) - s(i,j))*sgndif(i)
            diff(i,3) = (sc(i,j,3) - s(i,j))*sgndif(i)
            diff(i,2) = (sc(i,j,2) - s(i,j))*sgndif(i)
            diff(i,1) = (sc(i,j,1) - s(i,j))*sgndif(i)

            inc1 = merge(1,0,(diff(i,1) - eps) .ge. 0.0d0)
            inc2 = merge(1,0,(diff(i,2) - eps) .ge. 0.0d0)
            inc3 = merge(1,0,(diff(i,3) - eps) .ge. 0.0d0)
            inc4 = merge(1,0,(diff(i,4) - eps) .ge. 0.0d0)
            kdp(i) = inc1 + inc2 + inc3 + inc4
          enddo

          do k = 1,4
            do i = is-1,ie+1
              div = merge(1,kdp(i),kdp(i) .lt. 1)
	      choice1 = sumdif(i)*sgndif(i)/div
              redfac = merge(choice1,zero,diff(i,k).gt.eps)
              kdp(i) = merge(kdp(i) - 1,kdp(i),diff(i,k) .gt. eps)
              choice1 = sc(i,j,k) - smin(i,k)
              choice2 = smax(i,k) - sc(i,j,k)
              redmax = merge(choice1,choice2,sgndif(i) .ge. 0.0d0)
              redfac = min(redfac,redmax)
              sumdif(i) = sumdif(i) - redfac*sgndif(i)
              sc(i,j,k) = sc(i,j,k) - redfac*sgndif(i)
            enddo
          enddo
        enddo

        do i = is-1,ie+1
          slx(i,j) = half*( sc(i,j,4) + sc(i,j,3)&
                          -sc(i,j,1) - sc(i,j,2))/hx
          sly(i,j) = half*( sc(i,j,4) + sc(i,j,2)&
                          -sc(i,j,1) - sc(i,j,3))/hy
          slxy(i,j) = ( sc(i,j,4) + sc(i,j,1)&
                      -sc(i,j,2) - sc(i,j,3) ) / (hx*hy)
        enddo
      enddo

      return
    end subroutine bdsslope

      subroutine estate_bds(s, tforces, divu, DIMS(s),&
          xlo, xhi, sx, slxscr, stxlo, stxhi,&
          uedge, DIMS(uedge), xstate, DIMS(xstate),&

          ylo, yhi, sy, slyscr, stylo, styhi,&
          vedge, DIMS(vedge), ystate, DIMS(ystate),&

          DIMS(work),&
          bc,lo,hi,dt,dx,n,use_minion,iconserv) bind(C,name="estate_bds")
! c
! c     This subroutine computes edges states, right now it uses
! c     a lot of memory, but there becomes a trade off between
! c     simplicity-efficiency in the new way of computing states
! c     and complexity in the old way.  By eliminating loops over
! c     state components though, the new way uses much less memory.
! c
      implicit none
      integer i,j,n
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer is,js,ie,je
      REAL_T hx, hy, dt, dth, dthx, dthy
      REAL_T dx(SDIM)
      REAL_T eps,eps_for_bc
      parameter( eps        = 1.0D-6 )
      parameter( eps_for_bc = 1.0D-10 )

      integer DIMDEC(s)
      integer DIMDEC(work)
      integer DIMDEC(uedge)
      integer DIMDEC(xstate)
      integer DIMDEC(vedge)
      integer DIMDEC(ystate)

      REAL_T s(DIMV(s))
      REAL_T stxlo(DIM1(s)),stxhi(DIM1(s)),slxscr(DIM1(s),4)
      REAL_T stylo(DIM2(s)),styhi(DIM2(s)),slyscr(DIM2(s),4)

      REAL_T uedge(DIMV(uedge)), xstate(DIMV(xstate))
      REAL_T vedge(DIMV(vedge)), ystate(DIMV(ystate))

      REAL_T xlo(DIMV(work)), xhi(DIMV(work))
      REAL_T ylo(DIMV(work)), yhi(DIMV(work))
      REAL_T  sx(DIMV(work))
      REAL_T  sy(DIMV(work))
      REAL_T tforces(DIMV(work))
      REAL_T    divu(DIMV(work))

      REAL_T smin,smax

      REAL_T   sc(DIMV(work),4)
      REAL_T gamp(DIM1(work))
      REAL_T gamm(DIM1(work))
      REAL_T   xm(DIM1(work))
      REAL_T   ym(DIM1(work))
      REAL_T    c(DIM1(work),4)
      REAL_T dt3rd

      integer choose_hi(DIM1(work))

      integer use_minion, iconserv
      integer iup,isign,jup,jsign
      REAL_T vtrans,vmult,vdif,stem,vaddif

      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      hx   = dx(1)
      hy   = dx(2)
      is = lo(1)
      js = lo(2)
      ie = hi(1)
      je = hi(2)

      dt3rd = dt * third

      smax = s(is,js)
      smin = s(is,js)
      do j = js,je
        do i = is,ie
          smax = max(smax,s(i,j))
          smin = min(smin,s(i,j))
        enddo
      enddo

!C see comments in subroutine bdsslope
      call bdsslope(s,lo(1),lo(2),hi(1),hi(2),sx,sy,sc,dx)

      do j = js,je

        do i = is-1,ie
          choose_hi(i+1) = merge(0,1,uedge(i+1,j) .ge. 0.0d0)
        end do
        if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP&
               .or.bc(1,1).eq.REFLECT_EVEN) then
           choose_hi(is) = 1
        endif
        if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP&
               .or.bc(1,2).eq.REFLECT_EVEN) then
           choose_hi(ie+1) = 0
        endif

        do i = is-1,ie
          iup   = merge(i,i+1,choose_hi(i+1).eq.0)
          isign = merge(1,-1, choose_hi(i+1).eq.0)
          vtrans = vedge(iup,j+1)
          jup   = merge(j,j+1,vtrans .ge. 0.0d0)
          jsign = merge(1,-1,vtrans .ge. 0.0d0)
          vmult = merge(uedge(i+1,j),half*(uedge(i+1,j)+uedge(i+1,j+1)),vtrans .ge. 0.0d0)
          xm(i) = isign*half*hx - two3rd*dt*vmult
          xm(i) = merge(xm(i), isign*min(isign*xm(i),hx*half), vtrans .ge. 0.0d0)
          ym(i) = jsign*half*hy - dt3rd*vedge(iup,j+1)

          c(i,1) = sc(iup,jup,1)
          c(i,2) = sc(iup,jup,2)
          c(i,3) = sc(iup,jup,3)
          c(i,4) = sc(iup,jup,4)

        enddo

        call bilin(gamp,c,xm,ym,hx,hy,is-1,ie,lo(1),hi(1))

!c       Impose BCs on gamp.
        if (j .eq. je) then
          do i = is-1, ie
            if (bc(2,2).eq.EXT_DIR .and.(vedge(i,je+1).le.zero) ) then
               gamp(i) = s(i,je+1)
            else if (bc(2,2).eq.REFLECT_ODD) then
               gamp(i) = zero
            end if
          end do
        end if

        do i = is-1,ie
          iup   = merge(i,i+1,choose_hi(i+1).eq.0)
          isign = merge(1,-1, choose_hi(i+1).eq.0)
          vtrans = vedge(iup,j)
          jup   = merge(j-1,j,vtrans .ge. 0.0d0)
          jsign = merge(1,-1,vtrans .ge. 0.0d0)
          vmult = merge(half*(uedge(i+1,j)+uedge(i+1,j-1)),uedge(i+1,j),vtrans .ge. 0.0d0)
          xm(i) = isign*half*hx - two3rd*dt*vmult
          xm(i) = merge(isign*min(isign*xm(i),hx*half), xm(i), vtrans .ge. 0.0d0)
          ym(i) = jsign*half*hy - dt3rd*vedge(iup,j)

          c(i,1) = sc(iup,jup,1)
          c(i,2) = sc(iup,jup,2)
          c(i,3) = sc(iup,jup,3)
          c(i,4) = sc(iup,jup,4)

        enddo

        call bilin(gamm,c,xm,ym,hx,hy,is-1,ie,lo(1),hi(1))

!c       Impose BCs on gamm.
        if (j .eq. js) then
          do i = is-1, ie
            if (bc(2,1).eq.EXT_DIR .and.(vedge(i,js).ge.zero) ) then
               gamm(i) = s(i,js-1)
            else if (bc(2,1).eq.REFLECT_ODD) then
               gamm(i) = zero
            end if
          end do
        end if

        do i = is-1, ie
          iup   = merge(i,i+1,choose_hi(i+1).eq.0)
          isign = merge(1,-1, choose_hi(i+1).eq.0)
          vdif   = dth*(vedge(iup,j+1)*gamp(i) - vedge(iup,j  )*gamm(i) ) / hy
          stem   = s(iup,j) + (isign*hx - uedge(i+1,j)*dt)*half*sx(iup,j)

!c         vaddif = dth*stem*(uedge(iup+1,j) - uedge(iup,j))/hx
!c         vaddif = dth*stem*(divu(iup,j) - (vedge(iup,j+1)-vedge(iup,j))/hy)
!c         Change per JBB - dont use stem here!
          vaddif = dth*s(iup,j)*(divu(iup,j) - (vedge(iup,j+1)-vedge(iup,j))/hy)

          xstate(i+1,j) = stem - vdif - vaddif + dth*tforces(iup,j)

          if (iconserv .eq. 0) &
           xstate(i+1,j) = xstate(i+1,j) + dth*stem*divu(iup,j)

        enddo
      enddo

      do j = js-1,je

        do i = is,ie
          choose_hi(i) = merge(0,1,vedge(i,j+1) .ge. 0.0d0)
          if ( (j .eq. js-1) .and. &
            (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP&
                                .or.bc(2,1).eq.REFLECT_EVEN) )&
          choose_hi(i) = 1
          if ( (j .eq. je) .and. &
            (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP&
                                .or.bc(2,2).eq.REFLECT_EVEN) )&
          choose_hi(i) = 0
        end do

        do i = is,ie
          jup   = merge(j,j+1,choose_hi(i).eq.0)
          jsign = merge(1, -1,choose_hi(i).eq.0)
          vtrans = uedge(i+1,jup)
          iup   = merge(i,i+1,vtrans .ge. 0.0d0)
          isign = merge(1,-1,vtrans .ge. 0.0d0)
          vmult = merge(vedge(i,j+1),half*(vedge(i,j+1)+vedge(i+1,j+1)),vtrans .ge. 0.0d0)
          xm(i) = isign*half*hx - dt3rd*uedge(i+1,jup)
          ym(i) = jsign*half*hy - two3rd*dt*vmult
          ym(i) = merge(ym(i),jsign*min(jsign*ym(i),hy*half),vtrans .ge. 0.0d0)

          c(i,1) = sc(iup,jup,1)
          c(i,2) = sc(iup,jup,2)
          c(i,3) = sc(iup,jup,3)
          c(i,4) = sc(iup,jup,4)
        enddo

        call bilin(gamp,c,xm,ym,hx,hy,is,ie,lo(1),hi(1))

!c       Impose BCs on gamp.
        if (bc(1,2).eq.EXT_DIR .and.(uedge(ie+1,j).le.zero) ) then
           gamp(ie) = s(ie+1,j)
        else if (bc(1,2).eq.REFLECT_ODD) then
           gamp(ie) = zero
        end if

        do i = is,ie
          jup   = merge(j,j+1,choose_hi(i).eq.0)
          jsign = merge(1, -1,choose_hi(i).eq.0)
          vtrans = uedge(i,jup)
          iup   = merge(i-1,i,vtrans .ge. 0.0d0)
          isign = merge(1,-1,vtrans .ge. 0.0d0)
          vmult = merge(half*(vedge(i,j+1)+vedge(i-1,j+1)),vedge(i,j+1),vtrans .ge. 0.0d0)
          xm(i) = isign*half*hx - dt3rd*uedge(i,jup)
          ym(i) = jsign*half*hy - two3rd*dt*vmult
          ym(i) = merge(jsign*min(jsign*ym(i),hy*half),ym(i),vtrans .ge. 0.0d0)

          c(i,1) = sc(iup,jup,1)
          c(i,2) = sc(iup,jup,2)
          c(i,3) = sc(iup,jup,3)
          c(i,4) = sc(iup,jup,4)
        enddo

        call bilin(gamm,c,xm,ym,hx,hy,is,ie,lo(1),hi(1))

!c       Impose BCs on gamm.
        if (bc(1,1).eq.EXT_DIR .and.(uedge(is,j).ge.zero) ) then
           gamm(is) = s(is,j)
        else if (bc(1,1).eq.REFLECT_ODD) then
           gamm(is) = zero
        end if

        do i = is,ie
          jup   = merge(j,j+1,choose_hi(i).eq.0)
          jsign = merge(1, -1,choose_hi(i).eq.0)
          vdif   = dth*(uedge(i+1,jup)*gamp(i)-uedge(i,jup)*gamm(i))/hx
          stem   = s(i,jup) + (jsign*hy - vedge(i,j+1)*dt)*half*sy(i,jup)

!c         vaddif = dth*stem*(vedge(i,jup+1) - vedge(i,jup))/hy
!c         vaddif = dth*stem*(divu(i,jup) - (uedge(i+1,jup)-uedge(i,jup))/hx)
!c         Change per JBB - dont use stem here!
          vaddif = dth*s(i,jup)*(divu(i,jup) - (uedge(i+1,jup)-uedge(i,jup))/hx)

          ystate(i,j+1) = stem - vdif - vaddif + dth*tforces(i,jup)

          if (iconserv .eq. 0) &
           ystate(i,j+1) = ystate(i,j+1) + dth*stem*divu(i,jup)

        enddo
      enddo

!c
!c     Impose BCs on xstate.
!c
      do j = js,je
        if (bc(1,1).eq.EXT_DIR .and. uedge(is,j).ge.zero) then
           xstate(is,j) = s(is-1,j)
        else if (bc(1,1).eq.REFLECT_ODD) then
           xstate(is,j) = zero
        end if

        if (bc(1,2).eq.EXT_DIR .and. uedge(ie+1,j).le.zero) then
           xstate(ie+1,j) = s(ie+1,j)
        else if (bc(1,2).eq.REFLECT_ODD) then
           xstate(ie+1,j) = zero
        end if
      enddo

!c
!c     Impose BCs on ystate.
!c
      do i = is,ie
        if (bc(2,1).eq.EXT_DIR .and.(vedge(i,js).ge.zero) ) then
           ystate(i,js) = s(i,js-1)
        else if (bc(2,1).eq.REFLECT_ODD) then
            ystate(i,js) = zero
        end if

        if (bc(2,2).eq.EXT_DIR .and. (vedge(i,je+1).le.zero) ) then
           ystate(i,je+1) = s(i,je+1)
        else if (bc(2,2).eq.REFLECT_ODD) then
            ystate(i,je+1) = zero
        end if
      enddo

    end subroutine estate_bds

!C FIXME? check that the arrays are properly sized
      subroutine bilin(poly,c,x,y,hx,hy,istart,iend,lo_1,hi_1)

      implicit none
      integer lo_1,hi_1
      integer istart,iend
      REAL_T poly(lo_1-1:hi_1+1)
      REAL_T    x(lo_1-1:hi_1+1)
      REAL_T    y(lo_1-1:hi_1+1)
      REAL_T    c(lo_1-1:hi_1+1,4)
      REAL_T hx, hy
      REAL_T centx, centy
      REAL_T l1,l2,l3,l4

      integer i

      centx = hx/2
      centy = hy/2

      do i = istart,iend

        l1 = (-centx - x(i))*(-centy - y(i))
        l2 = (-centx - x(i))*( centy - y(i))
        l3 = ( centx - x(i))*(-centy - y(i))
        l4 = ( centx - x(i))*( centy - y(i))

        poly(i) = l2*c(i,3) + l3*c(i,2) - l1*c(i,4) - l4*c(i,1)
        poly(i) = -poly(i)/(hx*hy)
      enddo

    end subroutine bilin

 end module godunov_2d_module
