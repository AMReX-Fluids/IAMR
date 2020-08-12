
#include <AMReX_BC_TYPES.H>
#include <AMReX_FluxRegister.H>
#include <SyncRegister.H>
#include <AMReX_BLProfiler.H>
#include <NSB_K.H>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

SyncRegister::SyncRegister (const BoxArray& fine_boxes,
                            const DistributionMapping& dmap,
                            const IntVect&  ref_ratio)
    : ratio(ref_ratio)
{
    BL_ASSERT(grids.size() == 0);
    BL_ASSERT(fine_boxes.isDisjoint());

    grids = fine_boxes;
    grids.coarsen(ratio);

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        Orientation loface(dir, Orientation::low);
        Orientation hiface(dir, Orientation::high);

        BndryBATransformer lotrans(loface, IndexType::TheNodeType(), 0, 1, 0);
        BndryBATransformer hitrans(hiface, IndexType::TheNodeType(), 0, 1, 0);

        BoxArray loBA(grids, lotrans);
        BoxArray hiBA(grids, hitrans);

        bndry[loface].define(loBA,dmap,1);
        bndry_mask[loface].define(loBA,dmap,1);
        bndry[hiface].define(hiBA,dmap,1);
        bndry_mask[hiface].define(hiBA,dmap,1);
    }
}

SyncRegister::~SyncRegister () {}

void /* note that rhs is on a different BoxArray */
SyncRegister::InitRHS (MultiFab& rhs, const Geometry& geom, const BCRec& phys_bc)
{
    BL_PROFILE("SyncRegister::InitRHS()");

    BL_ASSERT(rhs.nComp() == 1);
 
    int ngrow = rhs.nGrow();

    rhs.setVal(0);

    for (OrientationIter face; face; ++face)
    {
        bndry[face()].copyTo(rhs,ngrow,0,0,bndry[face()].nComp(),geom.periodicity());
    }

    const int* phys_lo = phys_bc.lo();
    const int* phys_hi = phys_bc.hi();

    int outflow_dirs[AMREX_SPACEDIM-1]={-1};
    int nOutflow = 0;
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
      if (phys_lo[dir] == Outflow || phys_hi[dir] == Outflow)
      {
         outflow_dirs[nOutflow]=dir;
         nOutflow++;
      }
    }

    const Box& node_domain = amrex::surroundingNodes(geom.Domain());

    if (nOutflow > 0)
    {      
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (int n = 0; n < nOutflow; n++)
      {
         int dir = outflow_dirs[n];

         for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
         {
           const Box& vbx = mfi.validbox();
           auto const& rhs_arr = rhs.array(mfi);

           if (phys_lo[dir] == Outflow)
           {
             Box domlo(node_domain);
             domlo.setRange(dir,node_domain.smallEnd(dir),1);
             const Box& blo = vbx & domlo;

             if (blo.ok()) {
               amrex::ParallelFor(blo, [rhs_arr]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  rhs_arr(i,j,k) = 0.0;
               });
             }  
           }
           if (phys_hi[dir] == Outflow)
           {
             Box domhi(node_domain);
             domhi.setRange(dir,node_domain.bigEnd(dir),1);
             const Box& bhi = vbx & domhi;
             if (bhi.ok()) {
               amrex::ParallelFor(bhi, [rhs_arr]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  rhs_arr(i,j,k) = 0.0;
               });
             }  
           }
         }
      }
    }

    //
    // Set Up bndry_mask.
    //
    for (OrientationIter face; face; ++face)
    {
        bndry_mask[face()].setVal(0);
    }
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
       for (OrientationIter face_it; face_it; ++face_it)
       {
          FabSet& fs = bndry_mask[face_it()];

          FArrayBox tmpfab;
          std::vector< std::pair<int,Box> > isects;	    
          Vector<IntVect> pshifts(26);

          for (FabSetIter fsi(fs); fsi.isValid(); ++fsi)
          {
             FArrayBox& fab = fs[fsi];
             Box mask_cells = amrex::enclosedCells(amrex::grow(fab.box(),1));

             tmpfab.resize(mask_cells,1);
             Elixir tmpfab_i = tmpfab.elixir();
             const auto& tmp_arr = tmpfab.array();
             amrex::ParallelFor(mask_cells, [tmp_arr]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                  tmp_arr(i,j,k) = 0.0;
             });

             grids.intersections(mask_cells,isects);

             for (int is = 0, N = isects.size(); is < N; is++)
             {
                amrex::ParallelFor(isects[is].second, [tmp_arr]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                     tmp_arr(i,j,k) = 1.0;
                });
             }

             if (geom.isAnyPeriodic() && !geom.Domain().contains(mask_cells))
             {
                geom.periodicShift(geom.Domain(),mask_cells,pshifts);

                for (Vector<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end(); it != End; ++it)
                {
                   const IntVect& iv = *it;

                   grids.intersections(mask_cells+iv,isects);

                   for (int is = 0, N = isects.size(); is < N; is++)
                   {
                      Box& isect = isects[is].second;
                      isect     -= iv;
                      amrex::ParallelFor(isect, [tmp_arr]
                      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                      {
                           tmp_arr(i,j,k) = 1.0;
                      });
                   }
                }
             }
             const Box& bx = fab.box();
             auto const& mask = fab.array();  
             amrex::ParallelFor(bx, [mask,tmp_arr]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                  mask(i,j,k) = D_TERM(  tmp_arr(i  ,j  ,k  ) + tmp_arr(i-1,j  ,k  ),
                                       + tmp_arr(i  ,j-1,k  ) + tmp_arr(i-1,j-1,k  ),
                                       + tmp_arr(i  ,j  ,k-1) + tmp_arr(i-1,j  ,k-1)
                                       + tmp_arr(i  ,j-1,k-1) + tmp_arr(i-1,j-1,k-1));
             });
          }
       }
    }

    //
    // Here double the cell contributions if at a non-periodic physical bdry.
    //
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        if (!geom.isPeriodic(dir))
        {
            Box domlo(node_domain), domhi(node_domain);

            domlo.setRange(dir,node_domain.smallEnd(dir),1);
            domhi.setRange(dir,node_domain.bigEnd(dir),1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (OrientationIter face_it; face_it; ++face_it)
            {
                FabSet& fs = bndry_mask[face_it()];

                for (FabSetIter fsi(fs); fsi.isValid(); ++fsi)
                {
                    FArrayBox& fab = fs[fsi];
                    auto const& fab_arr = fs.array(fsi);

                    const Box& blo = fab.box() & domlo;

                    if (blo.ok()) {
                      amrex::ParallelFor(blo, [fab_arr]
                      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                      {
                         fab_arr(i,j,k) *= 2.0;
                      });
                    }

                    const Box& bhi = fab.box() & domhi;

                    if (bhi.ok()) {
                      amrex::ParallelFor(bhi, [fab_arr]
                      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                      {
                         fab_arr(i,j,k) *= 2.0;
                      });
                    }
                }
            }
        }
    }
    //
    // Here convert from sum of cell contributions to 0 or 1.
    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (OrientationIter face_it; face_it; ++face_it)
    {
        FabSet& fs = bndry_mask[face_it()];

        for (FabSetIter fsi(fs); fsi.isValid(); ++fsi)
        {
            FArrayBox& fab   = fs[fsi];
            const Box& bx    = fab.box();
            auto const& mask = fab.array();   
            const Real maxcount = D_TERM(AMREX_SPACEDIM,*AMREX_SPACEDIM,*AMREX_SPACEDIM) - 0.5;
            amrex::ParallelFor(bx, [mask,maxcount]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               mask(i,j,k) = (mask(i,j,k) > maxcount) ? 0.0 : 1.0;
            });
        }
    }

    // Multiply by Bndry Mask

    MultiFab tmp(rhs.boxArray(), rhs.DistributionMap(), 1, ngrow);

    for (OrientationIter face; face; ++face)
    {
        BL_ASSERT(bndry_mask[face()].nComp() == 1);
        tmp.setVal(1.0);
        bndry_mask[face()].copyTo(tmp, ngrow, 0, 0, 1);
        MultiFab::Multiply(rhs, tmp, 0, 0, 1, ngrow);
    }
}

void
SyncRegister::CrseInit (MultiFab& Sync_resid_crse, const Geometry& crse_geom, Real mult)
{
    BL_PROFILE("SyncRegister::CrseInit()");

    setVal(0.0);

    Sync_resid_crse.mult(mult);

    for (OrientationIter face; face; ++face)
    {
        bndry[face()].plusFrom(Sync_resid_crse,0,0,0,1,crse_geom.periodicity());
    }
}

void
SyncRegister::CompAdd (MultiFab& Sync_resid_fine,
                       const Geometry& fine_geom, const Geometry& crse_geom,
                       const BoxArray& Pgrids, Real mult)
{
    BL_PROFILE("SyncRegister::CompAdd()");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      Vector<IntVect> pshifts(27);
      std::vector< std::pair<int,Box> > isects;

      for (MFIter mfi(Sync_resid_fine); mfi.isValid(); ++mfi)
      {
         const Box& sync_box = mfi.validbox();

         Pgrids.intersections(sync_box,isects);
         const auto& sync_arr = Sync_resid_fine.array(mfi);

         for (int ii = 0, N = isects.size(); ii < N; ii++)
         {
             const Box& pbx = Pgrids[isects[ii].first];
             amrex::ParallelFor(isects[ii].second, [sync_arr]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                sync_arr(i,j,k) = 0.0;
             });
             fine_geom.periodicShift(sync_box, pbx, pshifts);

             for (Vector<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end(); it != End; ++it)
             {
                Box isect = pbx + *it;
                isect    &= sync_box;
                amrex::ParallelFor(isect, [sync_arr]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                   sync_arr(i,j,k) = 0.0;
                });
             }
         }
      }
    }

    FineAdd(Sync_resid_fine,crse_geom,mult);
}

void
SyncRegister::FineAdd (MultiFab& Sync_resid_fine, const Geometry& crse_geom, Real mult)
{
    BL_PROFILE("SyncRegister::FineAdd()");

    Sync_resid_fine.mult(mult);

    const Box& crse_node_domain = amrex::surroundingNodes(crse_geom.Domain());

    BoxArray cba = Sync_resid_fine.boxArray();
    cba.coarsen(ratio);

    MultiFab Sync_resid_crse(cba, Sync_resid_fine.DistributionMap(), 1, 0);
    Sync_resid_crse.setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
       FArrayBox cbndfab;
       for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
       { 
          for (MFIter mfi(Sync_resid_fine); mfi.isValid(); ++mfi)
          {
             FArrayBox& finefab = Sync_resid_fine[mfi];
             FArrayBox& crsefab = Sync_resid_crse[mfi];

             auto const& crsefab_a = crsefab.array();
             auto const& finefab_a = finefab.array();

             const Box& finebox  = finefab.box();
             const Box& crsebox  = crsefab.box();

             Box bndbox = crsebox;

             for (int side=0; side<2; ++side)
             {
                if (side == 0) {
                   bndbox.setRange(dir,crsebox.smallEnd(dir),1);
                } else {
                   bndbox.setRange(dir,crsebox.bigEnd(dir),1);
                }

                cbndfab.resize(bndbox, 1);
                auto const& cbndfab_a = cbndfab.array();
                Elixir cbndfab_e = cbndfab.elixir();
                amrex::GpuArray<int,AMREX_SPACEDIM> rratio = {D_DECL(ratio[0],ratio[1],ratio[2])};

                // TODO: this doesn't work on GPU. Will probably need to redo this completely ...
                //amrex::launch(bndbox, [finebox,dir,rratio,finefab_a,cbndfab_a]
                //AMREX_GPU_DEVICE (Box const& tbx)
                //{
                   srcrsereg_k(bndbox,finebox,dir,rratio,finefab_a, cbndfab_a);
                //});

                for (int n = 0; n < AMREX_SPACEDIM; ++n)
                {
                   if (!crse_geom.isPeriodic(n))
                   {
                       //
                       // Now points on the physical bndry must be doubled
                       // for any boundary but outflow or periodic
                       //
                       Box domlo(crse_node_domain);
                       domlo.setRange(n,crse_node_domain.smallEnd(n),1);
                       domlo &= bndbox;
                       if (domlo.ok()) {
                          amrex::ParallelFor(domlo, [cbndfab_a]
                          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                          {
                             cbndfab_a(i,j,k) *= 2.0;
                          });
                       }

                      Box domhi(crse_node_domain);
                      domhi.setRange(n,crse_node_domain.bigEnd(n),1);
                      domhi &= bndbox;
                      if (domhi.ok()) {
                          amrex::ParallelFor(domhi, [cbndfab_a]
                          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                          {
                             cbndfab_a(i,j,k) *= 2.0;
                          });
                      }
                   }
                }

                //crsefab += cbndfab;
                const auto& ovlp = crsefab.box() & cbndfab.box();
                if (ovlp.ok())
                {
                  amrex::ParallelFor(ovlp, [crsefab_a,cbndfab_a]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                  {
                     crsefab_a(i,j,k) += cbndfab_a(i,j,k);
                  });
                }
            }
         }
      }
   }

   for (OrientationIter face; face; ++face)
   {
      bndry[face()].plusFrom(Sync_resid_crse,0,0,0,1,crse_geom.periodicity());
   }
}

