//
// $Id: SyncRegister.cpp,v 1.24 1998-05-20 16:10:56 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstdlib>
#include <cstring>
#include <cstdio>
#else
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#endif

#include <BC_TYPES.H>
#include <SyncRegister.H>
#include <Tracer.H>

#include <NAVIERSTOKES_F.H>
#ifndef FORT_HGC2N
#include <PROJECTION_F.H>
#endif

#include <SYNCREG_F.H>

SyncRegister::SyncRegister ()
{
    fine_level = -1;
    ratio = IntVect::TheUnitVector(); ratio.scale(-1);
}

SyncRegister::SyncRegister (const BoxArray& fine_boxes,
                            const IntVect&  ref_ratio,
                            int             fine_lev)
{
    ratio = IntVect::TheUnitVector();
    ratio.scale(-1);
    define(fine_boxes,ref_ratio,fine_lev);
}

void
SyncRegister::define (const BoxArray& fine_boxes,
                      const IntVect&  ref_ratio,
                      int             fine_lev)
{
    TRACER("SyncRegister::define");

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
        assert(ratio[dir] == -1);
    assert(fine_boxes.isDisjoint());
    assert(!grids.ready());

    ratio      = ref_ratio;
    fine_level = fine_lev;
    grids.define(fine_boxes);
    grids.coarsen(ratio);

    const int ngrds = grids.length();

    for (OrientationIter face; face; ++face)
    {
        bndry[face()].resize(ngrds);
        bndry[face()].DefineGrids(grids);
        bndry[face()].DefineDistributionMap(grids);

        bndry_mask[face()].resize(ngrds);
        bndry_mask[face()].DefineGrids(grids);
        bndry_mask[face()].DefineDistributionMap(grids);
    }
    //
    // Construct disjoint "face" fabs that are node centered in all index
    // directions.  The "grow" function at the bottom of this loop shrinks
    // the domain in the index direction just allocated so that the fabs
    // in the other directions will not overlap.
    //
    const int myproc = ParallelDescriptor::MyProc();

    for (int k = 0; k < ngrds; k++)
    {
        Box ndbox(::surroundingNodes(grids[k]));

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            const int* blo = ndbox.loVect();
            Box nd_lo(ndbox);
            nd_lo.setRange(dir,blo[dir],1);
            FabSet& lo = bndry[Orientation(dir,Orientation::low)];
            lo.setBox(k, nd_lo);
            if (lo.DistributionMap()[k] == myproc)
            {
                //
                // Local.
                //
                assert(!lo.defined(k));
                lo.clear(k);
                lo.setFab(k,new FArrayBox(nd_lo,1));
            }
            FabSet& lo_mask = bndry_mask[Orientation(dir,Orientation::low)];
            lo_mask.setBox(k, nd_lo);
            if (lo_mask.DistributionMap()[k] == myproc)
            {
                //
                // Local.
                //
                assert(!lo_mask.defined(k));
                lo_mask.clear(k);
                lo_mask.setFab(k,new FArrayBox(nd_lo,1));
            }
            const int* bhi = ndbox.hiVect();
            Box nd_hi(ndbox);
            nd_hi.setRange(dir,bhi[dir],1);
            FabSet& hi = bndry[Orientation(dir,Orientation::high)];
            hi.setBox(k, nd_hi);
            if (hi.DistributionMap()[k] == myproc)
            {
                //
                // Local.
                //
                assert(!hi.defined(k));
                hi.clear(k);
                hi.setFab(k,new FArrayBox(nd_hi,1));
            }
            FabSet& hi_mask = bndry_mask[Orientation(dir,Orientation::high)];
            hi_mask.setBox(k, nd_hi);
            if (hi_mask.DistributionMap()[k] == myproc)
            {
                //
                // Local.
                //
                assert(!hi_mask.defined(k));
                hi_mask.clear(k);
                hi_mask.setFab(k,new FArrayBox(nd_hi,1));
            }
            assert(ndbox.shortside() > 0);
        }
    }
}

SyncRegister::~SyncRegister () {}

Real
SyncRegister::sum ()
{
    //
    // Sum values in all registers.
    // NOTE: overlap counted twice.
    //
    Box bb(grids[0]);
    for (int k = 1; k < grids.length(); k++)
        bb.minBox(grids[k]);
    bb.surroundingNodes();
    FArrayBox bfab(bb,1);
    bfab.setVal(0.0);

    for (OrientationIter face; face; ++face)
    {
        for (FabSetIterator fsi(bndry[face()]); fsi.isValid(); ++fsi)
        {
            bfab.copy(fsi());
        }
    }
    Real tempSum = bfab.sum(0,1);
    ParallelDescriptor::ReduceRealSum(tempSum);
    return tempSum;
}

void
SyncRegister::increment (const FArrayBox& src)
{
    for (OrientationIter face; face; ++face)
    {
        for (FabSetIterator fsi(bndry[face()]); fsi.isValid(); ++fsi)
        {
            fsi().plus(src);
        }
    }
}

void
SyncRegister::InitRHS (MultiFab&       rhs,
                       const Geometry& geom,
                       const BCRec*    phys_bc)
{
    cerr << "SyncRegister::InitRHS() not implemented in parallel.\n";
    if (ParallelDescriptor::NProcs() > 1)
        ParallelDescriptor::Abort("Bye");

    rhs.setVal(0.0);
    const BoxArray& rhs_boxes = rhs.boxArray();
    int nrhs                  = rhs_boxes.length();
    int nreg                  = grids.length();
    const Box& cell_domain    = geom.Domain();
    const Box& domain         = ::surroundingNodes(cell_domain);

    Array<IntVect> pshifts(27);
    //
    // FILL RHS FROM BNDRY REGISTERS
    //
    // If periodic, copy the values from sync registers onto the nodes
    // of the rhs which are not covered by sync registers through periodic
    // shifts.
    //
    if (geom.isAnyPeriodic())
    {
        for (int k = 0; k < nrhs; k++)
        {
            for (int j = 0; j < nreg; j++)
            {
                for (OrientationIter face; face; ++face)
                {
                    geom.periodicShift(domain,bndry[face()][j].box(),pshifts);

                    for (int iiv = 0; iiv < pshifts.length(); iiv++)
                    {
                        bndry[face()][j].shift(pshifts[iiv]);
                        rhs[k].copy(bndry[face()][j]);
                        bndry[face()][j].shift(-pshifts[iiv]);
                    }
                }
            }
        }
    }
    //
    // This will overwrite the above-set values on all nodes covered
    // by a sync register.
    //
    for (int k = 0; k < nrhs; k++)
    {
        for (int j = 0; j < nreg; j++)
        {
            for (OrientationIter face; face; ++face)
            {
                rhs[k].copy(bndry[face()][j]);
            }
        }
    }
    const int* dlo     = domain.loVect();
    const int* dhi     = domain.hiVect();
    const int* phys_lo = phys_bc->lo();
    const int* phys_hi = phys_bc->hi();

    FArrayBox tmp_rhs;

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        Box domlo(domain), domhi(domain);
        domlo.setRange(dir,dlo[dir],1);
        domhi.setRange(dir,dhi[dir],1);
        //
        // Any RHS point on the physical bndry must be multiplied
        // by two (only for ref-wall and inflow) and set to zero at outflow.
        //
        if (!geom.isPeriodic(dir))
        {
            for (int k = 0; k < nrhs; k++)
            {
                Box blo = rhs_boxes[k] & domlo;
                Box bhi = rhs_boxes[k] & domhi;

                if (blo.ok())
                    rhs[k].mult(2.0,blo,0,1);
                if (bhi.ok())
                    rhs[k].mult(2.0,bhi,0,1);
                if (blo.ok() && phys_lo[dir] == Outflow) 
                    rhs[k].setVal(0.0,blo,0,1);
                if (bhi.ok() && phys_hi[dir] == Outflow) 
                    rhs[k].setVal(0.0,bhi,0,1);
            }
        } 
    }
    //
    // SET UP BNDRY_MASK
    //
    for (int j = 0; j < nreg; j++)
    {
        for (OrientationIter face; face; ++face)
        {
            bndry_mask[face()][j].setVal(0);
        }
    }

    const int ngrds = grids.length();

    for (int j = 0; j < nreg; j++)
    {
        for (OrientationIter face; face; ++face)
        {
            FArrayBox& mask = bndry_mask[face()][j];

            Box mask_cells = ::enclosedCells(::grow(mask.box(),1));
            FArrayBox cellMask(mask_cells,1);
            cellMask.setVal(0);

            for (int n = 0; n < ngrds; n++)
            {
                Box intersect = mask_cells & grids[n];
                if (intersect.ok())
                    cellMask.setVal(1.0,intersect,0,1);
            }
 
            if (geom.isAnyPeriodic())
            {
                geom.periodicShift(cell_domain,mask_cells,pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    IntVect iv = pshifts[iiv];
                    mask_cells.shift(iv);
                    for (int n = 0; n < ngrds; n++)
                    {
                        Box intersect = mask_cells & grids[n];
                        if (intersect.ok())
                        {
                            intersect.shift(-iv);
                            cellMask.setVal(1.0,intersect,0,1);
                        }
                    }
                    mask_cells.shift(-iv);
                }
            }
            REAL* mask_dat = mask.dataPtr();
            const int* mlo = mask.loVect(); 
            const int* mhi = mask.hiVect();
            REAL* cell_dat = cellMask.dataPtr();
            const int* clo = cellMask.loVect(); 
            const int* chi = cellMask.hiVect();
        
            FORT_MAKEMASK(mask_dat,ARLIM(mlo),ARLIM(mhi),
                          cell_dat,ARLIM(clo),ARLIM(chi));
        }
    }

    Box node_domain = ::surroundingNodes(domain);
    const int* ndlo = node_domain.loVect();
    const int* ndhi = node_domain.hiVect();
    //
    // Here double the cell contributions if at a non-periodic physical bdry.
    //
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        if (!geom.isPeriodic(dir))
        {
            Box domlo(node_domain);
            Box domhi(node_domain);
            domlo.setRange(dir,ndlo[dir],1);
            domhi.setRange(dir,ndhi[dir],1);
            for (int j = 0; j < nreg; j++)
            {
                for (OrientationIter face; face; ++face)
                {
                    Box bndry_box_lo(bndry_mask[face()][j].box());
                    Box bndry_box_hi(bndry_mask[face()][j].box());
                    bndry_box_lo &= domlo;
                    bndry_box_hi &= domhi;
                    if (bndry_box_lo.ok()) 
                        bndry_mask[face()][j].mult(2.0,bndry_box_lo,0,1);
                    if (bndry_box_hi.ok()) 
                        bndry_mask[face()][j].mult(2.0,bndry_box_hi,0,1);
                }
            }
        }
    }
    //
    // Here convert from sum of cell contributions to 0. or 1.
    //
    for (int j = 0; j < nreg; j++)
    {
        for (OrientationIter face; face; ++face)
        {
            FArrayBox& mask = bndry_mask[face()][j];
            REAL* mask_dat  = mask.dataPtr();
            const int* mlo  = mask.loVect(); 
            const int* mhi  = mask.hiVect();

            FORT_CONVERTMASK(mask_dat,ARLIM(mlo),ARLIM(mhi));
        }
    }
    //
    // MULTIPLY RHS BY BNDRY_MASK 
    //
    for (int k = 0; k < nrhs; k++)
    {
        for (int j = 0; j < nreg; j++)
        {
            for (OrientationIter face; face; ++face)
            {
                rhs[k].mult(bndry_mask[face()][j]);
            }
        }
    }
}

void
SyncRegister::CrseDVInit (const MultiFab& U, 
                          const Geometry& geom, 
                          int             is_rz,
                          int**           crse_bc,
                          Real            mult)
{
    //
    // Zero all registers.
    //
    setVal(0);

    const int n_ghost = 1;

    MultiFab U_local(U.boxArray(),BL_SPACEDIM,n_ghost);
    //
    // Nullify ghost cells.
    //
    U_local.setBndry(0);
    //
    // First fill all the coarse cells, including ghost cells on periodic
    // and ext_dir edges, before worrying about zeroing out the ones under
    // fine grids.
    //
    for (MultiFabIterator mfi(U_local); mfi.isValid(false); ++mfi)
    {
        //
        // U and U_local have same BoxArray and hence same DistributionMapping.
        //
        DependentMultiFabIterator dmfi_U(mfi, U);

        assert(dmfi_U.validbox() == mfi.validbox());

        mfi().copy(dmfi_U(),dmfi_U.validbox(),0,dmfi_U.validbox(),0,BL_SPACEDIM);

        int* bc = crse_bc[mfi.index()];

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            int bc_index = 2*BL_SPACEDIM*dir + dir;
            //
            // Fill ghost cells outside of domain.
            // 
            if (dmfi_U.validbox().smallEnd(dir) == geom.Domain().smallEnd(dir))
            {
                Box sidelo(dmfi_U.validbox());
                sidelo.setRange(dir,sidelo.smallEnd(dir)-1,1);
                if (bc[bc_index] == EXT_DIR)
                    mfi().copy(dmfi_U(),sidelo,dir,sidelo,dir,1);
            }
            if (dmfi_U.validbox().bigEnd(dir) == geom.Domain().bigEnd(dir))
            {
                Box sidehi(dmfi_U.validbox());
                sidehi.setRange(dir,sidehi.bigEnd(dir)+1,1);
                if (bc[bc_index+BL_SPACEDIM] == EXT_DIR)
                    mfi().copy(dmfi_U(),sidehi,dir,sidehi,dir,1);
            }
        }
    }
    //
    // Enforce periodicity of U_local using extended valid boxes.
    //
    geom.FillPeriodicBoundary(U_local, true);
    //
    // Now do all the zeroing-out associated with the fine grids, on the
    // interior and on ghost cells on periodic and ext_dir edges.
    //
    Array<IntVect> pshifts(27);

    for (MultiFabIterator mfi(U_local); mfi.isValid(false); ++mfi)
    {
        assert(mfi.validbox() == U.boxArray()[mfi.index()]);

        int* bc = crse_bc[mfi.index()];

        for (int fine = 0; fine < grids.length(); fine++)
        {
            Box subbox = mfi.validbox() & grids[fine];

            if (subbox.ok())
            {
                mfi().setVal(0.0,subbox,0,BL_SPACEDIM);

                for (int dir = 0; dir < BL_SPACEDIM; dir++)
                {
                    int bc_index = 2*BL_SPACEDIM*dir + dir;
                    if (bc[bc_index] == EXT_DIR &&
                        grids[fine].smallEnd(dir) == mfi.validbox().smallEnd(dir))
                    {
                        Box finesidelo(subbox);
                        finesidelo.setRange(dir,finesidelo.smallEnd(dir)-1,1);
                        mfi().setVal(0.0,finesidelo,0,BL_SPACEDIM);
                    }
                    if (bc[bc_index+BL_SPACEDIM] == EXT_DIR &&
                        grids[fine].bigEnd(dir) == mfi.validbox().bigEnd(dir))
                    {
                        Box finesidehi(subbox);
                        finesidehi.setRange(dir,finesidehi.bigEnd(dir)+1,1);
                        mfi().setVal(0.0,finesidehi,0,BL_SPACEDIM);
                    }
                }
            }
            //
            // Now zero out under periodic translations of fine grids.
            //
            if (geom.isAnyPeriodic())
            {
                Box domain_plus(geom.Domain());

                for (int dir = 0; dir < BL_SPACEDIM; dir++) 
                    if (geom.isPeriodic(dir))
                        domain_plus.grow(dir,1);

                geom.periodicShift(domain_plus, grids[fine], pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    Box fine_shifted(grids[fine]);
                    fine_shifted.shift(pshifts[iiv]);
                    fine_shifted &= mfi().box();

                    if (fine_shifted.ok())
                    {
                        mfi().setVal(0.0,fine_shifted,0,BL_SPACEDIM);

                        for (int dir = 0; dir < BL_SPACEDIM; dir++)
                        {
                            int bc_index = 2*BL_SPACEDIM*dir + dir;

                            if (bc[bc_index] == EXT_DIR &&
                                fine_shifted.smallEnd(dir) == mfi.validbox().smallEnd(dir))
                            {
                                Box finesidelo(fine_shifted);
                                finesidelo.setRange(dir,fine_shifted.smallEnd(dir)-1,1);
                                mfi().setVal(0.0,finesidelo,0,BL_SPACEDIM);
                            }

                            if (bc[bc_index+BL_SPACEDIM] == EXT_DIR &&
                                grids[fine].bigEnd(dir) == mfi.validbox().bigEnd(dir))
                            {
                                Box finesidehi(fine_shifted);
                                finesidehi.setRange(dir,fine_shifted.bigEnd(dir)+1,1);
                                mfi().setVal(0.0,finesidehi,0,BL_SPACEDIM);
                            }
                        }
                    }
                }
            }
        }
    }
    //
    // Now compute node-centered divergence.
    //
    FArrayBox divu;

    for (MultiFabIterator mfi(U_local); mfi.isValid(false); ++mfi)
    {
        const int* ulo  = mfi().loVect();
        const int* uhi  = mfi().hiVect();
        Box ndbox       = ::surroundingNodes(mfi.validbox());
        const int* ndlo = ndbox.loVect();
        const int* ndhi = ndbox.hiVect();

        divu.resize(ndbox,1);

        FORT_SRDIVU(mfi().dataPtr(),ARLIM(ulo),ARLIM(uhi),
                    divu.dataPtr(),ARLIM(ndlo),ARLIM(ndhi),
                    ndlo,ndhi,geom.CellSize(),&mult,&is_rz);

        increment(divu);
    }
}

void
SyncRegister::FineDVAdd (const MultiFab& U, 
                         const Real*     dx_fine, 
                         const Geometry& crse_geom, 
                         int             is_rz,
                         int**           fine_bc,
                         Real            mult)
{
    const Box& crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    FArrayBox ufab;
    FArrayBox cfablo, cfabhi, ffablo, ffabhi;

    Array<IntVect> pshifts(27);

    for (ConstMultiFabIterator mfi(U); mfi.isValid(false); ++mfi)
    {
        Box ubox = ::grow(mfi.validbox(),1);

        ufab.resize(ubox,BL_SPACEDIM);
        ufab.setComplement(0,mfi.validbox(),0,BL_SPACEDIM);
        ufab.copy(mfi(),mfi.validbox(),0,mfi.validbox(),0,BL_SPACEDIM);

        const int* ulo = ubox.loVect();
        const int* uhi = ubox.hiVect();
        int* bc        = fine_bc[mfi.index()];

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            int bc_index = 2*BL_SPACEDIM*dir + dir;
            if (bc[bc_index] == EXT_DIR)
            {
                Box sidelo(mfi.validbox());
                sidelo.growLo(dir,1);
                const int* dlo = sidelo.loVect();
                sidelo.setRange(dir,dlo[dir],1);
                ufab.copy(mfi(),sidelo,dir,sidelo,dir,1);
            }
            if (bc[bc_index+BL_SPACEDIM] == EXT_DIR)
            {
                Box sidehi(mfi.validbox());
                sidehi.growHi(dir,1);
                const int* dhi = sidehi.hiVect();
                sidehi.setRange(dir,dhi[dir],1);
                ufab.copy(mfi(),sidehi,dir,sidehi,dir,1);
            }
        }
        //
        // Now compute node centered surrounding box.
        //
        Box ndbox       = ::surroundingNodes(mfi.validbox());
        const int* ndlo = ndbox.loVect();
        const int* ndhi = ndbox.hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            //
            // Determine region of interest, and size of tmp fabs.
            //
            Box tboxlo(ndbox), tboxhi(ndbox);
            tboxlo.setRange(dir,ndlo[dir],1);
            tboxhi.setRange(dir,ndhi[dir],1);
            Box cboxlo(tboxlo), cboxhi(tboxhi);
            Box reglo(tboxlo), reghi(tboxhi);
            for (int idir = 0; idir < BL_SPACEDIM; idir++)
            {
                if (idir < dir)
                {
                    //
                    // Previous direction have included these
                    // points in the stencil already, shrink.
                    //
                    reglo.grow(idir,-1);
                    reghi.grow(idir,-1);
                }
                if (idir != dir)
                {
                    //
                    // Need additional room for stencil calculation.
                    //
                    tboxlo.grow(idir,ratio[idir]-1);
                    tboxhi.grow(idir,ratio[idir]-1);
                }
            }
            //
            // Define fine grid tmp fabs.
            //
            const int* flo_lo = tboxlo.loVect();
            const int* flo_hi = tboxlo.hiVect();
            const int* fhi_lo = tboxhi.loVect();
            const int* fhi_hi = tboxhi.hiVect();
            ffablo.resize(tboxlo,1);
            ffabhi.resize(tboxhi,1);
            ffablo.setVal(0.0);
            ffabhi.setVal(0.0);
            //
            // Define coarsened tmp fabs.
            //
            cboxlo.coarsen(ratio);
            cboxhi.coarsen(ratio);
            cfablo.resize(cboxlo,1);
            cfabhi.resize(cboxhi,1);
            cfablo.setVal(0.0);
            cfabhi.setVal(0.0);
            //
            // Compute divu on fine grid edges in regions defined
            // by reglo and reghi.  Fabs are set to zero outside region.
            //
            FORT_SRDIVU(ufab.dataPtr(),ARLIM(ulo),ARLIM(uhi),
                        ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                        reglo.loVect(),reglo.hiVect(),
                        dx_fine,&mult,&is_rz);
            FORT_SRDIVU(ufab.dataPtr(),ARLIM(ulo),ARLIM(uhi),
                        ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                        reghi.loVect(),reghi.hiVect(),
                        dx_fine,&mult,&is_rz);
            //
            // Coarsen edge value.
            //
            const int* clo = cboxlo.loVect();
            const int* chi = cboxlo.hiVect();
            FORT_SRCRSEREG(ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                           cfablo.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            clo = cboxhi.loVect();
            chi = cboxhi.hiVect();
            FORT_SRCRSEREG(ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                           cfabhi.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            //
            // Intersect and add to registers.
            //
            increment(cfablo);
            increment(cfabhi);

            if (crse_geom.isAnyPeriodic())
            {
                crse_geom.periodicShift(crse_node_domain, cboxlo, pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    cfablo.shift(pshifts[iiv]);
                    increment(cfablo);
                    cfablo.shift(-pshifts[iiv]);
                }
                crse_geom.periodicShift(crse_node_domain, cboxhi, pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    cfabhi.shift(pshifts[iiv]);
                    increment(cfabhi);
                    cfabhi.shift(-pshifts[iiv]);
                }
            }
        }
    }
}

static
void
SetCenter (int             is_rz,
           Array<Real>&    rcen,
           const Geometry& geom,
           const Box&      dsdtbox)
{
    if (is_rz)
    {
        geom.GetCellLoc(rcen, dsdtbox, 0);
    }
    else
    {
        for (int i = 0; i < rcen.length(); i++)
            rcen[i] = 1.0;
    }
}

void
SyncRegister::CrseDsdtAdd (const MultiFab& dsdt,
                           const Geometry& geom,
                           int             is_rz,
                           int**           crse_bc, 
                           int             lowfix,
                           int             hifix,
                           Real            mult)
{
    FArrayBox dsdtfab, divu;

    for (ConstMultiFabIterator mfi(dsdt); mfi.isValid(false); ++mfi)
    {
        dsdtfab.resize(::grow(mfi.validbox(),1),1);
        dsdtfab.setComplement(0,mfi.validbox(),0,1);
        dsdtfab.copy(mfi(),mfi.validbox(),0,mfi.validbox(),0,1);

        int* bc = crse_bc[mfi.index()];

        for (int fine = 0; fine < grids.length(); fine++)
        {
            Box subbox = dsdtfab.box() & grids[fine];

            if (subbox.ok())
            {
                dsdtfab.setVal(0.0,subbox,0,1);

                for (int dir = 0; dir < BL_SPACEDIM; dir++)
                {
                    int bc_index = 2*BL_SPACEDIM*dir + dir;

                    if (bc[bc_index] == EXT_DIR &&
                        grids[fine].loVect()[dir] == mfi.validbox().loVect()[dir])
                    {
                        Box finesidelo(subbox);
                        finesidelo.growLo(dir,1);
                        finesidelo.setRange(dir,finesidelo.loVect()[dir],1);
                        dsdtfab.setVal(0.0,finesidelo,0,1);
                    }
                    if (bc[bc_index+BL_SPACEDIM] == EXT_DIR &&
                        grids[fine].hiVect()[dir] == mfi.validbox().hiVect()[dir])
                    {
                        Box finesidehi(subbox);
                        finesidehi.growHi(dir,1);
                        finesidehi.setRange(dir,finesidehi.hiVect()[dir],1);
                        dsdtfab.setVal(0.0,finesidehi,0,1);
                    }
                }
            }
        }
        //
        // Average dsdt to nodes.
        //
        divu.resize(::surroundingNodes(mfi.validbox()),1);

        const int* ndlo   = divu.box().loVect();
        const int* ndhi   = divu.box().hiVect();
        const int* dsdtlo = dsdtfab.box().loVect();
        const int* dsdthi = dsdtfab.box().hiVect();
        const int* rlo    = dsdtfab.box().loVect();
        const int* rhi    = dsdtfab.box().hiVect();
        const int* domlo  = geom.Domain().loVect();
        const int* domhi  = geom.Domain().hiVect();

        Array<Real> rcen(dsdtfab.box().length(0));

        SetCenter(is_rz, rcen, geom, dsdtfab.box());

#if (BL_SPACEDIM==2)
        int nghost         = 0;
        Real hx            = geom.CellSize()[0];
        int extrap_edges   = 0;
        int extrap_corners = 0;
        FORT_HGC2N(&nghost, ARLIM(dsdtlo), ARLIM(dsdthi), 
                   dsdtfab.dataPtr(),
                   rcen.dataPtr(), 
                   ARLIM(ndlo), ARLIM(ndhi), divu.dataPtr(), 
                   domlo, domhi, lowfix, hifix, &hx,
                   &extrap_edges, &extrap_corners, &is_rz);
#elif (BL_SPACEDIM==3)
        divu.setVal(0.0);
#endif
        divu.negate();
        divu.mult(mult);
        increment(divu);
    }
}

void
SyncRegister::FineDsdtAdd (const MultiFab& dsdt,
                           const Geometry& geom,
                           int             is_rz,
                           int**           fine_bc, 
                           int             lowfix,
                           int             hifix,
                           Real            mult)
{
    FArrayBox dsdtfab;
    FArrayBox cfablo, cfabhi, ffablo, ffabhi;

    for (ConstMultiFabIterator mfi(dsdt); mfi.isValid(false); ++mfi)
    {
        dsdtfab.resize(::grow(mfi.validbox(),1),1);
        dsdtfab.setComplement(0,mfi.validbox(),0,1);
        dsdtfab.copy(mfi(),mfi.validbox(),0,mfi.validbox(),0,1);

        const int* dsdtlo = dsdtfab.box().loVect();
        const int* dsdthi = dsdtfab.box().hiVect();
        //
        // Now compute node centered surrounding box.
        //
        Box ndbox       = ::surroundingNodes(mfi.validbox());
        const int* ndlo = ndbox.loVect();
        const int* ndhi = ndbox.hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            //
            // Determine region of interest, and size of tmp fabs.
            //
            Box tboxlo(ndbox), tboxhi(ndbox);
            tboxlo.setRange(dir,ndlo[dir],1);
            tboxhi.setRange(dir,ndhi[dir],1);
            Box cboxlo(tboxlo), cboxhi(tboxhi);
            Box reglo(tboxlo), reghi(tboxhi);
            for (int idir = 0; idir < BL_SPACEDIM; idir++)
            {
                if (idir < dir)
                {
                    //
                    // Previous direction have included these
                    // points in the stencil already, shrink.
                    //
                    reglo.grow(idir,-1);
                    reghi.grow(idir,-1);
                }
                if (idir != dir)
                {
                    //
                    // Need additional room for stencil calculation.
                    //
                    tboxlo.grow(idir,ratio[idir]-1);
                    tboxhi.grow(idir,ratio[idir]-1);
                }
            }
            //
            // Define fine grid tmp fabs.
            //
            const int* flo_lo = tboxlo.loVect();
            const int* flo_hi = tboxlo.hiVect();
            const int* fhi_lo = tboxhi.loVect();
            const int* fhi_hi = tboxhi.hiVect();
            ffablo.resize(tboxlo,1);
            ffabhi.resize(tboxhi,1);
            ffablo.setVal(0.0);
            ffabhi.setVal(0.0);
            //
            // Define coarsened tmp fabs.
            //
            cboxlo.coarsen(ratio);
            cboxhi.coarsen(ratio);
            cfablo.resize(cboxlo,1);
            cfabhi.resize(cboxhi,1);
            cfablo.setVal(0.0);
            cfabhi.setVal(0.0);
            //
            // Average dsdt to nodes on fine grid edges in regions defined
            // by reglo and reghi.  Fabs are set to zero outside region.
            //
            const int* dsdtlo = dsdtfab.box().loVect();
            const int* dsdthi = dsdtfab.box().hiVect();
            const int* rlo    = dsdtfab.box().loVect();
            const int* rhi    = dsdtfab.box().hiVect();
            const int* domlo  = geom.Domain().loVect();
            const int* domhi  = geom.Domain().hiVect();

            Array<Real> rcen(dsdtfab.box().length(0));

            SetCenter(is_rz, rcen, geom, dsdtfab.box());

            FArrayBox ffablo_tmp(reglo,1);

#if (BL_SPACEDIM==2)
            int nghost         = 0;
            int hi_fix         = 0;
            Real hx            = geom.CellSize()[0];
            int extrap_edges   = 0;
            int extrap_corners = 0;
            FORT_HGC2N(&nghost, ARLIM(dsdtlo), ARLIM(dsdthi), 
                       dsdtfab.dataPtr(),
                       rcen.dataPtr(), 
                       ARLIM(reglo.loVect()), ARLIM(reglo.hiVect()), 
                       ffablo_tmp.dataPtr(),
                       domlo, domhi, lowfix, hi_fix, &hx,
                       &extrap_edges, &extrap_corners, &is_rz);
#elif (BL_SPACEDIM==3)
            ffablo_tmp.setVal(0.0);
#endif
            ffablo_tmp.negate();
            ffablo_tmp.mult(mult);
            ffablo.copy(ffablo_tmp);
            FArrayBox ffabhi_tmp(reghi,1);

#if (BL_SPACEDIM==2)
            int low_fix = 0;
            FORT_HGC2N(&nghost, ARLIM(dsdtlo), ARLIM(dsdthi), 
                       dsdtfab.dataPtr(),
                       rcen.dataPtr(), 
                       ARLIM(reghi.loVect()), ARLIM(reghi.hiVect()), 
                       ffabhi_tmp.dataPtr(), 
                       domlo, domhi, low_fix, hifix, &hx,
                       &extrap_edges, &extrap_corners, &is_rz);
#elif (BL_SPACEDIM==3)
            ffabhi_tmp.setVal(0.0);
#endif
            ffabhi_tmp.negate();
            ffabhi_tmp.mult(mult);
            ffabhi.copy(ffabhi_tmp);
            //
            // Coarsen edge value.
            //
            const int* clo = cboxlo.loVect();
            const int* chi = cboxlo.hiVect();
            FORT_SRCRSEREG(ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                           cfablo.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            clo = cboxhi.loVect();
            chi = cboxhi.hiVect();
            FORT_SRCRSEREG(ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                           cfabhi.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            //
            // Intersect and add to registers.
            //
            increment(cfablo);
            increment(cfabhi);
        }
    }
}

void
SyncRegister::CompDVAdd (const MultiFab& U, 
                         const BoxArray& Pgrids,
                         const Real*     dx_fine, 
                         const Geometry& fine_geom, 
                         const Geometry& crse_geom, 
                         int             is_rz,
                         int**           fine_bc,
                         Real            mult)
{
    cerr << "SyncRegister::CompDVAdd() not implemented in parallel.\n";
    if (ParallelDescriptor::NProcs() > 1)
        ParallelDescriptor::Abort("Bye");

    const BoxArray& U_boxes     = U.boxArray();
    const int ngrds             = U_boxes.length();
    const Box& crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    Array<IntVect> pshifts(27);

    FArrayBox ufab;
    FArrayBox cfablo, cfabhi, ffablo, ffabhi;

    for (int k = 0; k < ngrds; k++)
    {
        Box ubox = ::grow(U_boxes[k],1);

        ufab.resize(ubox,BL_SPACEDIM);
        ufab.setVal(0.0);
        ufab.copy(U[k],U_boxes[k],0,U_boxes[k],0,BL_SPACEDIM);

        const int* ulo = ubox.loVect();
        const int* uhi = ubox.hiVect();
        int * bc       = fine_bc[k];

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            int bc_index = 2*BL_SPACEDIM*dir + dir;
            if (bc[bc_index] == EXT_DIR)
            {
                Box sidelo(U_boxes[k]);
                sidelo.growLo(dir,1);
                const int* dlo = sidelo.loVect();
                sidelo.setRange(dir,dlo[dir],1);
                ufab.copy(U[k],sidelo,dir,sidelo,dir,1);
            }
            if (bc[bc_index+BL_SPACEDIM] == EXT_DIR)
            {
                Box sidehi(U_boxes[k]);
                sidehi.growHi(dir,1);
                const int* dhi = sidehi.hiVect();
                sidehi.setRange(dir,dhi[dir],1);
                ufab.copy(U[k],sidehi,dir,sidehi,dir,1);
            }
        }
        //
        // Now compute node centered surrounding box.
        //
        Box ndbox       = ::surroundingNodes(U_boxes[k]);
        const int* ndlo = ndbox.loVect();
        const int* ndhi = ndbox.hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            //
            // Determine region of interest, and size of tmp fabs.
            //
            Box tboxlo(ndbox), tboxhi(ndbox);
            tboxlo.setRange(dir,ndlo[dir],1);
            tboxhi.setRange(dir,ndhi[dir],1);
            Box cboxlo(tboxlo), cboxhi(tboxhi);
            Box reglo(tboxlo), reghi(tboxhi);
            for (int idir = 0; idir < BL_SPACEDIM; idir++)
            {
                if (idir < dir)
                {
                    //
                    // Previous direction have included these
                    // points in the stencil already, shrink.
                    //
                    reglo.grow(idir,-1);
                    reghi.grow(idir,-1);
                }
                if (idir != dir)
                {
                    //
                    // Need additional room for stencil calculation.
                    //
                    tboxlo.grow(idir,ratio[idir]-1);
                    tboxhi.grow(idir,ratio[idir]-1);
                }
            }
            //
            // Define fine grid tmp fabs.
            //
            const int* flo_lo = tboxlo.loVect();
            const int* flo_hi = tboxlo.hiVect();
            const int* fhi_lo = tboxhi.loVect();
            const int* fhi_hi = tboxhi.hiVect();
            ffablo.resize(tboxlo,1);
            ffabhi.resize(tboxhi,1);
            ffablo.setVal(0.0);
            ffabhi.setVal(0.0);
            //
            // Compute divu on fine grid edges in regions defined
            // by reglo and reghi.  Fabs are set to zero outside region.
            //
            FORT_SRDIVU(ufab.dataPtr(),ARLIM(ulo),ARLIM(uhi),
                        ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                        reglo.loVect(),reglo.hiVect(),dx_fine,&mult,&is_rz);
            FORT_SRDIVU(ufab.dataPtr(),ARLIM(ulo),ARLIM(uhi),
                        ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                        reghi.loVect(),reghi.hiVect(),dx_fine,&mult,&is_rz);

            int n_comp   = 1;
            int set_comp = 0;
            for (int i = 0; i < Pgrids.length(); i++)
            {
                Box overlap_lo = reglo & Pgrids[i];
                if (overlap_lo.ok())
                    ffablo.setVal(0.,overlap_lo,set_comp,n_comp);

                fine_geom.periodicShift(reglo, Pgrids[i], pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    Box overlap_lo_per(Pgrids[i]);
                    overlap_lo_per.shift(pshifts[iiv]);
                    overlap_lo_per &= reglo;
                    ffablo.setVal(0.,overlap_lo_per,set_comp,n_comp);
                }

                Box overlap_hi = reghi & Pgrids[i];
                if (overlap_hi.ok())
                    ffabhi.setVal(0.,overlap_hi,set_comp,n_comp);

                fine_geom.periodicShift(reghi, Pgrids[i], pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    Box overlap_hi_per(Pgrids[i]);
                    overlap_hi_per.shift(pshifts[iiv]);
                    overlap_hi_per &= reghi;
                    ffabhi.setVal(0.,overlap_hi_per,set_comp,n_comp);
                }
            }
            //
            // Define coarsened tmp fabs.
            //
            cboxlo.coarsen(ratio);
            cboxhi.coarsen(ratio);
            cfablo.resize(cboxlo,1);
            cfabhi.resize(cboxhi,1);
            cfablo.setVal(0.0);
            cfabhi.setVal(0.0);
            //
            // Coarsen edge values.
            //
            const int* clo = cboxlo.loVect();
            const int* chi = cboxlo.hiVect();
            FORT_SRCRSEREG(ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                           cfablo.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            clo = cboxhi.loVect();
            chi = cboxhi.hiVect();
            FORT_SRCRSEREG(ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                           cfabhi.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            //
            // Intersect and add to registers.
            //
            increment(cfablo);
            increment(cfabhi);

            if (crse_geom.isAnyPeriodic())
            {
                crse_geom.periodicShift(crse_node_domain, cboxlo, pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    cfablo.shift(pshifts[iiv]);
                    increment(cfablo);
                    cfablo.shift(-pshifts[iiv]);
                }
                crse_geom.periodicShift(crse_node_domain, cboxhi, pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    cfabhi.shift(pshifts[iiv]);
                    increment(cfabhi);
                    cfabhi.shift(-pshifts[iiv]);
                }
            }
        }
    }
}

void
SyncRegister::CrseLPhiAdd (const MultiFab& Phi,
                           const MultiFab& Sigma,
                           const Geometry& geom,
                           int             is_rz,
                           Real            mult)
{
    const int n_ghost = 1;
    //
    // This code assumes Sigma and Phi have same processor distribution.
    //
    assert(Phi.boxArray().length() == Sigma.boxArray().length());

    MultiFab Sig_local(Sigma.boxArray(),1,n_ghost);
    MultiFab Phi_local(Phi.boxArray(),1,n_ghost);
    //
    // Nullify ghost cells.
    //
    Sig_local.setBndry(0);
    Phi_local.setBndry(0);
    //
    // Copy valid region of Phi into Phi_local
    //
    for (MultiFabIterator mfi(Phi_local); mfi.isValid(false); ++mfi)
    {
        mfi().copy(Phi[mfi.index()], mfi.validbox());
    }
    //
    // Copy valid region of Sigma into Sig_local.
    // Also, zero out region covered by fine grid.
    //
    for (MultiFabIterator mfi(Sig_local); mfi.isValid(false); ++mfi)
    {
        mfi().copy(Sigma[mfi.index()], mfi.validbox());

        for (int fine = 0; fine < grids.length(); fine++)
        {
            Box subbox = grids[fine] & mfi.validbox();

            if (subbox.ok())
                mfi().setVal(0.0,subbox,0,1);
        }
    }
    //
    // Enforce periodicity of Sig_local using extended valid boxes.
    //
    geom.FillPeriodicBoundary(Sig_local, 0, 1, true);
    //
    // Enforce periodicity of Phi_local.
    //
    geom.FillPeriodicBoundary(Phi_local, 0, 1);
    //
    // Now compute node centered div(Sigma*grad(PHI)).
    //
    FArrayBox divgp;

    for (MultiFabIterator pmfi(Phi_local); pmfi.isValid(false); ++pmfi)
    {
        DependentMultiFabIterator smfi(pmfi, Sig_local);

        const int* slo  = smfi().loVect();
        const int* shi  = smfi().hiVect();
        const int* p_lo = pmfi().loVect();
        const int* p_hi = pmfi().hiVect();
        divgp.resize(pmfi.validbox(),1);
        const int* glo  = divgp.loVect();
        const int* ghi  = divgp.hiVect();

        FORT_SRDGPHI(pmfi().dataPtr(),ARLIM(p_lo),ARLIM(p_hi),
                     smfi().dataPtr(),ARLIM(slo),ARLIM(shi),
                     divgp.dataPtr(),ARLIM(glo),ARLIM(ghi),
                     glo,ghi,geom.CellSize(),&mult,&is_rz);

        increment(divgp);
    }
}

void
SyncRegister::FineLPhiAdd (const MultiFab& Phi,
                           const MultiFab& Sigma,
                           const Real*     dx_fine,
                           const Geometry& crse_geom, 
                           int             is_rz,
                           Real            mult)
{
    //
    // This code assumes Sigma and Phi have same processor distribution.
    //
    assert(Phi.boxArray().length() == Sigma.boxArray().length());

    const Box& crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    FArrayBox pfab, sfab, cfablo, cfabhi, ffablo, ffabhi;

    Array<IntVect> pshifts(27);

    for (ConstMultiFabIterator pmfi(Phi); pmfi.isValid(false); ++pmfi)
    {
        const Box& ndbox = pmfi.validbox();

        pfab.resize(::grow(ndbox,1),1);
        pfab.setComplement(0.0,pmfi.validbox(),0,1);
        pfab.copy(pmfi(),pmfi.validbox());

        const int* pfab_lo = pfab.loVect();
        const int* pfab_hi = pfab.hiVect();

        ConstDependentMultiFabIterator smfi(pmfi, Sigma);

        sfab.resize(::grow(smfi.validbox(),1),1);
        sfab.setComplement(0.0,smfi.validbox(),0,1);
        sfab.copy(smfi(),smfi.validbox());

        const int* sfab_lo = sfab.loVect();
        const int* sfab_hi = sfab.hiVect();
        const int* ndlo    = ndbox.loVect();
        const int* ndhi    = ndbox.hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            //
            // Determine region of interest, and size of tmp fabs.
            //
            Box tboxlo(ndbox), tboxhi(ndbox);
            tboxlo.setRange(dir,ndlo[dir],1);
            tboxhi.setRange(dir,ndhi[dir],1);
            Box cboxlo(tboxlo), cboxhi(tboxhi);
            Box reglo(tboxlo), reghi(tboxhi);
            for (int idir = 0; idir < BL_SPACEDIM; idir++)
            {
                if (idir < dir)
                {
                    //
                    // Previous direction have included these
                    // points in the stencil already, shrink.
                    //
                    reglo.grow(idir,-1);
                    reghi.grow(idir,-1);
                }
                if (idir != dir)
                {
                    //
                    // Need additional room for stencil calculation.
                    //
                    tboxlo.grow(idir,ratio[idir]-1);
                    tboxhi.grow(idir,ratio[idir]-1);
                }
            }
            //
            // Define fine grid tmp fabs.
            //
            const int* flo_lo = tboxlo.loVect();
            const int* flo_hi = tboxlo.hiVect();
            const int* fhi_lo = tboxhi.loVect();
            const int* fhi_hi = tboxhi.hiVect();
            ffablo.resize(tboxlo,1);
            ffabhi.resize(tboxhi,1);
            ffablo.setVal(0.0);
            ffabhi.setVal(0.0);
            //
            // Define coarsened tmp fabs.
            //
            cboxlo.coarsen(ratio);
            cboxhi.coarsen(ratio);
            cfablo.resize(cboxlo,1);
            cfabhi.resize(cboxhi,1);
            cfablo.setVal(0.0);
            cfabhi.setVal(0.0);
            //
            // Compute divgp on fine grid edges in regions defined
            // by reglo and reghi.  Fabs are set to zero outside region.
            //
            FORT_SRDGPHI(pfab.dataPtr(),ARLIM(pfab_lo),ARLIM(pfab_hi),
                         sfab.dataPtr(),ARLIM(sfab_lo),ARLIM(sfab_hi),
                         ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                         reglo.loVect(),reglo.hiVect(),
                         dx_fine,&mult,&is_rz);
            FORT_SRDGPHI(pfab.dataPtr(),ARLIM(pfab_lo),ARLIM(pfab_hi),
                         sfab.dataPtr(),ARLIM(sfab_lo),ARLIM(sfab_hi),
                         ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                         reghi.loVect(),reghi.hiVect(),
                         dx_fine,&mult,&is_rz);
            //
            // Coarsen edge value.
            //
            const int* clo = cboxlo.loVect();
            const int* chi = cboxlo.hiVect();
            FORT_SRCRSEREG(ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                           cfablo.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            clo = cboxhi.loVect();
            chi = cboxhi.hiVect();
            FORT_SRCRSEREG(ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                           cfabhi.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            //
            // Intersect and add to registers.
            //
            increment(cfablo);
            increment(cfabhi);

            if (crse_geom.isAnyPeriodic())
            {
                crse_geom.periodicShift(crse_node_domain, cboxlo, pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    cfablo.shift(pshifts[iiv]);
                    increment(cfablo);
                    cfablo.shift(-pshifts[iiv]);
                }
                crse_geom.periodicShift(crse_node_domain, cboxhi, pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    cfabhi.shift(pshifts[iiv]);
                    increment(cfabhi);
                    cfabhi.shift(-pshifts[iiv]);
                }
            }
        }
    }
}

void
SyncRegister::CompLPhiAdd (const MultiFab& Phi,
                           const MultiFab& sigma,
                           const BoxArray& Pgrids, 
                           const Real*     dx_fine, 
                           const Geometry& fine_geom, 
                           const Geometry& crse_geom, 
                           int             is_rz,
                           Real            mult)
{
    cerr << "SyncRegister::CompLPhiAdd() not implemented in parallel.\n";
    if (ParallelDescriptor::NProcs() > 1)
        ParallelDescriptor::Abort("Bye");

    const Box& crse_node_domain = ::surroundingNodes(crse_geom.Domain());
    const BoxArray& Phi_boxes   = Phi.boxArray();
    const int ngrds             = Phi_boxes.length();
    const BoxArray& Sig_boxes   = sigma.boxArray();

    FArrayBox pfab, sfab;
    FArrayBox cfablo, cfabhi, ffablo, ffabhi;

    Array<IntVect> pshifts(27);

    for (int k = 0; k < ngrds; k++)
    {
        const Box& ndbox = Phi_boxes[k];
        Box pbox = ::grow(ndbox,1);
        pfab.resize(pbox,1);
        pfab.setVal(0.0);
        pfab.copy(Phi[k],Phi_boxes[k]);
        const int* pfab_lo = pfab.loVect();
        const int* pfab_hi = pfab.hiVect();

        const FArrayBox& sig = sigma[k];
        sfab.resize(::grow(Sig_boxes[k],1),1);
        sfab.setVal(0.0);
        sfab.copy(sig,Sig_boxes[k]);
        const int* sfab_lo = sfab.loVect();
        const int* sfab_hi = sfab.hiVect();

        const int* ndlo = ndbox.loVect();
        const int* ndhi = ndbox.hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            //
            // Determine region of interest, and size of tmp fabs.
            //
            Box tboxlo(ndbox), tboxhi(ndbox);
            tboxlo.setRange(dir,ndlo[dir],1);
            tboxhi.setRange(dir,ndhi[dir],1);
            Box cboxlo(tboxlo), cboxhi(tboxhi);
            Box reglo(tboxlo), reghi(tboxhi);
            for (int idir = 0; idir < BL_SPACEDIM; idir++)
            {
                if (idir < dir)
                {
                    //
                    // Previous direction have included these
                    // points in the stencil already, shrink.
                    //
                    reglo.grow(idir,-1);
                    reghi.grow(idir,-1);
                }
                if (idir != dir)
                {
                    //
                    // Need additional room for stencil calculation.
                    //
                    tboxlo.grow(idir,ratio[idir]-1);
                    tboxhi.grow(idir,ratio[idir]-1);
                }
            }
            //
            // Define fine grid tmp fabs.
            //
            const int* flo_lo = tboxlo.loVect();
            const int* flo_hi = tboxlo.hiVect();
            const int* fhi_lo = tboxhi.loVect();
            const int* fhi_hi = tboxhi.hiVect();
            ffablo.resize(tboxlo,1);
            ffabhi.resize(tboxhi,1);
            ffablo.setVal(0.0);
            ffabhi.setVal(0.0);
            //
            // Define coarsened tmp fabs.
            //
            cboxlo.coarsen(ratio);
            cboxhi.coarsen(ratio);
            cfablo.resize(cboxlo,1);
            cfabhi.resize(cboxhi,1);
            cfablo.setVal(0.0);
            cfabhi.setVal(0.0);
            //
            // Compute divgp on fine grid edges in regions defined
            // by reglo and reghi.  Fabs are set to zero outside region.
            //
            FORT_SRDGPHI(pfab.dataPtr(),ARLIM(pfab_lo),ARLIM(pfab_hi),
                         sfab.dataPtr(),ARLIM(sfab_lo),ARLIM(sfab_hi),
                         ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                         reglo.loVect(),reglo.hiVect(),
                         dx_fine,&mult,&is_rz);
            FORT_SRDGPHI(pfab.dataPtr(),ARLIM(pfab_lo),ARLIM(pfab_hi),
                         sfab.dataPtr(),ARLIM(sfab_lo),ARLIM(sfab_hi),
                         ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                         reghi.loVect(),reghi.hiVect(),
                         dx_fine,&mult,&is_rz);

            int n_comp   = 1;
            int set_comp = 0;

            for (int i = 0; i < Pgrids.length(); i++)
            {
                Box overlap_lo = reglo & Pgrids[i];
                if (overlap_lo.ok())
                    ffablo.setVal(0.,overlap_lo,set_comp,n_comp);

                fine_geom.periodicShift(reglo, Pgrids[i], pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    Box overlap_lo_per(Pgrids[i]);
                    overlap_lo_per.shift(pshifts[iiv]);
                    overlap_lo_per &= reglo;
                    ffablo.setVal(0.,overlap_lo_per,set_comp,n_comp);
                }

                Box overlap_hi = reghi & Pgrids[i];
                if (overlap_hi.ok())
                    ffabhi.setVal(0.,overlap_hi,set_comp,n_comp);

                fine_geom.periodicShift(reghi, Pgrids[i], pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    Box overlap_hi_per(Pgrids[i]);
                    overlap_hi_per.shift(pshifts[iiv]);
                    overlap_hi_per &= reghi;
                    ffabhi.setVal(0.,overlap_hi_per,set_comp,n_comp);
                }
            }
            //
            // Coarsen edge value.
            //
            const int* clo = cboxlo.loVect();
            const int* chi = cboxlo.hiVect();
            FORT_SRCRSEREG(ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                           cfablo.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            clo = cboxhi.loVect();
            chi = cboxhi.hiVect();
            FORT_SRCRSEREG(ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                           cfabhi.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            //
            // Intersect and add to registers.
            //
            increment(cfablo);
            increment(cfabhi);
            //
            // Intersect and add to registers.
            //
            if (crse_geom.isAnyPeriodic())
            {
                crse_geom.periodicShift(crse_node_domain, cboxlo, pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    cfablo.shift(pshifts[iiv]);
                    increment(cfablo);
                    cfablo.shift(-pshifts[iiv]);
                }
                crse_geom.periodicShift(crse_node_domain, cboxhi, pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    cfabhi.shift(pshifts[iiv]);
                    increment(cfabhi);
                    cfabhi.shift(-pshifts[iiv]);
                }
            }
        }
    }
}
