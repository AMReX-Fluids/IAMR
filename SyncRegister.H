#ifndef _SYNCREGISTER_H_
#define _SYNCREGISTER_H_

//
// $Id: SyncRegister.H,v 1.5 1997-10-01 01:03:16 car Exp $
//

#include <PArray.H>
#include <BoxArray.H>
#include <FArrayBox.H>
#include <MultiFab.H>
#include <Geometry.H>
#include <BndryRegister.H>
#include <BCRec.H>

class SyncRegister : public BndryRegister
{
public:
    SyncRegister();
    SyncRegister(const BoxArray& fine_boxes,
		 IntVect ref_ratio, int fine_lev);
    ~SyncRegister();

    void define(const BoxArray& fine_boxes,
		IntVect ref_ratio, int fine_lev);

      // access functions
    IntVect refRatio() const {
	return ratio;
    }
    int fineLevel() const {
	return fine_level;
    }
    int crseLevel() const {
	return fine_level-1;
    }
    const BoxArray& coarsenedBoxes() const {
	return grids;
    }

    void CrseDVInit(const MultiFab& V, const Geometry& geom, 
		    int is_rz, int ** crse_bc, REAL mult = 1.0);

    void FineDVAdd(const MultiFab& V, 
                   const REAL* dx_fine, const Geometry& crse_geom,
		   int is_rz, int ** fine_bc, REAL mult);

    void CrseDsdtAdd(const MultiFab& dsdt, const Geometry& geom,
			 int is_rz, int ** crse_bc, 
                         int lowfix, int hifix, REAL mult = 1.0);

    void FineDsdtAdd(const MultiFab& dsdt, const Geometry& geom,
		   int is_rz, int ** fine_bc, 
                   int lowfix, int hifix,REAL mult);

    void CompDVAdd(const MultiFab& V, 
                   const BoxArray& fine_boxes,
                   const REAL* dx_fine, 
                   const Geometry& fine_geom, 
                   const Geometry& crse_geom, 
                   int is_rz, int ** fine_bc, REAL mult);

    void CrseLPhiAdd(const MultiFab& Phi, const MultiFab& sigma,
		     const Geometry& geom, int is_rz, REAL mult = 1.0);

    void FineLPhiAdd(const MultiFab& Phi, const MultiFab& sigma,
		     const REAL* dx_fine, const Geometry& geom, 
                     int is_rz, REAL mult);

    void CompLPhiAdd(const MultiFab& Phi, const MultiFab& sigma,
		     const BoxArray& fine_boxes, 
                     const REAL* dx_fine, 
                     const Geometry& fine_geom, 
                     const Geometry& crse_geom, 
                     int is_rz, REAL mult);

    void InitRHS(MultiFab& rhs, const Geometry& geom, const BCRec* phys_bc);

      // compute sum of values (overlap values only counted once)
      // WARNING: this is a diagnostic and is not efficient
    REAL sum();

    void print(ostream& os);

protected:

    IntVect ratio;
    int fine_level;
    FabSet    bndry_mask[2*BL_SPACEDIM];

      // helper member functions
    void increment(const FArrayBox& fab);
};

#endif
