#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

   void compute_slopes (
      const int* lo, const int* hi, const int* nc,
      const amrex::Real* u, const int* ulo, const int* uhi, 
      D_DECL(amrex::Real* xslopes, amrex::Real* yslopes,
	     amrex::Real* zslopes),
      const int* slo, const int* shi, 
      const int* domlo, const int* domhi,
      D_DECL(const int* ubc, const int* vbc, const int* wbc));
      // const int* bc_ilo_type, const int* bc_ihi_type, 
      // const int* bc_jlo_type, const int* bc_jhi_type,
      // const int* bc_klo_type, const int* bc_khi_type,
      // const int* nghost );

   void compute_slopes_eb (
      const int* lo, const int* hi, const int* nc,
      const amrex::Real* u, const int* ulo, const int* uhi, 
      D_DECL(amrex::Real* xslopes, amrex::Real* yslopes,
	     amrex::Real* zslopes),
      const int* slo, const int* shi,
      const void* flags, const int* flo, const int* fhi,
      const int* domlo, const int* domhi,
      D_DECL(const int* ubc, const int* vbc, const int* wbc));
      // const int* bc_ilo_type, const int* bc_ihi_type, 
      // const int* bc_jlo_type, const int* bc_jhi_type,
      // const int* bc_klo_type, const int* bc_khi_type,
      // const int* nghost );

   void compute_velocity_at_faces (
      const int* lo, const int* hi,
      const amrex::Real* u_mac, const int* ulo, const int* uhi, 
      const amrex::Real* v_mac, const int* vlo, const int* vhi,
#if (BL_SPACEDIM == 3)
      const amrex::Real* w_mac, const int* wlo, const int* whi,
#endif
      const amrex::Real* vel   ,  const int* vello, const int* velhi, 
      D_DECL(const amrex::Real* xslopes, const amrex::Real* yslopes, 
	     const amrex::Real* zslopes),
      const int* slo, const int* shi, 
      D_DECL(const int* ubc, const int* vbc, const int* wbc),
      // const int* bc_ilo_type, const int* bc_ihi_type, 
      // const int* bc_jlo_type, const int* bc_jhi_type,
      // const int* bc_klo_type, const int* bc_khi_type,
      // const int* ng, 
      const int* domlo, const int* domhi );

   void compute_velocity_at_x_faces_eb (
      const int* lo, const int* hi,
      amrex::Real* u_mac,     const int* ulo,   const int* uhi, 
      const amrex::Real* vel   ,    const int* vello, const int* velhi, 
      const amrex::Real* slopes,    const int* slo,   const int* shi,
      const amrex::Real* areafrac,  const int* alo,   const int* ahi,
      const amrex::Real* cent,      const int* clo,   const int* chi,
      const void*        flags,     const int* flo,   const int* fhi,
      const int* ubc,
      //      const int* bc_ilo_type,       const int* bc_ihi_type, 
      //      const int* ng, 
      const int* domlo, const int* domhi );

   void compute_velocity_at_y_faces_eb (
      const int* lo, const int* hi,
            amrex::Real* v_mac,     const int* vlo,   const int* vhi, 
      const amrex::Real* vel   ,    const int* vello, const int* velhi, 
      const amrex::Real* slopes,    const int* slo,   const int* shi,
      const amrex::Real* areafrac,  const int* alo,   const int* ahi,
      const amrex::Real* cent,      const int* clo,   const int* chi,
      const void*        flags,     const int* flo,   const int* fhi,
      const int* vbc,
      //const int* bc_jlo_type,       const int* bc_jhi_type, 
      //      const int* ng, 
      const int* domlo, const int* domhi );

   void compute_velocity_at_z_faces_eb (
      const int* lo, const int* hi,
            amrex::Real* w_mac,     const int* wlo,   const int* whi, 
      const amrex::Real* vel   ,    const int* vello, const int* velhi, 
      const amrex::Real* slopes,    const int* slo,   const int* shi,
      const amrex::Real* areafrac,  const int* alo,   const int* ahi,
      const amrex::Real* cent,      const int* clo,   const int* chi,
      const void*        flags,     const int* flo,   const int* fhi,
      const int* wbc,
      //      const int* bc_klo_type,       const int* bc_khi_type, 
      //      const int* ng, 
      const int* domlo, const int* domhi );   

  void compute_ugradu (
	const int* lo, const int* hi,
	amrex::Real* rhs, const int* rlo, const int* rhi, 
	const amrex::Real* vel   ,  const int* vello, const int* velhi,
	const amrex::Real* u, const int* ulo, const int* uhi,
	const amrex::Real* v, const int* vlo, const int* vhi,	
	const amrex::Real* w, const int* wlo, const int* whi,
	const amrex::Real* xslopes, 
	const amrex::Real* yslopes, 
	const amrex::Real* zslopes, const int* slo, const int* shi, 
	const int* domlo, const int* domhi,
	// const int* bc_ilo_type, const int* bc_ihi_type, 
	// const int* bc_jlo_type, const int* bc_jhi_type,
	// const int* bc_klo_type, const int* bc_khi_type,
	const amrex::Real* dx, const int* ng
        );


   void compute_ugradu_eb (
	const int* lo, const int* hi,
	amrex::Real* ugradu, const int* glo, const int* ghi, 
	const amrex::Real* vel   ,  const int* vello, const int* velhi,
	const amrex::Real* u, const int* ulo, const int* uhi,
	const amrex::Real* v, const int* vlo, const int* vhi,	
	const amrex::Real* w, const int* wlo, const int* whi,
        const amrex::Real* afrac_x,const int*  axlo,const int*  axhi, 
        const amrex::Real* afrac_y,const int*  aylo,const int*  ayhi, 
        const amrex::Real* afrac_z,const int*  azlo,const int*  azhi, 
        const amrex::Real* cent_x,const int*   cxlo,const int*  cxhi, 
        const amrex::Real* cent_y,const int*   cylo,const int*  cyhi, 
        const amrex::Real* cent_z,const int*   czlo,const int*  czhi, 
        const void* flags, const int*  flo, const int*   fhi, 
        const amrex::Real* vfrac, const int* vflo, const int*  vfhi, 
        const amrex::Real* bcent, const int* blo, const int*  bhi, 
	const amrex::Real* xslopes, 
	const amrex::Real* yslopes, 
	const amrex::Real* zslopes, const int* slo, const int* shi, 
	const int* domlo, const int* domhi,
	// const int* bc_ilo_type, const int* bc_ihi_type, 
	// const int* bc_jlo_type, const int* bc_jhi_type,
	// const int* bc_klo_type, const int* bc_khi_type,
	const amrex::Real* dx, const int* ng
      );

    void compute_divuc (
	const int* lo, const int* hi,
	amrex::Real* rhs, const int* rlo, const int* rhi, 
	const amrex::Real* vel   ,  const int* vello, const int* velhi,
	const amrex::Real* u, const int* ulo, const int* uhi,
	const amrex::Real* v, const int* vlo, const int* vhi,	
	const amrex::Real* w, const int* wlo, const int* whi,
	const amrex::Real* xslopes, 
	const amrex::Real* yslopes, 
	const amrex::Real* zslopes, const int* slo, const int* shi, 
	const int* domlo, const int* domhi,
	// const int* bc_ilo_type, const int* bc_ihi_type, 
	// const int* bc_jlo_type, const int* bc_jhi_type,
	// const int* bc_klo_type, const int* bc_khi_type,
	const amrex::Real* dx, const int* ng
        );


   void compute_divuc_eb (
	const int* lo, const int* hi,
	amrex::Real* divuc, const int* glo, const int* ghi, 
	const amrex::Real* vel   ,  const int* vello, const int* velhi,
	const amrex::Real* u, const int* ulo, const int* uhi,
	const amrex::Real* v, const int* vlo, const int* vhi,	
	const amrex::Real* w, const int* wlo, const int* whi,
        const amrex::Real* afrac_x,const int*  axlo,const int*  axhi, 
        const amrex::Real* afrac_y,const int*  aylo,const int*  ayhi, 
        const amrex::Real* afrac_z,const int*  azlo,const int*  azhi, 
        const amrex::Real* cent_x,const int*   cxlo,const int*  cxhi, 
        const amrex::Real* cent_y,const int*   cylo,const int*  cyhi, 
        const amrex::Real* cent_z,const int*   czlo,const int*  czhi, 
        const void* flags, const int*  flo, const int*   fhi, 
        const amrex::Real* vfrac, const int* vflo, const int*  vfhi, 
        const amrex::Real* bcent, const int* blo, const int*  bhi, 
	const amrex::Real* xslopes, 
	const amrex::Real* yslopes, 
	const amrex::Real* zslopes, const int* slo, const int* shi, 
	const int* domlo, const int* domhi,
	// const int* bc_ilo_type, const int* bc_ihi_type, 
	// const int* bc_jlo_type, const int* bc_jhi_type,
	// const int* bc_klo_type, const int* bc_khi_type,
	const amrex::Real* dx, const int* ng
      );


//   void compute_velocity_at_faces (
//       const int* lo, const int* hi,
//       const amrex::Real* vel   ,  const int* vello, const int* velhi,
//       const amrex::Real* u_mac, const int* ulo, const int* uhi, 
//       const amrex::Real* xslopes, const int* slo, const int* shi,
//       const int* ubc, 
//       const amrex::Real* v_mac, const int* vlo, const int* vhi,
//       const amrex::Real* yslopes, const int* vbc, 
// # if (BL_SPACEDIM==3)
//       const amrex::Real* w_mac, const int* wlo, const int* whi, 
//       const amrex::Real* zslopes, const int* wbc,
// #endif
//       // bc arrays that don't require same BC all along boundary
//       // const int* bc_ilo_type, const int* bc_ihi_type, 
//       // const int* bc_jlo_type, const int* bc_jhi_type,
//       // const int* bc_klo_type, const int* bc_khi_type,
//       const int* ng, 
//       const int* domlo, const int* domhi );

//    void compute_velocity_at_x_faces_eb (
//       const int* lo, const int* hi,
//             amrex::Real* u_mac,     const int* ulo,   const int* uhi, 
//       const amrex::Real* vel   ,    const int* vello, const int* velhi, 
//       const amrex::Real* slopes,    const int* slo,   const int* shi,
//       const amrex::Real* areafrac,  const int* alo,   const int* ahi,
//       const amrex::Real* cent,      const int* clo,   const int* chi,
//       const void*        flags,     const int* flo,   const int* fhi,
//       const int* bc_ilo_type,       const int* bc_ihi_type, 
//       const int* ng, 
//       const int* domlo, const int* domhi );

//    void compute_velocity_at_y_faces_eb (
//       const int* lo, const int* hi,
//             amrex::Real* v_mac,     const int* vlo,   const int* vhi, 
//       const amrex::Real* vel   ,    const int* vello, const int* velhi, 
//       const amrex::Real* slopes,    const int* slo,   const int* shi,
//       const amrex::Real* areafrac,  const int* alo,   const int* ahi,
//       const amrex::Real* cent,      const int* clo,   const int* chi,
//       const void*        flags,     const int* flo,   const int* fhi,
//       const int* bc_jlo_type,       const int* bc_jhi_type, 
//       const int* ng, 
//       const int* domlo, const int* domhi );

// #if (BL_SPACEDIM == 3)
//    void compute_velocity_at_z_faces_eb (
//       const int* lo, const int* hi,
//             amrex::Real* w_mac,     const int* wlo,   const int* whi, 
//       const amrex::Real* vel   ,    const int* vello, const int* velhi, 
//       const amrex::Real* slopes,    const int* slo,   const int* shi,
//       const amrex::Real* areafrac,  const int* alo,   const int* ahi,
//       const amrex::Real* cent,      const int* clo,   const int* chi,
//       const void*        flags,     const int* flo,   const int* fhi,
//       const int* bc_klo_type,       const int* bc_khi_type, 
//       const int* ng, 
//       const int* domlo, const int* domhi );   
// #endif
  
//    void compute_slopes (
//       const int* lo, const int* hi,
//       const amrex::Real* u, const int* ulo, const int* uhi, 
//       amrex::Real* xslopes, amrex::Real* yslopes,
// #if (BL_SPACEDIM == 3)
//       amrex::Real* zslopes,
// #endif
//       const int* slo, const int* shi, 
//       const int* domlo, const int* domhi,
//       const int* bc,
//       // const int* bc_ilo_type, const int* bc_ihi_type, 
//       // const int* bc_jlo_type, const int* bc_jhi_type,
//       // const int* bc_klo_type, const int* bc_khi_type,
//       const int* nghost );

//    void compute_slopes_eb (
//       const int* lo, const int* hi,
//       const amrex::Real* u, const int* ulo, const int* uhi, 
//       amrex::Real* xslopes, amrex::Real* yslopes,
// #if (BL_SPACEDIM ==3) 
//       amrex::Real* zslopes,
// #endif
//       const int* slo, const int* shi,
//       const void* flags, const int* flo, const int* fhi,
//       const int* domlo, const int* domhi,
//       const int* bc,
//       // const int* bc_ilo_type, const int* bc_ihi_type, 
//       // const int* bc_jlo_type, const int* bc_jhi_type,
//       // const int* bc_klo_type, const int* bc_khi_type,
//       const int* nghost );


#ifdef __cplusplus
}
#endif
