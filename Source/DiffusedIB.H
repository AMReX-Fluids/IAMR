#ifndef DIFFUSEDIB_H
#define DIFFUSEDIB_H

#include <AMReX_Particles.H>
#include <AMReX_MultiFabUtil.H>

#include <AMReX_RealVect.H>
#include <fstream>

#include "Collision.H"
using namespace amrex;
// using deltaFuncType = std::function<AMREX_GPU_HOST_DEVICE void(Real, Real, Real, Real&)>;

AMREX_INLINE AMREX_GPU_DEVICE
Real nodal_phi_to_heavi(Real phi)
{
    if (phi <= 0.0) {
        return 0.0;
    }
    else {
        return 1.0;
    }
}

void nodal_phi_to_pvf(MultiFab& pvf, const MultiFab& phi_nodal);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                     particle and markers                      */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

enum P_ATTR{
    U_Marker = 0,
    V_Marker,
    W_Marker,
    Fx_Marker,
    Fy_Marker,
    Fz_Marker,
    Mx_Marker,
    My_Marker,
    Mz_Marker,
    numAttri
};

enum DELTA_FUNCTION_TYPE{
    FOUR_POINT_IB = 0,
    THREE_POINT_IB
};

struct kernel{
    int id;
    RealVect velocity{0.0,0.0,0.0};
    RealVect location{0.0,0.0,0.0};
    RealVect omega{0.0,0.0,0.0};

    RealVect velocity_old{0.0,0.0,0.0};
    RealVect location_old{0.0,0.0,0.0};
    RealVect omega_old{0.0,0.0,0.0};

    RealVect varphi{0.0,0.0,0.0};
    Real radius{0.0};
    Real rho{0.0};
    int ml{0};
    Real dv{0.0};
    Real Vp{0.0};
    RealVect ib_force{0.0,0.0,0.0};
    RealVect ib_moment{0.0, 0.0, 0.0};

    RealVect sum_u_new{0.0,0.0,0.0};
    RealVect sum_u_old{0.0,0.0,0.0};
    RealVect sum_t_new{0.0,0.0,0.0};
    RealVect sum_t_old{0.0,0.0,0.0};

    RealVect Fcp{0.0,0.0,0.0};
    RealVect Tcp{0.0,0.0,0.0};

    amrex::Gpu::DeviceVector<Real> phiK;
    amrex::Gpu::DeviceVector<Real> thetaK;
    
    IntVect TL{0, 0, 0}, RL{0, 0, 0};
};

class mParIter : public ParIter<0, 0, numAttri>{
public:
    using amrex::ParIter<0, 0, numAttri>::ParIter;
    using RealVector = amrex::ParIter<0, 0, numAttri>::ContainerType::RealVector;
 
    [[nodiscard]] const std::array<RealVector, numAttri>& GetAttribs () const {
        return GetStructOfArrays().GetRealData();
    }

    [[nodiscard]] const RealVector& GetAttribs (int comp) const {
        return GetStructOfArrays().GetRealData(comp);
    }

    std::array<RealVector, numAttri>& GetAttribs () {
        return GetStructOfArrays().GetRealData();
    }

    RealVector& GetAttribs (int comp) {
        return GetStructOfArrays().GetRealData(comp);
    }
};

void deltaFunction(Real xf, Real xp, Real h, Real& value);

using mParticleContainer = ParticleContainer<0, 0, numAttri>;

class mParticle
{
public:
    explicit mParticle() = default;

    void InitParticles(const Vector<Real>& x,
                       const Vector<Real>& y,
                       const Vector<Real>& z,
                       const Vector<Real>& rho_s,
                       const Vector<Real>& Vx,
                       const Vector<Real>& Vy,
                       const Vector<Real>& Vz,
                       const Vector<Real>& Ox,
                       const Vector<Real>& Oy,
                       const Vector<Real>& Oz,
                       const Vector<int>& TLX,
                       const Vector<int>& TLY,
                       const Vector<int>& TLZ,
                       const Vector<int>& RLX,
                       const Vector<int>& RLY,
                       const Vector<int>& RLZ,
                       const Vector<Real>& radius,
                       Real h,
                       Real gravity,
                       int verbose = 0);

    void InteractWithEuler(MultiFab &EulerVel, MultiFab &EulerForce, Real dt = 0.1, DELTA_FUNCTION_TYPE type = FOUR_POINT_IB);

    void WriteParticleFile(int index);

    void InitialWithLargrangianPoints(const kernel& current_kernel);

    void ResetLargrangianPoints();

    void VelocityInterpolation(amrex::MultiFab &Euler, DELTA_FUNCTION_TYPE type);

    void ComputeLagrangianForce(Real dt, const kernel& kernel);

    void ForceSpreading(amrex::MultiFab &Euler, kernel& kernel, DELTA_FUNCTION_TYPE type);

    void VelocityCorrection(amrex::MultiFab &Euler, amrex::MultiFab &EulerForce, Real dt) const;

    void UpdateParticles(int iStep, Real time, const MultiFab& Euler_old, const MultiFab& Euler, MultiFab& phi_nodal, MultiFab& pvf, Real dt);

    void DoParticleCollision(int model);

    static void WriteIBForceAndMoment(int step, amrex::Real time, amrex::Real dt, kernel& current_kernel);

    void RecordOldValue(kernel& kernel);

    Vector<kernel> particle_kernels;

    mParticleContainer *mContainer{nullptr};

    ParticleCollision m_Collision;

    int max_largrangian_num = 0;

    uint32_t ib_force_file_index = 0;

    RealVect m_gravity{0.0,0.0,0.0};

    int verbose = 0;

    Real spend_time;
};

class Particles{
public:
    static void create_particles(const Geometry &gm,
                                 const DistributionMapping & dm,
                                 const BoxArray & ba);

    static mParticle* get_particles();
    static void init_particle(Real gravity, Real h);
    static void Restart(Real gravity, Real h, int iStep);
    static void Initialize();
    static int ParticleFinestLevel();
    static void define_para(const Vector<Real>& x,
                            const Vector<Real>& y,
                            const Vector<Real>& z,
                            Real rho_s,
                            Real radius,
                            Real rho_f, 
                            int force_index, 
                            int velocity_index,
                            int finest_level);

    inline static bool isInitial{false};
private:
    inline static mParticle* particle = nullptr;
};

#endif
