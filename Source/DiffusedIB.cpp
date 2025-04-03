
#include <AMReX_Math.H>
#include <AMReX_Print.H>
#include <DiffusedIB.H>

#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#include <AMReX_Utility.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_FillPatchUtil.H>
#include <iamr_constants.H>

#include <filesystem>
#include <sstream>
namespace fs = std::filesystem;

#define GHOST_CELLS 2

using namespace amrex;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                     global variable                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define LOCAL_LEVEL 0

const Vector<std::string> direction_str{"X","Y","Z"};

namespace ParticleProperties{
    Vector<Real> _x{}, _y{}, _z{}, _rho{};
    Vector<Real> Vx{}, Vy{}, Vz{};
    Vector<Real> Ox{}, Oy{}, Oz{};
    Vector<Real> _radius;
    Real rd{0.0};
    Vector<int> TLX{}, TLY{},TLZ{},RLX{},RLY{},RLZ{};
    int euler_finest_level{0};
    int euler_velocity_index{0};
    int euler_force_index{0};
    Real euler_fluid_rho{0.0};
    int verbose{0};
    int loop_ns{2};
    int loop_solid{1};
    int Uhlmann{0};

    Vector<Real> GLO, GHI;
    int start_step{-1};
    int collision_model{0};

    int write_freq{1};
    bool init_particle_from_file{false};

    GpuArray<Real, 3> plo{0.0,0.0,0.0}, phi{0.0,0.0,0.0}, dx{0.0, 0.0, 0.0};
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                     other function                            */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void nodal_phi_to_pvf(MultiFab& pvf, const MultiFab& phi_nodal)
{

    // amrex::Print() << "In the nodal_phi_to_pvf\n";

    pvf.setVal(0.0);

    // Only set the valid cells of pvf
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(pvf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& pvffab   = pvf.array(mfi);
        auto const& pnfab = phi_nodal.array(mfi);
        amrex::ParallelFor(bx, [pvffab, pnfab]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real num = 0.0;
            for(int kk=k; kk<=k+1; kk++) {
                for(int jj=j; jj<=j+1; jj++) {
                    for(int ii=i; ii<=i+1; ii++) {
                        num += (-pnfab(ii,jj,kk)) * nodal_phi_to_heavi(-pnfab(ii,jj,kk));
                    }
                }
            }
            Real deo = 0.0;
            for(int kk=k; kk<=k+1; kk++) {
                for(int jj=j; jj<=j+1; jj++) {
                    for(int ii=i; ii<=i+1; ii++) {
                        deo += std::abs(pnfab(ii,jj,kk));
                    }
                }
            }
            pvffab(i,j,k) = num / (deo + 1.e-12);
        });
    }

}

void calculate_phi_nodal(MultiFab& phi_nodal, kernel& current_kernel)
{
    phi_nodal.setVal(0.0);

    amrex::Real Xp = current_kernel.location[0];
    amrex::Real Yp = current_kernel.location[1];
    amrex::Real Zp = current_kernel.location[2];
    amrex::Real Rp = current_kernel.radius;

    // Only set the valid cells of phi_nodal
    for (MFIter mfi(phi_nodal,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& pnfab = phi_nodal.array(mfi);
        auto dx = ParticleProperties::dx;
        auto plo = ParticleProperties::plo;
        amrex::ParallelFor(bx, [=]
            AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real Xn = i * dx[0] + plo[0];
                Real Yn = j * dx[1] + plo[1];
                Real Zn = k * dx[2] + plo[2];

                pnfab(i,j,k) = std::sqrt( (Xn - Xp)*(Xn - Xp)
                        + (Yn - Yp)*(Yn - Yp)  + (Zn - Zp)*(Zn - Zp)) - Rp;
                pnfab(i,j,k) = pnfab(i,j,k) / Rp;

            }
        );
    }
}

// May use ParReduce later, https://amrex-codes.github.io/amrex/docs_html/GPU.html#multifab-reductions
void CalculateSumU_cir (RealVect& sum,
                        const MultiFab& E,
                        const MultiFab& pvf,
                        int EulerVelIndex)
{
    auto const& E_data = E.const_arrays();
    auto const& pvf_data = pvf.const_arrays();
    const Real d = Math::powi<3>(ParticleProperties::dx[0]);
    amrex::GpuTuple<Real, Real, Real> tmpSum = ParReduce(TypeList<ReduceOpSum,ReduceOpSum,ReduceOpSum>{}, TypeList<Real, Real, Real>{},E, IntVect{0},
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> amrex::GpuTuple<Real, Real, Real>{
        auto E_ = E_data[box_no];
        auto pvf_ = pvf_data[box_no];
        return {
            E_(i, j, k, EulerVelIndex    ) * d * pvf_(i,j,k),
            E_(i, j, k, EulerVelIndex + 1) * d * pvf_(i,j,k),
            E_(i, j, k, EulerVelIndex + 2) * d * pvf_(i,j,k)
        };
    });
    sum[0] = amrex::get<0>(tmpSum);
    sum[1] = amrex::get<1>(tmpSum);
    sum[2] = amrex::get<2>(tmpSum);
}

void CalculateSumT_cir (RealVect& sum,
                        const MultiFab& E,
                        const MultiFab& pvf,
                        const RealVect pLoc,
                        int EulerVelIndex)
{
    auto plo = ParticleProperties::plo;
    auto dx = ParticleProperties::dx;

    auto const& E_data = E.const_arrays();
    auto const& pvf_data = pvf.const_arrays();
    const Real d = Math::powi<3>(ParticleProperties::dx[0]);
    amrex::GpuTuple<Real, Real, Real> tmpSum = ParReduce(TypeList<ReduceOpSum,ReduceOpSum,ReduceOpSum>{}, TypeList<Real, Real, Real>{},E, IntVect{0},
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> amrex::GpuTuple<Real, Real, Real>{
        auto E_ = E_data[box_no];
        auto pvf_ = pvf_data[box_no];

        Real x = plo[0] + i*dx[0] + 0.5*dx[0];
        Real y = plo[1] + j*dx[1] + 0.5*dx[1];
        Real z = plo[2] + k*dx[2] + 0.5*dx[2];

        Real vx = E_(i, j, k, EulerVelIndex    );
        Real vy = E_(i, j, k, EulerVelIndex + 1);
        Real vz = E_(i, j, k, EulerVelIndex + 2);

        RealVect tmp = RealVect(x - pLoc[0], y - pLoc[1], z - pLoc[2]).crossProduct(RealVect(vx, vy, vz));

        return {
            tmp[0] * d * pvf_(i, j, k),
            tmp[1] * d * pvf_(i, j, k),
            tmp[2] * d * pvf_(i, j, k)
        };
    });
    sum[0] = amrex::get<0>(tmpSum);
    sum[1] = amrex::get<1>(tmpSum);
    sum[2] = amrex::get<2>(tmpSum);
}

[[nodiscard]] AMREX_FORCE_INLINE
Real cal_momentum(Real rho, Real radius)
{
    return 8.0 * Math::pi<Real>() * rho * Math::powi<5>(radius) / 15.0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void deltaFunction(Real xf, Real xp, Real h, Real& value, DELTA_FUNCTION_TYPE type)
{
    Real rr = amrex::Math::abs(( xf - xp ) / h);

    switch (type) {
    case DELTA_FUNCTION_TYPE::FOUR_POINT_IB:
        if(rr >= 0 && rr < 1 ){
            value = 1.0 / 8.0 * ( 3.0 - 2.0 * rr + std::sqrt( 1.0 + 4.0 * rr - 4 * Math::powi<2>(rr))) / h;
        }else if (rr >= 1 && rr < 2) {
            value = 1.0 / 8.0 * ( 5.0 - 2.0 * rr - std::sqrt( -7.0 + 12.0 * rr - 4 * Math::powi<2>(rr))) / h;
        }else {
            value = 0;
        }
        break;
    case DELTA_FUNCTION_TYPE::THREE_POINT_IB:
        if(rr >= 0.5 && rr < 1.5){
            value = 1.0 / 6.0 * ( 5.0 - 3.0 * rr - std::sqrt( - 3.0 * Math::powi<2>( 1 - rr) + 1.0 )) / h;
        }else if (rr >= 0 && rr < 0.5) {
            value = 1.0 / 3.0 * ( 1.0 + std::sqrt( 1.0 - 3 * Math::powi<2>(rr))) / h;
        }else {
            value = 0;
        }
        break;
    default:
        break;
    }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                    mParticle member function                  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//loop all particels
void mParticle::InteractWithEuler(MultiFab &EulerVel, 
                                  MultiFab &EulerForce, 
                                  Real dt,
                                  DELTA_FUNCTION_TYPE type)
{
    if (verbose) amrex::Print() << "[Particle] mParticle::InteractWithEuler\n";
    // clear time , start record
    spend_time = 0;
    auto InteractWithEulerStart = ParallelDescriptor::second();

    MultiFab EulerForceTmp(EulerForce.boxArray(), EulerForce.DistributionMap(), 3, EulerForce.nGrow());
    //clean preStep's IB_porperties 
    for(auto& kernel : particle_kernels) {
        kernel.ib_force.scale(0.0);
        kernel.ib_moment.scale(0.0);
    }

    //for 1 -> Ns
    int loop = ParticleProperties::loop_ns;
    BL_ASSERT(loop > 0);
    while(loop > 0){
        if(verbose) amrex::Print() << "[Particle] Ns loop index : " << loop << "\n";
        
        EulerForce.setVal(0.0);

        for(kernel& kernel : particle_kernels){
            InitialWithLargrangianPoints(kernel); // Initialize markers for a specific particle
            ResetLargrangianPoints();
            EulerForceTmp.setVal(0.0);
            auto ib_force = kernel.ib_force;
            auto ib_moment = kernel.ib_moment;
            kernel.ib_force.scale(0.0); // clear kernel ib_force
            kernel.ib_moment.scale(0.0); // clear kernel ib_moment

            VelocityInterpolation(EulerVel, type);
            ComputeLagrangianForce(dt, kernel);
            ForceSpreading(EulerForceTmp, kernel, type);
            MultiFab::Add(EulerForce, EulerForceTmp, 0, 0, 3, EulerForce.nGrow());

            kernel.ib_force += ib_force;
            kernel.ib_moment += ib_moment;
        }
        VelocityCorrection(EulerVel, EulerForce, dt);
        loop--;
    }
    spend_time += ParallelDescriptor::second() - InteractWithEulerStart;
}

void mParticle::InitParticles(const Vector<Real>& x,
                              const Vector<Real>& y,
                              const Vector<Real>& z,
                              const Vector<Real>& rho_s,
                              const Vector<Real>& Vx,
                              const Vector<Real>& Vy,
                              const Vector<Real>& Vz,
                              const Vector<Real>& Ox,
                              const Vector<Real>& Oy,
                              const Vector<Real>& Oz,
                              const Vector<int>& TLXt,
                              const Vector<int>& TLYt,
                              const Vector<int>& TLZt,
                              const Vector<int>& RLXt,
                              const Vector<int>& RLYt,
                              const Vector<int>& RLZt,
                              const Vector<Real>& radius,
                              Real h,
                              Real gravity,
                              int _verbose)
{
    verbose = _verbose;
    if (verbose) amrex::Print() << "[Particle] mParticle::InitParticles\n";

    m_gravity[2] = gravity;

    //pre judge
    if(!((x.size() == y.size()) && (x.size() == z.size()))){
        Print() << "particle's position container are all different size";
        return;
    }

    //all the particles have different radius
    for(int index = 0; index < x.size(); index++){
        int real_index;
        // if initial with input file, initialized by [0] data
        if(ParticleProperties::init_particle_from_file){
            real_index = 0;
        }else{
            real_index = index;
        }

        kernel mKernel;
        mKernel.id = index + 1;
        mKernel.location[0] = x[index];
        mKernel.location[1] = y[index];
        mKernel.location[2] = z[index];
        mKernel.velocity[0] = Vx[real_index];
        mKernel.velocity[1] = Vy[real_index];
        mKernel.velocity[2] = Vz[real_index];
        mKernel.omega[0] = Ox[real_index];
        mKernel.omega[1] = Oy[real_index];
        mKernel.omega[2] = Oz[real_index];

        // use current state to initialize old state
        mKernel.location_old = mKernel.location;
        mKernel.velocity_old = mKernel.velocity;
        mKernel.omega_old = mKernel.omega;

        mKernel.TL[0] = TLXt[real_index];
        mKernel.TL[1] = TLYt[real_index];
        mKernel.TL[2] = TLZt[real_index];
        mKernel.RL[0] = RLXt[real_index];
        mKernel.RL[1] = RLYt[real_index];
        mKernel.RL[2] = RLZt[real_index];
        mKernel.rho = rho_s[real_index];
        mKernel.radius = radius[real_index];
        mKernel.Vp = Math::pi<Real>() * 4 / 3 * Math::powi<3>(radius[real_index]);

        //int Ml = static_cast<int>( Math::pi<Real>() / 3 * (12 * Math::powi<2>(mKernel.radius / h)));
        //Real dv = Math::pi<Real>() * h / 3 / Ml * (12 * mKernel.radius * mKernel.radius + h * h);
        int Ml = static_cast<int>((amrex::Math::powi<3>(mKernel.radius - (ParticleProperties::rd - 0.5) * h)
               - amrex::Math::powi<3>(mKernel.radius - (ParticleProperties::rd + 0.5) * h))/(3.*h*h*h/4./Math::pi<Real>()));
        Real dv = (amrex::Math::powi<3>(mKernel.radius - (ParticleProperties::rd - 0.5) * h)
               - amrex::Math::powi<3>(mKernel.radius - (ParticleProperties::rd + 0.5) * h))/(3.*Ml/4./Math::pi<Real>());
        mKernel.ml = Ml;
        mKernel.dv = dv;
        if( Ml > max_largrangian_num ) max_largrangian_num = Ml;

        Real phiK = 0;
        for(int marker_index = 0; marker_index < Ml; marker_index++){
            Real Hk = -1.0 + 2.0 * (marker_index) / ( Ml - 1.0);
            Real thetaK = std::acos(Hk);    
            if(marker_index == 0 || marker_index == (Ml - 1)){
                phiK = 0;
            }else {
                phiK = std::fmod( phiK + 3.809 / std::sqrt(Ml) / std::sqrt( 1 - Math::powi<2>(Hk)) , 2 * Math::pi<Real>());
            }
            mKernel.phiK.push_back(phiK);
            mKernel.thetaK.push_back(thetaK);
        }

        particle_kernels.emplace_back(mKernel);

        if (verbose) amrex::Print() << "h: " << h << ", Ml: " << Ml << ", D: " << Math::powi<3>(h) << " gravity : " << gravity << "\n"
                                    << "Kernel : " << index << ": Location (" << x[index] << ", " << y[index] << ", " << z[index] 
                                    << "), Velocity : (" << mKernel.velocity[0] << ", " << mKernel.velocity[1] << ", "<< mKernel.velocity[2] 
                                    << "), Radius: " << mKernel.radius << ", Ml: " << Ml << ", dv: " << dv << ", Rho: " << mKernel.rho << "\n";
    }
    //collision box generate
    m_Collision.SetGeometry(RealVect(ParticleProperties::GLO), RealVect(ParticleProperties::GHI),particle_kernels[0].radius, h);
}

void mParticle::InitialWithLargrangianPoints(const kernel& current_kernel){

    if (verbose) amrex::Print() << "mParticle::InitialWithLargrangianPoints\n";
    for(mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        if(np == 0) continue;
        auto *particles = pti.GetArrayOfStructs().data();

        const auto location = current_kernel.location;
        const auto radius = current_kernel.radius;
        const auto* phiK = current_kernel.phiK.dataPtr();
        const auto* thetaK = current_kernel.thetaK.dataPtr();

        amrex::ParallelFor( np, [=]
            AMREX_GPU_DEVICE (int i) noexcept {
                auto id = particles[i].id();
                particles[i].pos(0) = location[0] + radius * std::sin(thetaK[id - 1]) * std::cos(phiK[id - 1]);
                particles[i].pos(1) = location[1] + radius * std::sin(thetaK[id - 1]) * std::sin(phiK[id - 1]);
                particles[i].pos(2) = location[2] + radius * std::cos(thetaK[id - 1]);
            }
        );
    }
    // Redistribute the markers after updating their locations
    mContainer->Redistribute();
    if (verbose) {
        amrex::Print() << "[particle] : particle num :" << mContainer->TotalNumberOfParticles() << "\n";
        mContainer->WriteAsciiFile(amrex::Concatenate("particle", 1));
    }
}

template <typename P = Particle<numAttri>>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void VelocityInterpolation_cir(int p_iter, P const& p, Real& Up, Real& Vp, Real& Wp,
                               Array4<Real const> const& E, int EulerVIndex,
                               const int *lo, const int *hi, 
                               GpuArray<Real, AMREX_SPACEDIM> const& plo,
                               GpuArray<Real, AMREX_SPACEDIM> const& dx,
                               DELTA_FUNCTION_TYPE type)
{

    //std::cout << "lo " << lo[0] << " " << lo[1] << " "<< lo[2] << "\n";
    //std::cout << "hi " << hi[0] << " " << hi[1] << " "<< hi[2] << "\n";

    const Real d = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

    const Real lx = (p.pos(0) - plo[0]) / dx[0]; // x
    const Real ly = (p.pos(1) - plo[1]) / dx[1]; // y
    const Real lz = (p.pos(2) - plo[2]) / dx[2]; // z

    int i = static_cast<int>(Math::floor(lx)); // i
    int j = static_cast<int>(Math::floor(ly)); // j
    int k = static_cast<int>(Math::floor(lz)); // k

    //std::cout << "p_iter " << p_iter << " p.pos(0): " << p.pos(0) << " p.pos(1): " << p.pos(1) << " p.pos(2): " << p.pos(2) << "\n";

    // std::cout << "d: " << d << "\n"
    //         << "lx: " << lx << ", ly: " << ly << ", lz: " << lz << "\n"
    //         << "i: " << i << ", j: " << j << ", k: " << k << std::endl;

    Up = 0;
    Vp = 0;
    Wp = 0;
    //Euler to largrangian
    for(int ii = -2; ii < 3; ii++){
        for(int jj = -2; jj < 3; jj++){
            for(int kk = -2; kk < 3; kk ++){
                Real tU, tV, tW;
                const Real xi = plo[0] + (i + ii) * dx[0] + dx[0]/2;
                const Real yj = plo[1] + (j + jj) * dx[1] + dx[1]/2;
                const Real kz = plo[2] + (k + kk) * dx[2] + dx[2]/2;
                deltaFunction( p.pos(0), xi, dx[0], tU, type);
                deltaFunction( p.pos(1), yj, dx[1], tV, type);
                deltaFunction( p.pos(2), kz, dx[2], tW, type);
                const Real delta_value = tU * tV * tW;
                Up += delta_value * E(i + ii, j + jj, k + kk, EulerVIndex    ) * d;
                Vp += delta_value * E(i + ii, j + jj, k + kk, EulerVIndex + 1) * d;
                Wp += delta_value * E(i + ii, j + jj, k + kk, EulerVIndex + 2) * d;
            }
        }
    }
}

void mParticle::VelocityInterpolation(MultiFab &EulerVel,
                                      DELTA_FUNCTION_TYPE type)//
{
    if (verbose) amrex::Print() << "\tmParticle::VelocityInterpolation\n";

    //amrex::Print() << "euler_finest_level " << euler_finest_level << std::endl;
    const auto& gm = mContainer->GetParGDB()->Geom(LOCAL_LEVEL);
    auto plo = gm.ProbLoArray();
    auto dx = gm.CellSizeArray();
    // attention
    // velocity ghost cells will be up-to-date
    EulerVel.FillBoundary(ParticleProperties::euler_velocity_index, 3, gm.periodicity());

    const int EulerVelocityIndex = ParticleProperties::euler_velocity_index;

    for(mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        
        const Box& box = pti.validbox();
        
        auto& particles = pti.GetArrayOfStructs();
        auto *p_ptr = particles.data();
        const Long np = pti.numParticles();

        auto& attri = pti.GetAttribs();
        auto* Up = attri[P_ATTR::U_Marker].data();
        auto* Vp = attri[P_ATTR::V_Marker].data();
        auto* Wp = attri[P_ATTR::W_Marker].data();
        const auto& E = EulerVel.array(pti);

        amrex::ParallelFor(np, [=] 
        AMREX_GPU_DEVICE (int i) noexcept{
            VelocityInterpolation_cir(i, p_ptr[i], Up[i], Vp[i], Wp[i], E, EulerVelocityIndex, box.loVect(), box.hiVect(), plo, dx, type);
        });
    }
    if (verbose) mContainer->WriteAsciiFile(amrex::Concatenate("particle", 2));
    //amrex::Abort("stop here!");
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ForceSpreading_cic (P const& p,
                         Real Px,
                         Real Py,
                         Real Pz,
                         ParticleReal& fxP,
                         ParticleReal& fyP,
                         ParticleReal& fzP,
                         ParticleReal& mxP,
                         ParticleReal& myP,
                         ParticleReal& mzP,
                         Array4<Real> const& E,
                         int EulerForceIndex,
                         Real dv,
                         GpuArray<Real,AMREX_SPACEDIM> const& plo,
                         GpuArray<Real,AMREX_SPACEDIM> const& dx,
                         DELTA_FUNCTION_TYPE type)
{
    //const Real d = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    //plo to ii jj kk
    Real lx = (p.pos(0) - plo[0]) / dx[0];
    Real ly = (p.pos(1) - plo[1]) / dx[1];
    Real lz = (p.pos(2) - plo[2]) / dx[2];

    int i = static_cast<int>(Math::floor(lx));
    int j = static_cast<int>(Math::floor(ly));
    int k = static_cast<int>(Math::floor(lz));
    fxP *= dv;
    fyP *= dv;
    fzP *= dv;
    RealVect moment = RealVect((p.pos(0) - Px), (p.pos(1) - Py), (p.pos(2) - Pz)).crossProduct(RealVect(fxP, fyP, fzP));
    mxP = moment[0];
    myP = moment[1];
    mzP = moment[2];
    //lagrangian to Euler
    for(int ii = -2; ii < +3; ii++){
        for(int jj = -2; jj < +3; jj++){
            for(int kk = -2; kk < +3; kk ++){
                Real tU, tV, tW;
                const Real xi =plo[0] + (i + ii) * dx[0] + dx[0]/2;
                const Real yj =plo[1] + (j + jj) * dx[1] + dx[1]/2;
                const Real kz =plo[2] + (k + kk) * dx[2] + dx[2]/2;
                deltaFunction( p.pos(0), xi, dx[0], tU, type);
                deltaFunction( p.pos(1), yj, dx[1], tV, type);
                deltaFunction( p.pos(2), kz, dx[2], tW, type);
                Real delta_value = tU * tV * tW;
                Gpu::Atomic::AddNoRet(&E(i + ii, j + jj, k + kk, EulerForceIndex  ), delta_value * fxP);
                Gpu::Atomic::AddNoRet(&E(i + ii, j + jj, k + kk, EulerForceIndex+1), delta_value * fyP);
                Gpu::Atomic::AddNoRet(&E(i + ii, j + jj, k + kk, EulerForceIndex+2), delta_value * fzP);
            }
        }
    }
}

void mParticle::ForceSpreading(MultiFab & EulerForce,
                               kernel& kernel,
                               DELTA_FUNCTION_TYPE type)
{
    if (verbose) amrex::Print() << "\tmParticle::ForceSpreading\n";
    const auto& gm = mContainer->GetParGDB()->Geom(LOCAL_LEVEL);
    auto plo = gm.ProbLoArray();
    auto dxi = gm.CellSizeArray();
    for(mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        const auto& particles = pti.GetArrayOfStructs();
        auto Uarray = EulerForce[pti].array();
        auto& attri = pti.GetAttribs();

        auto *const fxP_ptr = attri[P_ATTR::Fx_Marker].data();
        auto *const fyP_ptr = attri[P_ATTR::Fy_Marker].data();
        auto *const fzP_ptr = attri[P_ATTR::Fz_Marker].data();
        auto *const mxP_ptr = attri[P_ATTR::Mx_Marker].data();
        auto *const myP_ptr = attri[P_ATTR::My_Marker].data();
        auto *const mzP_ptr = attri[P_ATTR::Mz_Marker].data();
        const auto *const p_ptr = particles().data();

        auto loc_ptr = kernel.location;
        auto dv = kernel.dv;
        auto force_index = ParticleProperties::euler_force_index;
        amrex::ParallelFor(np, [=]
        AMREX_GPU_DEVICE (int i) noexcept{
            ForceSpreading_cic(p_ptr[i], loc_ptr[0], loc_ptr[1], loc_ptr[2],
                               fxP_ptr[i], fyP_ptr[i], fzP_ptr[i], 
                               mxP_ptr[i], myP_ptr[i], mzP_ptr[i], 
                               Uarray, force_index, dv, plo, dxi, type);
        });
    }
    //barrier for sync;
    amrex::ParallelDescriptor::Barrier();

    using pc = mParticleContainer::SuperParticleType;
    // Each Processor
    auto fx = amrex::ReduceSum( *mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{return p.rdata(P_ATTR::Fx_Marker);});
    auto fy = amrex::ReduceSum( *mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{return p.rdata(P_ATTR::Fy_Marker);});
    auto fz = amrex::ReduceSum( *mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{return p.rdata(P_ATTR::Fz_Marker);});
    auto mx = amrex::ReduceSum( *mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{return p.rdata(P_ATTR::Mx_Marker);});
    auto my = amrex::ReduceSum( *mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{return p.rdata(P_ATTR::My_Marker);});
    auto mz = amrex::ReduceSum( *mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{return p.rdata(P_ATTR::Mz_Marker);});
    // MPI sum reduce -> current particle all IB force and moment
    amrex::ParallelAllReduce::Sum(fx, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(fy, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(fz, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(mx, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(my, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(mz, ParallelDescriptor::Communicator());

    kernel.ib_force = {fx, fy, fz};
    kernel.ib_moment = {mx, my, mz};

    EulerForce.SumBoundary(ParticleProperties::euler_force_index, 3, gm.periodicity());

    // if (false) {
    //     // Check the Multifab
    //     // Open a file stream for output
    //     std::ofstream outFile("EulerForce.txt");

    //     // Check the Multifab
    //     // for (MFIter mfi(EulerForce, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    //     for (MFIter mfi(EulerForce, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    //     {
    //         const Box& bx = mfi.validbox();
    //         outFile << "Box: " << bx << "\n"
    //                 << "From: (" << bx.smallEnd(0) << ", " << bx.smallEnd(1) << ", " << bx.smallEnd(2) << ") "
    //                 << "To: (" << bx.bigEnd(0) << ", " << bx.bigEnd(1) << ", " << bx.bigEnd(2) << ")\n";

    //         Array4<Real> const& a = EulerForce[mfi].array();

    //         // CPU context or illustrative purposes only
    //         for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
    //             for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
    //                 for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
    //                     // This print statement is for demonstration and should not be used in actual GPU code.
    //                     outFile << "Processing i: " << i << ", j: " << j << ", k: " << k << " " << a(i,j,k,0) << " " << a(i,j,k,1) << " " << a(i,j,k,2) << "\n";
    //                 }
    //             }
    //         }
    //     }

    //     // Close the file when done
    //     outFile.close();
    // }

}

void mParticle::ResetLargrangianPoints()
{
    if (verbose) amrex::Print() << "\tmParticle::ResetLargrangianPoints\n";

    for(mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        auto& attri = pti.GetAttribs();

        auto *const vUP_ptr = attri[P_ATTR::U_Marker].data();
        auto *const vVP_ptr = attri[P_ATTR::V_Marker].data();
        auto *const vWP_ptr = attri[P_ATTR::W_Marker].data();
        auto *const fxP_ptr = attri[P_ATTR::Fx_Marker].data();
        auto *const fyP_ptr = attri[P_ATTR::Fy_Marker].data();
        auto *const fzP_ptr = attri[P_ATTR::Fz_Marker].data();
        auto *const mxP_ptr = attri[P_ATTR::Mx_Marker].data();
        auto *const myP_ptr = attri[P_ATTR::My_Marker].data();
        auto *const mzP_ptr = attri[P_ATTR::Mz_Marker].data();
        amrex::ParallelFor(np, [=]
        AMREX_GPU_DEVICE (int i) noexcept{
            vUP_ptr[i] = 0.0;
            vVP_ptr[i] = 0.0;
            vWP_ptr[i] = 0.0;
            fxP_ptr[i] = 0.0;
            fyP_ptr[i] = 0.0;
            fzP_ptr[i] = 0.0;
            mxP_ptr[i] = 0.0;
            myP_ptr[i] = 0.0;
            mzP_ptr[i] = 0.0;
        });
    }
}

void mParticle::UpdateParticles(int iStep,
                                Real time,
                                const MultiFab& Euler_old, 
                                const MultiFab& Euler,
                                MultiFab& phi_nodal, 
                                MultiFab& pvf, 
                                Real dt)
{
    if (verbose) amrex::Print() << "mParticle::UpdateParticles\n";
    // start record
    auto UpdateParticlesStart = ParallelDescriptor::second();

    //Particle Collision calculation
    DoParticleCollision(ParticleProperties::collision_model);
    
    MultiFab AllParticlePVF(pvf.boxArray(), pvf.DistributionMap(), pvf.nComp(), pvf.nGrow());
    AllParticlePVF.setVal(0.0);
    
    //continue condition 6DOF
    for(auto& kernel : particle_kernels){

        calculate_phi_nodal(phi_nodal, kernel);
        nodal_phi_to_pvf(pvf, phi_nodal);

        // // fixed particle
        // if( ( kernel.TL.sum() == 0 ) &&
        //     ( kernel.RL.sum() == 0 ) ) {
        //     amrex::Print() << "Particle (" << kernel.id << ") is fixed\n";
        //     MultiFab::Add(AllParticlePVF, pvf, 0, 0, 1, 0); // do not copy ghost cell values
        //     continue;
        // }

        int ncomp = pvf.nComp();
        int ngrow = pvf.nGrow();
        MultiFab pvf_old(pvf.boxArray(), pvf.DistributionMap(), ncomp, ngrow);
        MultiFab::Copy(pvf_old, pvf, 0, 0, ncomp, ngrow);

        // bool at_least_one_free_trans_motion = ( kernel.TL[0] == 2 ) || 
        //                                       ( kernel.TL[1] == 2 ) ||
        //                                       ( kernel.TL[2] == 2 );
        // bool at_least_one_free_rot_motion   = ( kernel.RL[0] == 2 ) || 
        //                                       ( kernel.RL[1] == 2 ) ||
        //                                       ( kernel.RL[2] == 2 );

        int loop = ParticleProperties::loop_solid;

        while (loop > 0 && iStep > ParticleProperties::start_step) {

            // if(at_least_one_free_trans_motion) {
                kernel.sum_u_new.scale(0.0);
                kernel.sum_u_old.scale(0.0);
                // sum U
                CalculateSumU_cir(kernel.sum_u_new, Euler, pvf, ParticleProperties::euler_velocity_index);
                CalculateSumU_cir(kernel.sum_u_old, Euler_old, pvf_old, ParticleProperties::euler_velocity_index);
                amrex::ParallelAllReduce::Sum(kernel.sum_u_new.dataPtr(), 3, amrex::ParallelDescriptor::Communicator());
                amrex::ParallelAllReduce::Sum(kernel.sum_u_old.dataPtr(), 3, amrex::ParallelDescriptor::Communicator());
            // }

            // if(at_least_one_free_rot_motion) {
                kernel.sum_t_new.scale(0.0);
                kernel.sum_t_old.scale(0.0);
                // sum T
                CalculateSumT_cir(kernel.sum_t_new, Euler, pvf, kernel.location, ParticleProperties::euler_velocity_index);
                CalculateSumT_cir(kernel.sum_t_old, Euler_old, pvf_old, kernel.location, ParticleProperties::euler_velocity_index);
                amrex::ParallelAllReduce::Sum(kernel.sum_t_new.dataPtr(), 3, amrex::ParallelDescriptor::Communicator());
                amrex::ParallelAllReduce::Sum(kernel.sum_t_old.dataPtr(), 3, amrex::ParallelDescriptor::Communicator());
            // }

            // 6DOF
            if(ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()){

                for(auto idir : {0,1,2})
                {
                    //TL
                    if (kernel.TL[idir] == 0) {
                        kernel.velocity[idir] = 0.0;
                    }
                    else if (kernel.TL[idir] == 1) {
                        kernel.location[idir] = kernel.location_old[idir] + (kernel.velocity[idir] + kernel.velocity_old[idir]) * dt * 0.5;
                    }
                    else if (kernel.TL[idir] == 2) {
                        if(!ParticleProperties::Uhlmann){
                            kernel.velocity[idir] = kernel.velocity_old[idir]
                                                + ((kernel.sum_u_new[idir] - kernel.sum_u_old[idir]) * ParticleProperties::euler_fluid_rho / dt 
                                                - kernel.ib_force[idir] * ParticleProperties::euler_fluid_rho
                                                + m_gravity[idir] * (kernel.rho - ParticleProperties::euler_fluid_rho) * kernel.Vp 
                                                + kernel.Fcp[idir]) * dt / kernel.rho / kernel.Vp ;
                        }else{
                            //Uhlmann
                            kernel.velocity[idir] = kernel.velocity_old[idir]
                                                + (ParticleProperties::euler_fluid_rho / kernel.Vp /(ParticleProperties::euler_fluid_rho - kernel.rho)*kernel.ib_force[idir]
                                                + m_gravity[idir]) * dt;
                        }
                        kernel.location[idir] = kernel.location_old[idir] + (kernel.velocity[idir] + kernel.velocity_old[idir]) * dt * 0.5;
                    }
                    else {
                        amrex::Print() << "Particle (" << kernel.id << ") has wrong TL"<< direction_str[idir] <<" value\n";
                        amrex::Abort("Stop here!");
                    }
                    //RL
                    if (kernel.RL[idir] == 0) {
                        kernel.omega[idir] = 0.0;
                    }
                    else if (kernel.RL[idir] == 1) {
                    }
                    else if (kernel.RL[idir] == 2) {
                        if(!ParticleProperties::Uhlmann){
                            kernel.omega[idir] = kernel.omega_old[idir]
                                            + ((kernel.sum_t_new[idir] - kernel.sum_t_old[idir]) * ParticleProperties::euler_fluid_rho / dt
                                            - kernel.ib_moment[idir] * ParticleProperties::euler_fluid_rho
                                            + kernel.Tcp[idir]) * dt / cal_momentum(kernel.rho, kernel.radius);
                        }else{
                            //Uhlmann
                            kernel.omega[idir] = kernel.omega_old[idir]
                                            + ParticleProperties::euler_fluid_rho /(ParticleProperties::euler_fluid_rho - kernel.rho) * kernel.ib_moment[idir] * kernel.dv
                                            / cal_momentum(kernel.rho, kernel.radius) * kernel.rho * dt;
                        }
                    }
                    else {
                        amrex::Print() << "Particle (" << kernel.id << ") has wrong RL"<< direction_str[idir] <<" value\n";
                        amrex::Abort("Stop here!");
                    }

                }
            }
            ParallelDescriptor::Bcast(&kernel.location[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.location_old[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.velocity[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.velocity_old[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.omega[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.omega_old[0],3,ParallelDescriptor::IOProcessorNumber());
        
            loop--;

            if (loop > 0) {
                calculate_phi_nodal(phi_nodal, kernel);
                nodal_phi_to_pvf(pvf, phi_nodal);
            }

        }
        
        RecordOldValue(kernel);
        MultiFab::Add(AllParticlePVF, pvf, 0, 0, 1, 0); // do not copy ghost cell values
    }
    // calculate the pvf based on the information of all particles
    MultiFab::Copy(pvf, AllParticlePVF, 0, 0, 1, pvf.nGrow());

    spend_time += ParallelDescriptor::second() - UpdateParticlesStart;
    ParallelDescriptor::ReduceRealMax(spend_time);
    amrex::Print() << "[DIBM] IB and update particle, step : "<< iStep <<", time : " << spend_time << "\n";

    int particle_write_freq = ParticleProperties::write_freq;
    if (iStep % particle_write_freq == 0) {
        for(auto kernel: particle_kernels) 
            WriteIBForceAndMoment(iStep, time, dt, kernel);
    }

    if (verbose) mContainer->WriteAsciiFile(amrex::Concatenate("particle", 4));
}

void mParticle::DoParticleCollision(int model)
{
    if(particle_kernels.size() < 2 ) return ;

    if (verbose) amrex::Print() << "\tmParticle::DoParticleCollision\n";
    
    if(ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()){
        for(const auto& kernel : particle_kernels){
            m_Collision.InsertParticle(kernel.location, kernel.velocity, kernel.radius, kernel.rho);
        }
        
        m_Collision.takeModel(model);

        for(auto & particle_kernel : particle_kernels){
            particle_kernel.Fcp = m_Collision.Particles.front().preForece 
                                * particle_kernel.Vp * particle_kernel.rho * m_gravity.vectorLength();
            m_Collision.Particles.pop_front();
        }
    }
    for(auto& kernel : particle_kernels){
        ParallelDescriptor::Bcast(kernel.Fcp.dataPtr(), 3, ParallelDescriptor::IOProcessorNumber());
    }
}

void mParticle::ComputeLagrangianForce(Real dt, 
                                       const kernel& kernel)
{
    
    if (verbose) amrex::Print() << "\tmParticle::ComputeLagrangianForce\n";

    Real Ub = kernel.velocity[0];
    Real Vb = kernel.velocity[1];
    Real Wb = kernel.velocity[2];
    Real Px = kernel.location[0];
    Real Py = kernel.location[1];
    Real Pz = kernel.location[2];

    for(mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        auto& attri = pti.GetAttribs();
        auto const* p_ptr = pti.GetArrayOfStructs().data();

        auto* Up = attri[P_ATTR::U_Marker].data();
        auto* Vp = attri[P_ATTR::V_Marker].data();
        auto* Wp = attri[P_ATTR::W_Marker].data();
        auto *FxP = attri[P_ATTR::Fx_Marker].data();
        auto *FyP = attri[P_ATTR::Fy_Marker].data();
        auto *FzP = attri[P_ATTR::Fz_Marker].data();

        amrex::ParallelFor(np,
        [=] AMREX_GPU_DEVICE (int i) noexcept{
            auto Ur = (kernel.omega).crossProduct(RealVect(p_ptr[i].pos(0) - Px, p_ptr[i].pos(1) - Py, p_ptr[i].pos(2) - Pz));
            FxP[i] = (Ub + Ur[0] - Up[i])/dt; //
            FyP[i] = (Vb + Ur[1] - Vp[i])/dt; //
            FzP[i] = (Wb + Ur[2] - Wp[i])/dt; //
        });
    }
    if (verbose) mContainer->WriteAsciiFile(amrex::Concatenate("particle", 3));
}

void mParticle::VelocityCorrection(amrex::MultiFab &Euler, amrex::MultiFab &EulerForce, Real dt) const
{
    if(verbose) amrex::Print() << "\tmParticle::VelocityCorrection\n";
    MultiFab::Saxpy(Euler, dt, EulerForce, ParticleProperties::euler_force_index, ParticleProperties::euler_velocity_index, 3, 0); //VelocityCorrection
}

void mParticle::RecordOldValue(kernel& kernel)
{
    kernel.location_old = kernel.location;
    kernel.velocity_old = kernel.velocity;
    kernel.omega_old = kernel.omega;
}

void mParticle::WriteParticleFile(int index)
{
    mContainer->WriteAsciiFile(amrex::Concatenate("particle", index));
}

void mParticle::WriteIBForceAndMoment(int step, amrex::Real time, amrex::Real dt, kernel& current_kernel)
{
    
    if(amrex::ParallelDescriptor::MyProc() != ParallelDescriptor::IOProcessorNumber()) return; 

    std::string file("IB_Particle_" + std::to_string(current_kernel.id) + ".csv");
    std::ofstream out_ib_force;

    std::string head;
    if(!fs::exists(file)){
        head = "iStep,time,X,Y,Z,Vx,Vy,Vz,Rx,Ry,Rz,Fx,Fy,Fz,Mx,My,Mz,Fcpx,Fcpy,Fcpz,Tcpx,Tcpy,Tcpz,SumUx,SumUy,SumUz,SumTx,SumTy,SumTz\n";
    }else{
        head = "";
    }

    out_ib_force.open(file, std::ios::app);
    if(!out_ib_force.is_open()){
        amrex::Print() << "[Particle] write particle file error , step: " << step;
    }else{
        out_ib_force << head << step << "," << time << "," 
                     << current_kernel.location[0] << "," << current_kernel.location[1] << "," << current_kernel.location[2] << ","
                     << current_kernel.velocity[0] << "," << current_kernel.velocity[1] << "," << current_kernel.velocity[2] << ","
                     << current_kernel.omega[0] << "," << current_kernel.omega[1] << "," << current_kernel.omega[2] << ","
                     << current_kernel.ib_force[0] << "," << current_kernel.ib_force[1] << "," << current_kernel.ib_force[2] << "," 
                     << current_kernel.ib_moment[0] << "," << current_kernel.ib_moment[1] << "," << current_kernel.ib_moment[2] << ","
                     << current_kernel.Fcp[0] << "," << current_kernel.Fcp[1] << "," << current_kernel.Fcp[2] << ","
                     << current_kernel.Tcp[0] << "," << current_kernel.Tcp[1] << "," << current_kernel.Tcp[2] << ","
                     << (current_kernel.sum_u_new[0] - current_kernel.sum_u_old[0])/dt << ","
                     << (current_kernel.sum_u_new[1] - current_kernel.sum_u_old[1])/dt << ","
                     << (current_kernel.sum_u_new[2] - current_kernel.sum_u_old[2])/dt << ","
                     << (current_kernel.sum_t_new[0] - current_kernel.sum_t_old[0])/dt << ","
                     << (current_kernel.sum_t_new[1] - current_kernel.sum_t_old[1])/dt << ","
                     << (current_kernel.sum_t_new[2] - current_kernel.sum_t_old[2])/dt << "\n";
    }
    out_ib_force.close();
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                    Particles member function                  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void Particles::create_particles(const Geometry &gm,
                                 const DistributionMapping & dm,
                                 const BoxArray & ba)
{
    amrex::Print() << "[Particle] : create Particle Container\n";
    if(particle->mContainer != nullptr){
        delete particle->mContainer;
        particle->mContainer = nullptr;
    }
    particle->mContainer = new mParticleContainer(gm, dm, ba);

    //get particle tile
    std::pair<int, int> key{0,0};
    auto& particleTileTmp = particle->mContainer->GetParticles(0)[key];
    //insert markers
    if ( ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber() ) {
        //insert particle's markers
        //Real phiK = 0;
        for(int marker_index = 0; marker_index < particle->particle_kernels[0].ml; marker_index++){
            //insert code
            mParticleContainer::ParticleType markerP;
            markerP.id() = marker_index + 1;
            markerP.cpu() = ParallelDescriptor::MyProc();
            markerP.pos(0) = particle->particle_kernels[0].location[0];
            markerP.pos(1) = particle->particle_kernels[0].location[1];
            markerP.pos(2) = particle->particle_kernels[0].location[2];

            std::array<ParticleReal, numAttri> Marker_attr;
            Marker_attr[U_Marker] = 0.0;
            Marker_attr[V_Marker] = 0.0;
            Marker_attr[W_Marker] = 0.0;
            Marker_attr[Fx_Marker] = 0.0;
            Marker_attr[Fy_Marker] = 0.0;
            Marker_attr[Fz_Marker] = 0.0;

            particleTileTmp.push_back(markerP);
            particleTileTmp.push_back_real(Marker_attr);
        }
    }
    particle->mContainer->Redistribute(); // Still needs to redistribute here! 

    ParticleProperties::plo = gm.ProbLoArray();
    ParticleProperties::phi = gm.ProbHiArray();
    ParticleProperties::dx = gm.CellSizeArray();
}

mParticle* Particles::get_particles()
{
    return particle;
}


void Particles::init_particle(Real gravity, Real h)
{
    amrex::Print() << "[Particle] : create Particle's kernel\n";
    particle = new mParticle;
    if(particle != nullptr){
        isInitial = true;
        particle->InitParticles(
            ParticleProperties::_x, 
            ParticleProperties::_y, 
            ParticleProperties::_z,
            ParticleProperties::_rho,
            ParticleProperties::Vx,
            ParticleProperties::Vy,
            ParticleProperties::Vz,
            ParticleProperties::Ox,
            ParticleProperties::Oy,
            ParticleProperties::Oz,
            ParticleProperties::TLX,
            ParticleProperties::TLY,
            ParticleProperties::TLZ,
            ParticleProperties::RLX,
            ParticleProperties::RLY,
            ParticleProperties::RLZ,
            ParticleProperties::_radius,
            h,
            gravity,
            ParticleProperties::verbose);
    }

}

void Particles::Restart(Real gravity, Real h, int iStep)
{
    amrex::Print() << "[Particle] : restart Particle's kernel, step :" << iStep << "\n"
                   << "\tstart read particle csv file , default name is IB_Particle_x.csv\n" 
                   << "\tdo not delete those file before \"restart\"\n\n";
    delete particle;
    particle = new mParticle;
            particle->InitParticles(
            ParticleProperties::_x, 
            ParticleProperties::_y, 
            ParticleProperties::_z,
            ParticleProperties::_rho,
            ParticleProperties::Vx,
            ParticleProperties::Vy,
            ParticleProperties::Vz,
            ParticleProperties::Ox,
            ParticleProperties::Oy,
            ParticleProperties::Oz,
            ParticleProperties::TLX,
            ParticleProperties::TLY,
            ParticleProperties::TLZ,
            ParticleProperties::RLX,
            ParticleProperties::RLY,
            ParticleProperties::RLZ,
            ParticleProperties::_radius,
            h,
            gravity,
            ParticleProperties::verbose);
    //deal in IO processor
    //start read csv file
    for(auto& kernel : particle->particle_kernels){
        //filename
        if(amrex::ParallelDescriptor::MyProc() == amrex::ParallelDescriptor::IOProcessorNumber()){
            std::string fileName = "IB_Particle_" + std::to_string(kernel.id) + ".csv";
            std::string tmpfile = "tmp" + fileName;
            //file stream
            std::ifstream particle_data(fileName);
            std::ofstream particle_file(tmpfile);
            // open state
            if(!particle_data.is_open() || !particle_file.is_open()){
                amrex::Abort("\tCan not open particle file : " + fileName);
            }
            std::string lineData;
            int line{0};
            while(std::getline(particle_data, lineData)){
                line++;
                particle_file << lineData << "\n";
                if(line <= iStep) {
                    continue;
                }
                //old location
                //iStep,time,X,Y,Z,Vx,Vy,Vz,Rx,Ry,Rz,Fx,Fy,Fz,Mx,My,Mz,Fcpx,Fcpy,Fcpz,Tcpx,Tcpy,Tcpz
                if(line == iStep + 1){
                    std::stringstream ss(lineData);
                    std::string data;
                    std::vector<amrex::Real> dataStruct;
                    while(std::getline(ss, data, ',')){
                        dataStruct.emplace_back(std::stod(data));
                    }
                    kernel.location[0] = dataStruct[2];
                    kernel.location[1] = dataStruct[3];
                    kernel.location[2] = dataStruct[4];
                    kernel.velocity[0] = dataStruct[5];
                    kernel.velocity[1] = dataStruct[6];
                    kernel.velocity[2] = dataStruct[7];
                    kernel.omega[0] = dataStruct[8];
                    kernel.omega[1] = dataStruct[9];
                    kernel.omega[2] = dataStruct[10];

                    kernel.location_old = kernel.location;
                    kernel.velocity_old = kernel.velocity;
                    kernel.omega_old    = kernel.omega;
                    break;
                }
                else
                    break;
            }
            particle_data.close();
            particle_file.close();
            std::remove(fileName.c_str());
            std::rename(tmpfile.c_str(), fileName.c_str());
        }
        ParallelDescriptor::Bcast(&kernel.location[0], 3, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&kernel.velocity[0], 3,ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&kernel.omega[0], 3,ParallelDescriptor::IOProcessorNumber());

        ParallelDescriptor::Bcast(&kernel.location_old[0], 3, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&kernel.velocity_old[0], 3,ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&kernel.omega_old[0], 3,ParallelDescriptor::IOProcessorNumber());
    }

    isInitial = true;
}

void Particles::Initialize()
{
    ParmParse pp("particle");

    std::string particle_inputfile;
    std::string particle_init_file;
    pp.get("input",particle_inputfile);
    
    if(!particle_inputfile.empty()){
        ParmParse p_file(particle_inputfile);
        p_file.query("init", particle_init_file);
        p_file.getarr("x",          ParticleProperties::_x);
        p_file.getarr("y",          ParticleProperties::_y);
        p_file.getarr("z",          ParticleProperties::_z);
        p_file.getarr("rho",        ParticleProperties::_rho);
        p_file.getarr("velocity_x", ParticleProperties::Vx);
        p_file.getarr("velocity_y", ParticleProperties::Vy);
        p_file.getarr("velocity_z", ParticleProperties::Vz);
        p_file.getarr("omega_x",    ParticleProperties::Ox);
        p_file.getarr("omega_y",    ParticleProperties::Oy);
        p_file.getarr("omega_z",    ParticleProperties::Oz);
        p_file.getarr("TLX",        ParticleProperties::TLX);
        p_file.getarr("TLY",        ParticleProperties::TLY);
        p_file.getarr("TLZ",        ParticleProperties::TLZ);
        p_file.getarr("RLX",        ParticleProperties::RLX);
        p_file.getarr("RLY",        ParticleProperties::RLY);
        p_file.getarr("RLZ",        ParticleProperties::RLZ);
        p_file.getarr("radius",     ParticleProperties::_radius);
        p_file.query("RD",          ParticleProperties::rd);
        p_file.query("LOOP_NS",     ParticleProperties::loop_ns);
        p_file.query("LOOP_SOLID",  ParticleProperties::loop_solid);
        p_file.query("verbose",     ParticleProperties::verbose);
        p_file.query("start_step",  ParticleProperties::start_step);
        p_file.query("Uhlmann",     ParticleProperties::Uhlmann);
        p_file.query("collision_model", ParticleProperties::collision_model);
        p_file.query("write_freq",  ParticleProperties::write_freq);
        
        ParmParse ns("ns");
        ns.get("fluid_rho",      ParticleProperties::euler_fluid_rho);
        
        ParmParse level_parse("amr");
        level_parse.get("max_level", ParticleProperties::euler_finest_level);

        ParmParse geometry_parse("geometry");
        geometry_parse.getarr("prob_lo", ParticleProperties::GLO);
        geometry_parse.getarr("prob_hi", ParticleProperties::GHI);
        amrex::Print() << "[Particle] : Reading partilces cfg file : " << particle_inputfile << "\n"
                       << "             Particle's level : " << ParticleProperties::euler_finest_level << "\n";

        if(!particle_init_file.empty()){
            ParticleProperties::init_particle_from_file = true;
            //clear particle position container
            ParticleProperties::_x.clear();
            ParticleProperties::_y.clear();
            ParticleProperties::_z.clear();
            // parse particle's location data
            std::ifstream init_particle(particle_init_file);
            std::string line_data;
            while(std::getline(init_particle, line_data)){
                // id x_location y_location z_location
                std::istringstream line(line_data);
                std::vector<std::string> str_tokne;
                std::string token;
                while(line >> token){
                    str_tokne.push_back(token);
                }

                ParticleProperties::_x.push_back(std::stod(str_tokne[0]));
                ParticleProperties::_y.push_back(std::stod(str_tokne[1]));
                ParticleProperties::_z.push_back(std::stod(str_tokne[2]));
            }
            ParticleProperties::_x.shrink_to_fit();
            ParticleProperties::_y.shrink_to_fit();
            ParticleProperties::_z.shrink_to_fit();
            amrex::Print() << "             intial Particle by file : " << particle_init_file
                           << "             particle's size : " << ParticleProperties::_x.size() << "\n";
        }

    }else {
        amrex::Abort("[Particle] : can't read particles settings, pls check your config file \"particle.input\"");
    }
}

int Particles::ParticleFinestLevel()
{
    return ParticleProperties::euler_finest_level;
}
