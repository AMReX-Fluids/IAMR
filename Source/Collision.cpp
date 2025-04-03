#include "Collision.H"
#include <AMReX_Print.H>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                    useful function                            */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                    Collision member function                  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void ParticleCollision::SetGeometry(amrex::RealVect gm_lo, amrex::RealVect gm_hi, amrex::Real size, amrex::Real l)
{
    lo = gm_lo;
    hi = gm_hi;

    Nx = (int)amrex::Math::floor(gm_hi[0] / size);
    Ny = (int)amrex::Math::floor(gm_hi[1] / size);
    Nz = (int)amrex::Math::floor(gm_hi[2] / size);
    //cell size : particle D
    cell_size = size;
    //mesh size : euler mesh size
    mesh_size = l;
    // Cells = new std::vector<CollisionCell>(Nx * Ny * Nz);

    // amrex::Print() << "[Collision] : size (" << Nx << "," << Ny << "," << Nz << ")\n";   
}

void ParticleCollision::InsertParticle(amrex::RealVect location, amrex::RealVect velocity, amrex::Real radius, amrex::Real rho)
{
    CollisionParticle p;
    p.location = location;
    p.velocity = velocity;
    p.radius = radius;
    p.rho = rho;
    p.type = COLLISION_PARTICLE;

    int i = (int)amrex::Math::floor(location[0] / cell_size) - 1;
    int j = (int)amrex::Math::floor(location[1] / cell_size) - 1;
    int k = (int)amrex::Math::floor(location[2] / cell_size) - 1;

    if ( i > Nx || j > Ny || k > Nz) {
        // amrex::Abort("[Collision Error] : this kernel not inside calculate domain");
    }

    Particles.emplace_back(p);
    // (*Cells)[k * Nx * Ny + j * Nx + i].hasParticle = true;
    // (*Cells)[k * Nx * Ny + j * Nx + i].collectParticle.push_back(&(Particles.back()));
}

void ParticleCollision::InsertCollision(CollisionParticle* p1, CollisionParticle* p2)
{
    CollisionPair pair;
    pair.solid1 = p1;
    pair.solid2 = p2;
    CollisionCollector.emplace_back(pair);
}


void ParticleCollision::GenerateCollisionPairs()
{
    CollisionCollector.clear();

    for(auto p : Particles){
        //cal particles's i,j,k
        int i = (int)amrex::Math::floor(p.location[0] / cell_size) - 1;
        int j = (int)amrex::Math::floor(p.location[1] / cell_size) - 1;
        int k = (int)amrex::Math::floor(p.location[2] / cell_size) - 1;

        for(auto ii : {i - 1, i , i + 1}){
            for(auto jj : { j - 1, j, j + 1}){
                for(auto kk : {k - 1, k, k + 1}){
                    if(ii > Nx || jj > Ny || kk > Nz) continue;
                    if(ii < 0 || jj < 0 || kk < 0) continue;
                    for(auto *another : (*Cells)[kk * Nx * Ny + jj * Nx + ii].collectParticle){
                        if(p.location == another->location) continue;
                        InsertCollision(&p, another);
                    }
                }
            }
        }
    }
}

void ParticleCollision::ResolveCollisionPairs()
{
    for(auto pair : CollisionCollector){
        auto* p1 = pair.solid1;
        const auto* p2 = pair.solid2;
        //Dij
        auto dij = p1->location - p2->location;
        //D + dc
        auto d2s = p1->radius + p2->radius + mesh_size;
        //||Dij||
        auto dis = dij.vectorLength();
        if(dis < d2s){
            //without rho_p V_p ||g||
            p1->preForece = dij.scale( amrex::Math::powi<2>((dis - d2s) / mesh_size) * 1E4);
        }
    }
}

void ParticleCollision::DKTModel(){
    amrex::Print() << "\t ParticleCollision : do DKT calculation\n";
    auto& s1 = Particles.front();
    auto& s2 = Particles.back();
    amrex::RealVect dis{s1.location[0] - s2.location[0], s1.location[1] - s2.location[1], s1.location[2] - s2.location[2]};
    double d2s = dis.vectorLength();
    double judge = s1.radius + s2.radius + mesh_size;
    amrex::Print() <<"judge : " << judge <<", distance of s1, s2 : " << d2s << ", dc : " << mesh_size << "\n";
    if(d2s < judge){
        auto mid = std::pow((d2s - judge) / mesh_size, 2);
        s1.preForece = 1e4 * mid / d2s * dis;
        s2.preForece = - s1.preForece;
        amrex::Print() <<"do Collision : s1 :" << s1.preForece << ", s2 :" << s2.preForece << "\n";
    }else{
        s1.preForece.scale(0.);
        s2.preForece.scale(0.);
    }
}

void ParticleCollision::takeModel(int model){
    switch (model) {
    case 1:
        //two sphere DKT Collision model
        DKTModel();
        break;
    default:
        break;
    }
}