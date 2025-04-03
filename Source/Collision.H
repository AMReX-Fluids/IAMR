#ifndef COLLISION_H
#define COLLISION_H

#include <AMReX_RealVect.H>
#include <array>
#include <list>

//collision object type
enum CollisionType{
    COLLISION_PARTICLE = 0,
    COLLISION_PLANE,
};

//particle
typedef struct{
    CollisionType type{COLLISION_PARTICLE};
    amrex::RealVect location{0., 0., 0.};
    amrex::RealVect velocity{0., 0., 0.};
    amrex::RealVect preForece{0., 0., 0.};
    amrex::RealVect fc{0., 0., 0.};
    amrex::Real radius{0.};
    amrex::Real rho{0.};
}CollisionParticle;

// mesh cell to store particle
typedef struct{
    bool hasParticle{false};
    std::list<CollisionParticle*> collectParticle;
}CollisionCell;

// collision paire to calculate
typedef struct{
    CollisionParticle* solid1;
    CollisionParticle* solid2;
}CollisionPair;

class ParticleCollision{
public:
    void takeModel(int model);

    void DKTModel();

    explicit ParticleCollision() = default;
    //generate background meshs
    void SetGeometry(amrex::RealVect gm_lo, amrex::RealVect gm_hi, amrex::Real size, amrex::Real l);
    //initial Collision particle;
    void InsertParticle(amrex::RealVect location, amrex::RealVect velocity, amrex::Real radius, amrex::Real rho);
    //insert collision pair to CollisionCollector
    void InsertCollision(CollisionParticle* p1, CollisionParticle* p2);
    //generator collision pair
    void GenerateCollisionPairs();
    //Resolve collision
    void ResolveCollisionPairs();
    //all the Geometry split mesh
    std::vector<CollisionCell> *Cells;
    std::list<CollisionParticle> Particles;
    //mesh x,y,z cell size, loop for i,j,k
    int Nx, Ny, Nz;
    amrex::Real cell_size, mesh_size;
    //lo,hi
    amrex::RealVect lo, hi;
    //store collision type
    std::list<CollisionPair> CollisionCollector;

};

#endif