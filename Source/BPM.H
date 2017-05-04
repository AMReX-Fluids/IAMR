#ifndef _BPM_H_
#define _BPM_H_

#include <map>
#include <deque>
#include <fstream>
#include <iostream>
#include <string>

#include <AMReX_REAL.H>
#include <AMReX_Particles.H>

//
// Spring
//
struct Spring
{
  int p1, p2;                   // particle indexes
  amrex::Real l;                       // equilibrium/rest length
  amrex::Real k;                       // spring constant
};

typedef std::deque<Spring> Springs;   // list of springs
typedef std::multimap<int,int> ParticleSpringMap;
typedef std::pair<ParticleSpringMap::iterator,ParticleSpringMap::iterator> ParticleSpringRange;

//
// Particle mesh for BPM method
//
class BPM
{
  Springs springs;
  ParticleSpringMap particle_to_spring_map;

  int verbose;

public:

  amrex::Real diameter;

  BPM(int verbose=0) : verbose(verbose) { }

  void InitFromAsciiFile(const std::string& file);
  void ComputeForces(amrex::ParticleContainer<2*BL_SPACEDIM>& particles,
                     amrex::MultiFab& force, amrex::MultiFab& divforce, amrex::Geometry& geom,
                     int lev);

};


#endif
