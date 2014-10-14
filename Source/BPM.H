#ifndef _BPM_H_
#define _BPM_H_

#include <map>
#include <deque>
#include <fstream>
#include <iostream>
#include <string>

#include <REAL.H>
#include <Particles.H>

//
// Spring
//
struct Spring
{
  int p1, p2;                   // particle indexes
  Real l;                       // equilibrium/rest length
  Real k;                       // spring constant
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
  ParticleSpringMap pmap;

  int verbose;

public:
  BPM(int verbose=0) : verbose(verbose) { }

  void InitFromAsciiFile(const std::string& file);
  void ComputeForces(ParticleContainer<2*BL_SPACEDIM>& particles, MultiFab& force, int lev, const int fcomp);

};


#endif
