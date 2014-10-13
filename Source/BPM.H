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

//
// Particle mesh for BPM method
//
class BPMMesh
{
  std::deque<Spring> springs;
  std::multimap<int,int> pmap;

public:
  void InitFromAsciiFile(const std::string& file);
};

#endif
