
#include <algorithm>
#include <cmath>
#include <iostream>

#include <BPM.H>
#include <BPM_F.H>

using namespace std;

typedef ParticleContainer<2*BL_SPACEDIM>::PBox PBox;
typedef ParticleContainer<2*BL_SPACEDIM>::PMap PMap;
typedef ParticleContainer<2*BL_SPACEDIM>::ParticleType ParticleType;


//
// Compute spring forces given particle locations.
//
void BPM::ComputeForces(ParticleContainer<2*BL_SPACEDIM>& particles, MultiFab& f, MultiFab& divf, Geometry& geom, int lev)
{
  BL_PROFILE("BPM::ComputeForces()");

  const PMap& pmap = particles.GetParticles(lev);

  //
  // Get all particle IDs, locations, and velocities.  These are
  // global (parallel) operations.
  //
  Array<int> ids;
  Array<Real> locs;
  Array<Real> vels;
  Array<Real> frcs;

  particles.GetParticleIDs(ids);
  particles.GetParticleLocations(locs);
  particles.GetParticleVelocities(vels);
  frcs.resize(locs.size());
  fill(frcs.begin(), frcs.end(), 0.0);

  //
  // Build mapping from particle ID to index into the particle
  // location array.
  //
  // Building this mapping every time we call this routine seems
  // expensive, but any trickery to avoid doing so might be
  // problematic.
  //
  map<int,int> idmap;
  for (int i=0; i<ids.size(); i++)
    idmap.insert(pair<int,int>(ids[i],i));

  //
  // Loop through particles and accumulate forces.
  //
  for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it) {
    const PBox& pbox = pmap_it->second;
    const int   n    = pbox.size();

    for (int i = 0; i < n; i++) {
      const ParticleType& p = pbox[i];

      const int id = p.m_id;
      if (id <= 0) continue;

      ParticleSpringRange range = particle_to_spring_map.equal_range(p.m_id);
      for (ParticleSpringMap::iterator it=range.first; it!=range.second; ++it) {
        if (it->first != id) continue;
        Spring spring = this->springs.at(it->second);

        // get indexes into particle location array
        int i1 = (spring.p1 == id) ? idmap[spring.p1] : idmap[spring.p2];
        int i2 = (spring.p1 == id) ? idmap[spring.p2] : idmap[spring.p1];

        bpm_compute_forces_2d(&locs[BL_SPACEDIM*i1], &vels[BL_SPACEDIM*i1], &frcs[BL_SPACEDIM*i1],
                              &locs[BL_SPACEDIM*i2], &vels[BL_SPACEDIM*i2], &frcs[BL_SPACEDIM*i2],
                              spring.k, spring.l);

      }
    }
  }

  //
  // Loop through particles and spread forces.
  //
  for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it) {
    const int   grid = pmap_it->first;
    const PBox& pbox = pmap_it->second;
    const int   n    = pbox.size();

    FArrayBox& ffab = f[grid];
    FArrayBox& divffab = divf[grid];

    for (int i = 0; i < n; i++) {
      const ParticleType& p = pbox[i];

      const int id = p.m_id;
      if (id <= 0) continue;

      int i1 = idmap[id];

      bpm_parameters params;
      params.type = 1;
      params.diameter = 0.05;

      bpm_spread_forces_2d(&locs[BL_SPACEDIM*i1], &frcs[BL_SPACEDIM*i1], &params,
                           ffab.loVect(), ffab.hiVect(), 4, geom.ProbLo(), geom.ProbHi(), geom.CellSize()[0],
                           ffab.dataPtr(), divffab.dataPtr());

    }
  }

}

//
// Read springs from file.
//
// Each line represents a spring, format is: id, id, resting length, spring constant
//
void BPM::InitFromAsciiFile(const string& file)
{
  // XXX: Use NReaders...

  ifstream ifs;
  ifs.open(file.c_str(), ios::in);

  if (!ifs.good())
    BoxLib::FileOpenFailed(file);

  int p1, p2;
  Real l, k;

  while (ifs >> p1 >> p2 >> l >> k) {
    Spring spring = { p1, p2, l, k };
    this->springs.push_back(spring);
    this->particle_to_spring_map.insert(pair<int,int>(p1, this->springs.size()-1));
    this->particle_to_spring_map.insert(pair<int,int>(p2, this->springs.size()-1));
  }

  cout << "BPM::InitFromAsciiFile: " << this->springs.size() << " springs read." << endl;
}
