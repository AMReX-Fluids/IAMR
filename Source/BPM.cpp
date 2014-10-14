
#include <BPM.H>

#include <iostream>

using namespace std;

typedef ParticleContainer<2*BL_SPACEDIM>::PBox PBox;
typedef ParticleContainer<2*BL_SPACEDIM>::PMap PMap;
typedef ParticleContainer<2*BL_SPACEDIM>::ParticleType ParticleType;

//
// Compute spring forces given particle locations.
//
void BPM::ComputeForces(ParticleContainer<2*BL_SPACEDIM>& particles, MultiFab& force, int lev, const int fcomp)
{
  BL_PROFILE("BPM::ComputeForces()");

  const PMap& pmap = particles.GetParticles(lev);

  //
  // Get all particle IDs and locations.  This is a global (parallel)
  // operation.
  //
  Array<int> ids; particles.GetParticleIDs(ids);
  Array<Real> locs; particles.GetParticleLocations(locs);

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
  // Loop through particles and compute forces...
  //
  for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it) {
    const int   grid = pmap_it->first;
    const PBox& pbox = pmap_it->second;
    const int   n    = pbox.size();

    FArrayBox& ffab = force[grid];

    for (int i = 0; i < n; i++) {
      const ParticleType& p = pbox[i];

      const int id = p.m_id;
      if (id <= 0) continue;

      ParticleSpringRange range = this->pmap.equal_range(p.m_id);
      for (ParticleSpringMap::iterator it=range.first; it!=range.second; ++it) {
        if (it->first != id) continue;
        Spring spring = this->springs.at(it->second);

        // get indexes into particle location array
        int i1 = (spring.p1 == id) ? idmap[spring.p1] : idmap[spring.p2];
        int i2 = (spring.p1 == id) ? idmap[spring.p2] : idmap[spring.p1];

        // compute spring force; we're computing the force on particle i1 from particle i2
        Real dvect[BL_SPACEDIM];
        for (int d=0; d < BL_SPACEDIM; d++) {
          dvect[d] = locs[BL_SPACEDIM*i2+d] - locs[BL_SPACEDIM*i1+d];
        }

        Real dnorm = sqrt(D_TERM(dvect[0]*dvect[0], + dvect[1]*dvect[1], + dvect[2]*dvect[2]));
        Real force = spring.k * (spring.l - dnorm);

        Real fvect[BL_SPACEDIM];
        for (int d=0; d < BL_SPACEDIM; d++) {
          fvect[d] = dvect[d] / dnorm * force;
        }

        if (this->verbose > 1) {
          cout << "BPM::ComputeForces: force on " << ids[i1] << " from " << ids[i2] << " is: ";
          for (int d=0; d < BL_SPACEDIM; d++) {
            cout << fvect[d] << " ";
          }
          cout << endl;
        }

        // spread force into fab...


      }
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
    this->pmap.insert(pair<int,int>(p1, this->springs.size()-1));
    this->pmap.insert(pair<int,int>(p2, this->springs.size()-1));
  }

  cout << "BPM::InitFromAsciiFile: " << this->springs.size() << " springs read." << endl;
}
