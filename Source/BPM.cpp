
#include <cmath>
#include <iostream>

#include <BPM.H>

using namespace std;

typedef ParticleContainer<2*BL_SPACEDIM>::PBox PBox;
typedef ParticleContainer<2*BL_SPACEDIM>::PMap PMap;
typedef ParticleContainer<2*BL_SPACEDIM>::ParticleType ParticleType;

extern "C" {
  void bpm_compute_forces_2d(double const *x1, double const *x2,
                             double sconstant, double slength,
                             double blob_diam, int blob_patch_size,
                             int const *lo, int const *hi, double const *xlo, double const *xhi, double dx,
                             double *f, double *divf);
}

//
// Compute spring forces given particle locations.
//
void BPM::ComputeForces(ParticleContainer<2*BL_SPACEDIM>& particles, MultiFab& force, Geometry& geom, int lev, int fcomp)
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

  MultiFab divf(force.boxArray(), 1, 2);

  //
  // Loop through particles and compute forces...
  //
  for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it) {
    const int   grid = pmap_it->first;
    const PBox& pbox = pmap_it->second;
    const int   n    = pbox.size();

    FArrayBox& ffab = force[grid];
    FArrayBox& divffab = divf[grid];

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

        double blob_diameter = 123.0;
        int blob_patch_size = 2;

        // as it is now, this doesn't accumulate forces and then apply
        // the bpm, it applies the bpm to each end of each spring...

        bpm_compute_forces_2d(&locs[BL_SPACEDIM*i1], &locs[BL_SPACEDIM*i2],
                              spring.k, spring.l, blob_diameter, blob_patch_size,
                              ffab.loVect(), ffab.hiVect(), geom.ProbLo(), geom.ProbHi(), geom.CellSize()[0],
                              ffab.dataPtr(), divffab.dataPtr());

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
