#include <AMReX_Vector.H>
#include <set>

#include <depRand.H>

using namespace amrex;

//
// BoxLib Interface to Mersenne Twistor
//

/* A C-program for MT19937: Real number version (1999/10/28)    */
/*   genrand() generates one pseudorandom real number (double)  */
/* which is uniformly distributed on [0,1]-interval, for each   */
/* call. sgenrand(seed) sets initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be       */
/* called once. (seed is any 32-bit integer.)                   */
/* Integer generator is obtained by modifying two lines.        */
/*   Coded by Takuji Nishimura, considering the suggestions by  */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.            */

/* This library is free software under the Artistic license:       */
/* see the file COPYING distributed together with this code.       */
/* For the verification of the code, its output sequence file      */
/* mt19937-1.out is attached (2001/4/2)                           */

/* Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura. */
/* Any feedback is very welcome. For any question, comments,       */
/* see http://www.math.keio.ac.jp/matumoto/emt.html or email       */
/* matumoto@math.keio.ac.jp                                        */

/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */

unsigned long DepRand::mt19937::init_seed;
unsigned long DepRand::mt19937::mt[DepRand::mt19937::N];
int           DepRand::mt19937::mti;

//
// initializing with a NONZERO seed.
//
void
DepRand::mt19937::sgenrand(unsigned long seed)
{
    mt[0]= seed & 0xffffffffUL;
    for ( mti=1; mti<N; mti++ ) 
    {
        mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30L)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;       /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void 
DepRand::mt19937::sgenrand(unsigned long init_key[], int key_length)
{
    int i, j, k;
    sgenrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for ( ; k; k-- ) 
    {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for ( k=N-1; k; k-- ) 
    {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

void
DepRand::mt19937::reload()
{
    unsigned long y;
    int kk;

    const int M = 397;

#define MATRIX_A    0x9908B0DFUL // Constant vector a
#define UPPER_MASK  0x80000000UL // Most significant w-r bits
#define LOWER_MASK  0x7FFFFFFFUL // least significant r bits
    //
    // mag01[x] = x * MATRIX_A  for x=0,1
    //
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    for ( kk=0; kk<N-M; kk++ )
    {
	y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
	mt[kk] = mt[kk+M] ^ (y >> 1L) ^ mag01[y & 0x1UL];
    }
    for ( ; kk<N-1; kk++ )
    {
	y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
	mt[kk] = mt[kk+(M-N)] ^ (y >> 1L) ^ mag01[y & 0x1UL];
    }
    y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
    mt[N-1] = mt[M-1] ^ (y >> 1L) ^ mag01[y & 0x1UL];

    mti = 0;

#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK
}

unsigned long
DepRand::mt19937::igenrand()
{
    //
    // Generate N words at one time.
    //
    if ( mti >= N ) reload();

    unsigned long y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

DepRand::mt19937::mt19937(unsigned long seed)
{
    init_seed = seed;
    mti = N;
    sgenrand(seed);
}

DepRand::mt19937::mt19937(unsigned long seed, int numprocs)
{
#ifdef _OPENMP
#pragma omp parallel
  {
    init_seed = seed + omp_get_thread_num() * numprocs;
    mti = N;
    sgenrand(init_seed);
  }
#else
    init_seed = seed;
    mti = N;
    sgenrand(init_seed);
#endif
}

DepRand::mt19937::mt19937 (unsigned long seed_array[], int len)
{
    sgenrand(seed_array, len);
}

void
DepRand::mt19937::rewind()
{
    sgenrand(init_seed);
}

void
DepRand::mt19937::reset(unsigned long seed)
{
    sgenrand(seed);
}

//
// [0,1] random numbers
//
double
DepRand::mt19937::d_value()
{
    return double(igenrand()) * (1.0/4294967295.0);  // divided by 2^32-1
}

//
// [0,1) random numbers
//
double
DepRand::mt19937::d1_value()
{
    return double(igenrand()) * (1.0/4294967296.0);  // divided by 2^32
}

//
// (0,1) random numbers
//
double
DepRand::mt19937::d2_value()
{
    return (double(igenrand()) + .5) * (1.0/4294967296.0);  // divided by 2^32
}

long
DepRand::mt19937::l_value()
{
    return (long)(igenrand()>>1);
}

unsigned long
DepRand::mt19937::u_value()
{
    return igenrand();
}

void
DepRand::mt19937::save (Vector<unsigned long>& state) const
{
    state.resize(N+2);
    state[0] = init_seed;
    for (int i = 0; i < N; i++)
        state[i+1] = mt[i];
    state[N+1] = mti;
}

int
DepRand::mt19937::RNGstatesize () const
{
    return N+2;
}

void
DepRand::mt19937::restore (const Vector<unsigned long>& state)
{
    if (state.size() != N+2)
        Error("mt19937::restore(): incorrectly sized state vector");

    init_seed = state[0];
    for (int i = 0; i < N; i++)
        mt[i] = state[i+1];
    mti = state[N+1];

    if (mti < 0 || mti > N)
        Error("mt19937::restore(): mti out-of-bounds");
}

namespace
{
    DepRand::mt19937 the_generator;
}

void
DepRand::InitRandom (unsigned long seed)
{
    the_generator = mt19937(seed);
}

void
DepRand::InitRandom (unsigned long seed, int numprocs)
{
    the_generator = mt19937(seed, numprocs);
}

void DepRand::ResetRandomSeed(unsigned long seed)
{
    the_generator.reset(seed);
}

double
DepRand::Random ()
{
    return the_generator.d_value();
}

double
DepRand::Random1 ()
{
    return the_generator.d1_value();
}

double
DepRand::Random2 ()
{
    return the_generator.d2_value();
}

unsigned long
DepRand::Random_int(unsigned long n)
{
  const unsigned long umax = 4294967295UL; // 2^32-1
  BL_ASSERT( n > 0 && n <= umax ); 
  unsigned long scale = umax/n;
  unsigned long r;
  do {
    r = the_generator.u_value() / scale;
  } while (r >= n);
  return r;
}

void
DepRand::SaveRandomState (Vector<unsigned long>& state)
{
    the_generator.save(state);
}

int
DepRand::sizeofRandomState ()
{
    return the_generator.RNGstatesize();
}

void
DepRand::RestoreRandomState (const Vector<unsigned long>& state)
{
    the_generator.restore(state);
}

void
DepRand::UniqueRandomSubset (Vector<int> &uSet, int setSize, int poolSize,
                            bool printSet)
{
  if(setSize > poolSize) {
    Abort("**** Error in UniqueRandomSubset:  setSize > poolSize.");
  }
  std::set<int> copySet;
  Vector<int> uSetTemp;
  while(copySet.size() < (unsigned long)setSize) {
    int r(DepRand::Random_int(poolSize));
    if(copySet.find(r) == copySet.end()) {
      copySet.insert(r);
      uSetTemp.push_back(r);
    }
  }
  uSet = uSetTemp;
  if(printSet) {
    for(int i(0); i < uSet.size(); ++i) {
      std::cout << "uSet[" << i << "]  = " << uSet[i] << std::endl;
    }
  }
}


#if 0
//
// Fortran entry points for DepRand::Random().
//

BL_FORT_PROC_DECL(BLUTILINITRAND,blutilinitrand)(const int* sd)
{
    unsigned long seed = *sd;
    DepRand::InitRandom(seed);
}

BL_FORT_PROC_DECL(BLINITRAND,blinitrand)(const int* sd)
{
    unsigned long seed = *sd;
    DepRand::InitRandom(seed);
}

BL_FORT_PROC_DECL(BLUTILRAND,blutilrand)(Real* rn)
{
    *rn = DepRand::Random();
}
#endif
