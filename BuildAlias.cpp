//
// $Id: BuildAlias.cpp,v 1.3 1997-09-24 17:45:44 lijewski Exp $
//

#include <BuildAlias.H>

static void error(char* message)
{
  BoxLib::Error(message);
}

static int compatible(const level_mesh& mesh, const BoxArray& bs)
{
  int retval = 1;
  if (mesh.ngrids() != bs.length()) {
    retval = 0;
  }
  else {
    for (int i = 0; i < mesh.ngrids(); i++) {
      Box tmp = bs[i];
      if (tmp.convert(IntVect::TheCellVector()) != mesh[i]) {
	  cerr << "compatable: NOT: tmp = " << tmp
	       << " mesh[i] = " << mesh[i] << endl;
	  retval = 0;
      }
      
    }
  }
  return retval;
}

level_mesh make_level_mesh(const Box& domain, const BoxArray& bs)
{
  if (!bs.ready())
    BoxLib::Error("make_level_mesh---bad BoxArray");

  level_mesh lmesh(domain, bs.length());
  for (int i = 0; i < bs.length(); i++) {
    lmesh.set_box(i, bs[i]);
  }
  lmesh.check();
  return lmesh;
}

amr_mesh make_amr_mesh(Box domain, int nlev, const int rat[],
		       const BoxArray *bs[])
{
  if (nlev <= 0)
    BoxLib::Error("make_amr_mesh---need nlev > 0");

  amr_mesh amesh(nlev);
  amesh.set_level(0, make_level_mesh(domain, *bs[0]));
  for (int i = 1; i < nlev; i++) {
    domain.refine(rat[i-1]);
    amesh.set_level(i, make_level_mesh(domain, *bs[i]));
  }
  amesh.check();
  return amesh;
}

level_real make_level_real(const level_mesh& mesh,
			   MultiFab& mf, int icomp)
{
  if (mesh.null() || !compatible(mesh, mf.boxArray()))
    BoxLib::Error("make_level_real---inputs don't match");

  level_real lreal(mesh, mf.boxArray()[0].type(), mf.nGrow(), 0);
  for (int i = 0; i < mesh.ngrids(); i++) {
    lreal.set_grid(i, grid_real(mf[i]), icomp);
  }
  lreal.check();
  return lreal;
}

amr_real make_amr_real(const amr_mesh& mesh, MultiFab mf[])
{
  if (mesh.null() || mesh.nlevels() == 0)
    BoxLib::Error("make_amr_real---need nlev > 0");

  level_real lreal = make_level_real(mesh[0], mf[0]);
  amr_real areal(mesh, lreal.type(), lreal.border(), 0);
  areal.set_level(0, lreal);
  for (int i = 1; i < mesh.nlevels(); i++) {
    areal.set_level(i, make_level_real(mesh[i], mf[i]));
  }
  areal.check();
  return areal;
}

amr_real make_amr_real(const amr_mesh& mesh, MultiFab& mf)
{
  if (mesh.null() || mesh.nlevels() == 0)
    BoxLib::Error("make_amr_real---need nlev > 0");

  level_real lreal = make_level_real(mesh[0], mf);
  amr_real areal(mesh, lreal.type(), lreal.border(), 0);
  areal.set_level(0, lreal);
  areal.check();
  return areal;
}

// ###################################################################
// ##### BLD_AMR_REAL
// ###################################################################
void bldAmrReal(amr_real& areal, amr_mesh mesh, MultiFab& mf,
                int level, int comp)
{
    level_real lreal = make_level_real(mesh[level],mf,comp);
    areal.alloc(mesh,lreal.type(),lreal.border(),0);
    areal.set_level(level,lreal);
}

// ###################################################################
// ##### BLD_AMR_REAL
// ###################################################################
void bldAmrReal(amr_real& areal, amr_mesh mesh, MultiFab* mf[],
                int levc, int levf, int comp)
{
    level_real lreal = make_level_real(mesh[levc],*mf[levc],comp);
    areal.alloc(mesh,lreal.type(),lreal.border(),0);
    areal.set_level(levc,lreal);
    for (int lev = levc+1; lev <= levf; lev++) {
	level_real lreal = make_level_real(mesh[lev],*mf[lev],comp);
	areal.set_level(lev,lreal);  
    }
}

// ###################################################################
// ##### BLD_AMR_REAL
// ###################################################################
void bldAmrReal(amr_real& areal, amr_mesh mesh, MultiFab * mf_crse,
                MultiFab * mf_fine, int levc, int levf, int comp)
{

    assert(levf==levc+1);
    {
    level_real lreal = make_level_real(mesh[levc],*mf_crse,comp);
    areal.alloc(mesh,lreal.type(),lreal.border(),0);
    areal.set_level(levc,lreal);
    }

    {
    level_real lreal = make_level_real(mesh[levf],*mf_fine,comp);
    areal.set_level(levf,lreal);
    }
}

