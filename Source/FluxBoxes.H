#ifndef _FLUXBOXES_H_
#define _FLUXBOXES_H_

#include <AMReX_AmrLevel.H>

class FluxBoxes
{
public:

    FluxBoxes () : data(0) {}

    explicit FluxBoxes (const amrex::AmrLevel* amr_level, int nvar=1, int nghost=0)
	: data(0) { define(amr_level, nvar, nghost); }

    ~FluxBoxes () { clear(); }

    amrex::MultiFab** define (const amrex::AmrLevel* amr_level, int nvar=1, int nghost=0);

    void clear ();

    amrex::MultiFab** get () {return data;}

private:

    amrex::MultiFab** data;

};

#endif
