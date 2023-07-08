#ifndef IAMR_FLUXBOXES_H_
#define IAMR_FLUXBOXES_H_

#include <AMReX_AmrLevel.H>

class FluxBoxes
{
public:

    FluxBoxes () = default;

    explicit FluxBoxes (const amrex::AmrLevel* amr_level, int nvar=1, int nghost=0)
        { define(amr_level, nvar, nghost); }

    ~FluxBoxes () { clear(); }

    FluxBoxes (FluxBoxes const&) = delete;
    FluxBoxes (FluxBoxes &&) = delete;
    FluxBoxes& operator= (FluxBoxes const&) = delete;
    FluxBoxes& operator= (FluxBoxes &&) = delete;

    amrex::MultiFab** define (const amrex::AmrLevel* amr_level, int nvar=1, int nghost=0);

    void clear ();

    amrex::MultiFab** get () {return data;}

private:

    amrex::MultiFab** data = nullptr;

};

#endif
