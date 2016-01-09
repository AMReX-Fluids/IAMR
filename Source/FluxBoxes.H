#ifndef _FLUXBOXES_H_
#define _FLUXBOXES_H_

#include <AmrLevel.H>

class FluxBoxes
{
public:

    FluxBoxes () : data(0) {}

    explicit FluxBoxes (const AmrLevel* amr_level, int nvar=1, int nghost=0)
	: data(define(amr_level, nvar, nghost)) {}

    ~FluxBoxes () { clear(); };

    MultiFab** define (const AmrLevel* amr_level, int nvar=1, int nghost=0);

    void clear ();

    MultiFab** get () {return data;}

private:

    MultiFab** data;

};

#endif
