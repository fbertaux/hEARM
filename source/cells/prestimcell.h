#ifndef PRESTIMCELL_H
#define PRESTIMCELL_H

#include "genemrnasim.h"
#include "genetoprotparams.h"
#include "vector"

class PreStimCell
{
public:

    static Doub AbsoluteTime;

    // state fields
    VecDoub GeneMrnaState;
    VecDoub ProtState;
    Doub CellAge;

    // knowing parameters
    GeneToProtParams* GeneToProtParameters ;

    // knowing a simulator
    GeneMrnaSim* GeneMrna_Sim ;

    // CONSTRUCTORS
    PreStimCell (GeneToProtParams* gtp_pars,GeneMrnaSim* gm_sim) ;
    PreStimCell (PreStimCell &toCopy); // for division ( copy has CellAge 0, mother CellAge not copied)


    // CELL SIMULATION METHODS
    void simulate (Doub duration) ;

};


#endif // PRESTIMCELL_H
