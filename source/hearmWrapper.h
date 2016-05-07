#ifndef HEARMWRAPPER_H
#define HEARMWRAPPER_H


#include "cells/prestimcell.h"
#include "cells/poststimcell.h"

struct HearmWrapper
{
    HearmWrapper () : isInit(false) {}

    // regarding initialization
    bool isInit ;
    void initialize (int seed=0) ;
    void initializeParams () ;
    void initializeSimulators (int seed) ;

    // pointers to params
    GeneToProtParams* GTP_Params ;
    TrailApoptosisParams* TrailApopParams ;

    // pointers to simulators
    GeneMrnaSim* GeneMrna_Sim ;
    HybridOdeSim* HybridOde_Sim ;

    // methods for performing simulations
    PreStimCell* sampleOnePreTrailCell () ;
    PostStimCell* simulateTrailPulseOnCell ( double duration_in_minutes , PreStimCell* initial_cell , double trail_ng_per_mL ) ;
    void simulateTrailPulseOnCell ( double duration_in_minutes , PostStimCell* initial_cell , double trail_ng_per_mL ) ;

};


#endif // HEARMWRAPPER_H
