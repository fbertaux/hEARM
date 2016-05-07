#include "hearmWrapper.h"
#include "common.h"

void
HearmWrapper::initialize (int seed)
{
    if (isInit) { cout << " DONT INIT TWICE THE HEARM WRAPPER !!!!!!!! I quit" << endl ; exit(1) ; }

    initializeParams () ;
    initializeSimulators (seed) ;
}

void
HearmWrapper::initializeParams ()
{
    GTP_Params = new GeneToProtParams () ;
    TrailApopParams = new TrailApoptosisParams () ;
}


void
HearmWrapper::initializeSimulators (int seed)
{
    // create needed simulators
    GeneMrna_Sim = new GeneMrnaSim ( GTP_Params , seed ) ;
//    cout << "gene to mrna simulator created and set." << endl ;
    HybridOde_Sim = new HybridOdeSim ( TrailApopParams , GTP_Params ) ;
//    cout << "hybrid ode simulater created and set." << endl ;
    isInit = true ;
}


PreStimCell*
HearmWrapper::sampleOnePreTrailCell ()
{
    PreStimCell* pre_cell = new PreStimCell (GTP_Params,GeneMrna_Sim) ;
    pre_cell->simulate ( 10 * 24 * 3600 ) ; // 10 days largely enough to reach steady-state distribution
    return pre_cell ;
}

PostStimCell*
HearmWrapper::simulateTrailPulseOnCell (double duration_in_minutes, PreStimCell *initial_cell, double trail_ng_per_mL)
{
    PostStimCell* cell = new PostStimCell (TrailApopParams,HybridOde_Sim,"cPARP",100000.,initial_cell ) ;
    simulateTrailPulseOnCell (duration_in_minutes,cell,trail_ng_per_mL) ;
    return cell ;
}


void
HearmWrapper::simulateTrailPulseOnCell (double duration_in_minutes, PostStimCell *cell, double trail_ng_per_mL)
{
    cell->addTrail ( trail_ng_per_mL*60 ) ;
    double advised_time_step = 180. ;
    double time_step = ( advised_time_step < duration_in_minutes ) ? advised_time_step : duration_in_minutes ;
    cell->simulate ( duration_in_minutes*60 , advised_time_step*60 ) ;
}
