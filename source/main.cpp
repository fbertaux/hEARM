// EASY CELL, Francois BERTAUX ///


#include "hearmWrapper.h"


int
main ()
{
    cout << endl << "****** hEARM simulation engine ******" << endl;

    // random generator seed
    int seed = 0. ;

    // create and initialize the main object class to deal with hEARM model simulations
    HearmWrapper hearm ;
    hearm.initialize ( seed ) ;

    // parameters
    double duration_minutes =  60. * 8. ;
    double mrna_storage_duration_minutes = 10. ;
    double trail_dose_ng_per_mL = 250 ;
    int num_cells = 500 ;

    // simulate TRAIL treatment on several cells
    int survived = 0 ;
    for ( int n=0 ; n<num_cells ; n++ )
    {
        // create an initial cell state
        PreStimCell* pre_cell = hearm.sampleOnePreTrailCell () ;
        PostStimCell* cell = new PostStimCell (hearm.TrailApopParams,hearm.HybridOde_Sim,"cPARP",100000.,pre_cell) ;

        // simulate treatment
        cell->addTrail ( trail_dose_ng_per_mL * 60. ) ;
        cell->simulate ( duration_minutes*60. , mrna_storage_duration_minutes*60. ) ;
        if (!cell->DidMOMP) survived++ ;
        
        // free memory
        delete pre_cell ;
        delete cell ;
    }

    cout << endl << "Fraction surviving cells = " << (double) survived / (double) num_cells << endl;

    cout << endl << "****** hEARM says goodbye ******" << endl;
    return 0 ;
}
