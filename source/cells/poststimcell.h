#ifndef POSTSTIMCELL_H
#define POSTSTIMCELL_H

#define CONTENT_OUTPUT 0

#include "prestimcell.h"
#include "hybridodesim.h"
#include "trailapoptosisparams.h"

#define TRAJ_OUTPUT 0

class PostStimCell
{
public:

    //! basic constructor
    PostStimCell ( GeneToProtParams* gtp_pars,GeneMrnaSim* gm_sim,
                    TrailApoptosisParams* trail_pars,HybridOdeSim* hyb_sim,
                    string death_crit, double death_tresh_fact ) ;

    //! to init from a PreStimCell
    PostStimCell (TrailApoptosisParams* trail_pars,HybridOdeSim* hyb_sim,
                   string death_crit, double death_tresh_fact,
                   PreStimCell* psCell) ;

    //! copy constructor, mimicking division ( new cell has CellAge = 0 )
    PostStimCell(PostStimCell &toCopy) ;


    // knowing parameters
    GeneToProtParams* GeneToProtParameters ;
    TrailApoptosisParams* TrailApopParameters ;

    // knowing simulators
    GeneMrnaSim* GeneMrna_Sim ;
    HybridOdeSim* HybridOde_Sim ;

    // knowing when to die
    string DeathCriterion ;
    Doub DeathThreshFactor ;

    // state fields : molecular content
    VecDoub GeneMrnaState;
    VecDoub AllProtState;

    // state fields : phenotype
    Doub AbsoluteTime;
    Bool DidMOMP;
    Doub MOMPTime;

    // technical field
    Bool IsReadyForFutureTrajectoryComputation;


    // simulation methods
    void addTrail (Doub Trail_Dose) ;
    void computeFutureGeneMrnaTrajectory(Doub duration);
    void simulate(Doub duration, Doub timestep);
    void checkMompOrDeath () ;

    // "reset" methods, to do some nice tests
    void resetTranscriptionalState () ;
    void resetActiveForms () ;
    void resetNativeForms () ;

#if CONTENT_OUTPUT == 1
    /*to compute some outputs*/
    VecDoub fon;
    Bool storedStateOfInterest;
    Bool outputDone;
    VecDoub stateOfInterest;
    /**/
#endif

#if TRAJ_OUTPUT == 1
	MatDoub ProtTraj;
	VecDoub TimePoints;
	static const Int NumTimePointsMax = 10000;
	Int NumTimePointsDone;
#endif


};

#endif // POSTSTIMCELL_H
