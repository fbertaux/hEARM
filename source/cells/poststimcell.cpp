#include "poststimcell.h"

#define DEBUG 0



//! CONSTRUCTOR FROM NOTHING
PostStimCell::PostStimCell (GeneToProtParams* gtp_pars,GeneMrnaSim* gm_sim,
                              TrailApoptosisParams* trail_pars,HybridOdeSim* hyb_sim,
                              string death_crit, double death_tresh_fact)
    : GeneToProtParameters (gtp_pars) , GeneMrna_Sim (gm_sim) ,
      TrailApopParameters (trail_pars) , HybridOde_Sim ( hyb_sim) ,
      DeathCriterion (death_crit) , DeathThreshFactor ( death_tresh_fact) ,
      IsReadyForFutureTrajectoryComputation(true), AbsoluteTime(0.), DidMOMP(false), MOMPTime(0.)
{
    // init state vectors
    GeneMrnaState = VecDoub(3*GeneToProtParameters->NumGenes,0.);
    AllProtState = VecDoub(TrailApopParameters->num_species+1,0.);
    for (Int i=0;i<GeneToProtParameters->NumGenes;i++)
    {
        GeneMrnaState[3*i] = 1; GeneMrnaState[3*i+1] = 0; // gene is on for all genes
        GeneMrnaState[3*i+2] = (Int) GeneToProtParameters->giveMrnaMean(i); // Mrna to mean level
    }

    // iterate on native forms and set the mean
    typedef map<string, int>::iterator it_type ;
    for ( it_type iterator = TrailApopParameters->nativeNamesToIndex.begin() ; iterator != TrailApopParameters->nativeNamesToIndex.end() ; iterator++ )
    {
        Int genIdx = iterator->second ;
        Int natIdx = TrailApopParameters->nativeNamesToNativeIndex [ iterator->first ] ;
        AllProtState[ genIdx ] = GeneToProtParameters->giveProtMean (natIdx) ;
    }


#if CONTENT_OUTPUT == 1
    /*to compute some outputs*/
    fon = VecDoub(17,0.);
    storedStateOfInterest = false;
    outputDone = false;
    stateOfInterest = VecDoub(TrailApopParameters->num_species,0.);
    /**/
#endif

#if TRAJ_OUTPUT == 1
    ProtTraj = MatDoub(TrailApopParameters->num_species,NumTimePointsMax,0.);
    TimePoints = VecDoub(NumTimePointsMax,0.);
    NumTimePointsDone = 0;
#endif

}

//! CONSTRUCTOR FROM A PRE STIM CELL
PostStimCell::PostStimCell (TrailApoptosisParams* trail_pars,HybridOdeSim* hyb_sim,
                              string death_crit, double death_tresh_fact,
                              PreStimCell *psCell)
    : GeneToProtParameters (psCell->GeneToProtParameters) , GeneMrna_Sim (psCell->GeneMrna_Sim) ,
      TrailApopParameters (trail_pars) , HybridOde_Sim ( hyb_sim) ,
      DeathCriterion (death_crit) , DeathThreshFactor ( death_tresh_fact) ,
      IsReadyForFutureTrajectoryComputation(true), AbsoluteTime(0.), DidMOMP(false), MOMPTime(0.)
{
//    cout << "creating a post stim cell from a prestim cell" << endl ;
    // init state vectors
    GeneMrnaState = VecDoub(3*GeneToProtParameters->NumGenes,0.);
    AllProtState = VecDoub(TrailApopParameters->num_species+1,0.);
    for (Int i=0;i<GeneToProtParameters->NumGenes;i++)
    {
        GeneMrnaState[3*i] = psCell->GeneMrnaState[3*i]; GeneMrnaState[3*i+1] = psCell->GeneMrnaState[3*i+1];
        GeneMrnaState[3*i+2] = psCell->GeneMrnaState[3*i+2];
        AllProtState[TrailApopParameters->nativeIndexToIndex[i]] = psCell->ProtState[i];
    }
//    cout << "state initialized" << endl ;

#if CONTENT_OUTPUT == 1
    /*to compute some outputs*/
    fon = VecDoub(17,0.);
    storedStateOfInterest = false;
    outputDone = false;
    stateOfInterest = VecDoub(TrailApopParameters->num_species,0.);
    /**/
#endif

#if TRAJ_OUTPUT == 1
    ProtTraj = MatDoub(TrailApopParameters->num_species,NumTimePointsMax,0.);
    TimePoints = VecDoub(NumTimePointsMax,0.);
    NumTimePointsDone = 0;
#endif
}

//! COPY CONSTRUCTOR
PostStimCell::PostStimCell(PostStimCell &toCopy)
    : GeneToProtParameters (toCopy.GeneToProtParameters) , GeneMrna_Sim (toCopy.GeneMrna_Sim) ,
      TrailApopParameters (toCopy.TrailApopParameters) , HybridOde_Sim (toCopy.HybridOde_Sim) ,
      DeathCriterion (toCopy.DeathCriterion) , DeathThreshFactor (toCopy.DeathThreshFactor) ,
      IsReadyForFutureTrajectoryComputation(true), AbsoluteTime(toCopy.AbsoluteTime), DidMOMP(false), MOMPTime(0.)
{
    if (toCopy.DidMOMP) { cout << "Don't divide a MOMP cell.. I quit." << endl; exit(4); }
    // init state vectors
    GeneMrnaState = VecDoub(3*GeneToProtParameters->NumGenes,0.);
    AllProtState = VecDoub(TrailApopParameters->num_species+1,0.);
    for (Int i=0;i<GeneToProtParameters->NumGenes;i++)
    {
        GeneMrnaState[3*i] = toCopy.GeneMrnaState[3*i]; GeneMrnaState[3*i+1] = toCopy.GeneMrnaState[3*i+1];
        GeneMrnaState[3*i+2] = toCopy.GeneMrnaState[3*i+2];
    }
    for (Int j=0;j<TrailApopParameters->num_species;j++)
    {
        AllProtState[j] = toCopy.AllProtState[j];
    }

#if CONTENT_OUTPUT == 1
    /*to compute some outputs*/
    fon = VecDoub(17,0.);
    storedStateOfInterest = false;
    outputDone = false;
    stateOfInterest = VecDoub(TrailApopParameters->num_species,0.);
    /**/
#endif

#if TRAJ_OUTPUT == 1
    ProtTraj = MatDoub(TrailApopParameters->num_species,NumTimePointsMax,0.);
    TimePoints = VecDoub(NumTimePointsMax,0.);
    NumTimePointsDone = 0;
#endif
}

// simulation methods
void
PostStimCell::computeFutureGeneMrnaTrajectory(Doub duration)
{


    if (!IsReadyForFutureTrajectoryComputation) { cout << "You didn't 'cleared' previous future trajectory computation. I quit." << endl; exit(10); }
    // prepare mrna sim for this cell
    GeneMrna_Sim->prepareForSteps(GeneMrnaState);
//    cout << "in mrna traj : gene mrna sim prepared"<< endl;
    Doub dt,nt;
    Doub t = 0;
    HybridOde_Sim->EventObtained = 1;
    // store init state
    //cout << "in mrna traj : just before init state loop"<< endl;
    for (int i=0;i<GeneToProtParameters->NumGenes;i++)
    {
        //cout << "in mrna traj : in init state loop"<< endl;
        HybridOde_Sim->FutureGeneTable[i][0] = GeneMrnaState[3*i];
        HybridOde_Sim->FutureMrnaTable[i][0] = GeneMrnaState[3*i+2];
        HybridOde_Sim->FutureGeneMrnaEventsTable[0] = 0.;
    }
    // simulate
//    cout << "in mrna traj : just before simulating " << endl;
//    cout << "should simulate during " << duration << endl;
    //        cout << " R prot  : " << AllProtState[1] << endl;
//    cout << " R mrna  : " << GeneMrnaState[2] << endl;
    //    cout << " R gene  on : " << GeneMrnaState[0] << endl;
    //    cout << " R gene  off : " << GeneMrnaState[1] << endl;
    while (t<duration)
    {
        nt = GeneMrna_Sim->step(GeneMrnaState,duration);
        dt = nt - t;
//        cout << "event dt : " << dt << endl;
        if (dt==0) { cout << "dt = 0.." << endl; exit(3);}
//        cout << "num events : " << HybridOdeSim->EventObtained << endl;
        t = GeneMrna_Sim->t;
        // store
        for (int i=0;i<GeneToProtParameters->NumGenes;i++)
        {
            HybridOde_Sim->FutureGeneTable[i][HybridOde_Sim->EventObtained] = GeneMrnaState[3*i];
            HybridOde_Sim->FutureMrnaTable[i][HybridOde_Sim->EventObtained] = GeneMrnaState[3*i+2];
            HybridOde_Sim->FutureGeneMrnaEventsTable[HybridOde_Sim->EventObtained] = t;
        }
        HybridOde_Sim->EventObtained++;
        if (HybridOde_Sim->EventObtained>HybridOde_Sim->MaxEvents) { cout << "Mrna Table To short. Change Max Events ?" << endl; exit(42); }
    }
    IsReadyForFutureTrajectoryComputation = false;
//    cout << " mrna traj computed" << endl;

}

void
PostStimCell::simulate(Doub duration, Doub timestepOde)
{
    // cout << "timestep ODE = " << timestepOde << endl ;
    if (DidMOMP) { return; }

    computeFutureGeneMrnaTrajectory(duration) ;
    HybridOde_Sim->stochRhs->setMrnaTables(&HybridOde_Sim->FutureMrnaTable,&HybridOde_Sim->FutureGeneMrnaEventsTable,HybridOde_Sim->EventObtained);

    // debug
    // Doub AbsBeforeCall = AbsoluteTime;

    Doub t = 0;
    while (t<duration)
    {
        // cout << "step." << endl ;
        if (t+timestepOde < duration)
        {
            if (!DidMOMP)
            {
                HybridOde_Sim->OdeIntStoch->integrate(AllProtState,t,t+timestepOde);
                checkMompOrDeath ();
                if (DidMOMP) { MOMPTime = AbsoluteTime; return;}
            }
            t += timestepOde;
            AbsoluteTime += timestepOde;
        }
        else
        {
            if (!DidMOMP)
            {
                HybridOde_Sim->OdeIntStoch->integrate(AllProtState,t,duration);
                checkMompOrDeath ();
                if (DidMOMP) { MOMPTime = AbsoluteTime; return;}
            }
            AbsoluteTime += duration - t;
            t = duration;
            //            if (AbsoluteTime!=AbsBeforeCall+duration) { cout << "smthg wrong 1" << endl; exit(34); }
        }

    }
    // cout << "-------------------" << endl << endl ;

    int indexInTableToTake = HybridOde_Sim->EventObtained - 1 ;

    for (Int i=0;i<GeneToProtParameters->NumGenes;i++)
    {
        GeneMrnaState[3*i] = HybridOde_Sim->FutureGeneTable[i][indexInTableToTake];
        GeneMrnaState[3*i+1] = 1 - GeneMrnaState[3*i];
        GeneMrnaState[3*i+2] = HybridOde_Sim->FutureMrnaTable[i][indexInTableToTake];
    }
    IsReadyForFutureTrajectoryComputation = true;
}

void
PostStimCell::addTrail(Doub Trail_Dose)
{
    AllProtState[0] += Trail_Dose;
}


void
PostStimCell::checkMompOrDeath ()
{
    if (DidMOMP) { cout << "DONT CHECK DEATH ON A DEAD CELL !" << endl; exit(121); }

    // 28 = smac_m ; 53 = smac_m : M*
    // 29 = smac_r ; 30 = smac ; 57 = smac : XIAP

    if(DeathCriterion=="MOMP")
    {
        DidMOMP =  AllProtState[28]+AllProtState[53] < DeathThreshFactor*(AllProtState[29]+AllProtState[30]+AllProtState[57]);
        if (DidMOMP&&DeathThreshFactor==0.)
        {
            cout << "\t\tinside : " << AllProtState[28]+AllProtState[53]  << endl;
            cout << "\t  Smac_m  " << AllProtState[28] << endl;
            cout << "\t  Smac_m:M*  " << AllProtState[53] << endl;
            cout << "\t\toutside : " << AllProtState[29]+AllProtState[30]+AllProtState[57] << endl;
            cout << "\t  smac_r  " << AllProtState[29] << endl;
            cout << "\t  smac  " << AllProtState[30] << endl;
            cout << "\t  smac:XIAP " << AllProtState[57] << endl;
        }
    }
    else if(DeathCriterion=="cPARP")
    {
        DidMOMP = ( AllProtState[13] > DeathThreshFactor ) ;
    }

    else if (DeathCriterion=="PARP")
    {
        DidMOMP = ( AllProtState[12] < DeathThreshFactor ) ;
//        cout << " PARP DEATH " << endl ;
    }
    else if ( DeathCriterion == "RatioCleavedPARP" )
    {
        DidMOMP = (  ( AllProtState[13] / AllProtState[12] ) > DeathThreshFactor ) ;
    }
    else if (DeathCriterion=="None")
    {

    }
    else
    {
        cout << "Unknown death crit ?" << endl ;
        exit (3) ;
    }
    //    if (true/*DidMOMP*/)
    //    {
    //        //        cout << "\t\tinside : " << AllProtState[28]+AllProtState[53]  << endl;
    //        //        cout << "\t  Smac_m  " << AllProtState[28] << endl;
    //        //        cout << "\t  Smac_m:M*  " << AllProtState[53] << endl;
    //        //        cout << "\t\toutside : " << AllProtState[29]+AllProtState[30]+AllProtState[57] << endl;
    //        //        cout << "\t  smac_r  " << AllProtState[29] << endl;
    //        //        cout << "\t  smac  " << AllProtState[30] << endl;
    //        //        cout << "\t  smac:XIAP " << AllProtState[57] << endl;
    //        //        cout << "trail : " << AllProtState[0] << endl;
    //        //        cout << "generation of creation : " << GenerationOfCreation << endl;
    //        // exit (43);
    //    }
}


void
PostStimCell::resetTranscriptionalState ()
{
    // create a preStimCell , simulate long enough
    PreStimCell* psCell = new PreStimCell (GeneToProtParameters,GeneMrna_Sim) ;
    double timeToReachEq = 20 * 24 * 60 * 60 ; // 20 days in seconds
    psCell->simulate (timeToReachEq) ;

    // copy the transcriptional state into current post stim cell
    for (Int i=0;i<GeneToProtParameters->NumGenes;i++)
    {
        GeneMrnaState[3*i] = psCell->GeneMrnaState[3*i]; GeneMrnaState[3*i+1] = psCell->GeneMrnaState[3*i+1];
        GeneMrnaState[3*i+2] = psCell->GeneMrnaState[3*i+2];
    }

    // delete the psCell
    delete psCell ;
}

void
PostStimCell::resetNativeForms ()
{
    // create a preStimCell , simulate long enough
    PreStimCell* psCell = new PreStimCell (GeneToProtParameters,GeneMrna_Sim) ;
    double timeToReachEq = 20 * 24 * 60 * 60 ; // 20 days in seconds
    psCell->simulate (timeToReachEq) ;

    // copy the (native) protein state into current post stim cell
    for (Int i=0;i<GeneToProtParameters->NumGenes;i++)
    {
        AllProtState[TrailApopParameters->nativeIndexToIndex[i]] = psCell->ProtState[i] ;
    }

    // delete the psCell
    delete psCell ;
}

void
PostStimCell::resetActiveForms ()
{
    // iterate on the species except TRAIL
    for (Int j=1;j<TrailApopParameters->num_species;j++)
    {
        // reset only if not native
        if ( !TrailApopParameters->isNativeFromIndex (j) ) AllProtState[j] = 0. ;
    }
}

