#include "hybridodesim.h"
#include "trailapoptosisparams.h"

HybridOdeSim::HybridOdeSim(TrailApoptosisParams *TrailApop_Params, GeneToProtParams *GtP_Params)
{
    // numerical parameters
    Doub AbsTol = 1e-3; // for the ode solver
    Doub RelTol = 1e-3; // for the ode solver


    // init future state vectors
    FutureGeneTable = MatDoub(GtP_Params->NumGenes,MaxEvents,0.);
    FutureMrnaTable = MatDoub(GtP_Params->NumGenes,MaxEvents,0.);
    FutureGeneMrnaEventsTable = VecDoub(MaxEvents,0.);

    stochRhs = new StochRhs (TrailApop_Params,GtP_Params) ;

    // prepare simulator
    CellOutput = new Output();

     OdeIntStoch = new Odeint<StepperSie<StochRhs> > (58,0.,0.,AbsTol,RelTol,0.1,0.,*CellOutput,*stochRhs);
    //    OdeIntStoch = new Odeint<StepperDopr5<StochRhs> > (58,0.,0.,AbsTol,RelTol,0.1,0.,*CellOutput,*stochRhs);
//   OdeIntStoch = new Odeint<StepperDopr853<StochRhs> > (58,0.,0.,AbsTol,RelTol,0.1,0.,*CellOutput,*stochRhs);

}


