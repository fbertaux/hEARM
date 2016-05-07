#include "detodesim.h"

DetOdeSim::DetOdeSim(TrailApoptosisParams *TrailApop_Params, GeneToProtParams *GtP_Params)
{
    // numerical parameters
    Doub AbsTol = 1e-3; // for the ode solver
    Doub RelTol = 1e-3; // for the ode solver

    detRhs = new DetRhs(TrailApop_Params,GtP_Params);

    // prepare simulator
    CellOutput = new Output();
//    OdeIntDet = new Odeint<StepperDopr853<DetRhs> > (58,0.,0.,AbsTol,RelTol,0.1,0.,*CellOutput,*detRhs);
    OdeIntDet = new Odeint<StepperSie<DetRhs> > (58,0.,0.,AbsTol,RelTol,0.1,0.,*CellOutput,*detRhs);
}
