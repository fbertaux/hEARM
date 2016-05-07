#ifndef HYBRIDODESIM_H
#define HYBRIDODESIM_H

#include "genemrnasim.h"
#include "stochrhs.h"
#include "hearmodeint.h"

using namespace hearmOdeInt;

class HybridOdeSim
{
public:

    HybridOdeSim(TrailApoptosisParams *TrailApop_Params, GeneToProtParams *GtP_Params);


    // fields
    static const Int MaxEvents = 1e6;
    MatDoub FutureMrnaTable;
    MatDoub FutureGeneTable;
    VecDoub FutureGeneMrnaEventsTable;
    Int EventObtained;

    // containing params and equations
    StochRhs* stochRhs; // needed to have hybrid ODE-gillespie

    hearmOdeInt::Output* CellOutput;
    hearmOdeInt::Odeint<StepperSie<StochRhs> >* OdeIntStoch;
//    Odeint<StepperDopr5<StochRhs> >* OdeIntStoch;
//    hearmOdeInt::Odeint<StepperDopr853<StochRhs> >* OdeIntStoch;


};

#endif // HYBRIDODESIM_H
