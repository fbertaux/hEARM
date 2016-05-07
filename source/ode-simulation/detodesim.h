#ifndef DETODESIM_H
#define DETODESIM_H

#include "detrhs.h"
#include "hearmodeint.h"

using namespace hearmOdeInt ;

class DetOdeSim
{
public:

    DetOdeSim(TrailApoptosisParams *TrailApop_Params, GeneToProtParams *GtP_Params);

    // containing params and equations
    DetRhs* detRhs;

    hearmOdeInt::Output* CellOutput;
//    Odeint<StepperDopr853<DetRhs> >* OdeIntDet;
    hearmOdeInt::Odeint<StepperSie<DetRhs> >* OdeIntDet;

};

#endif // DETODESIM_H
