#ifndef DETRHS_H
#define DETRHS_H

#include "trailapoptosisparams.h"
#include "genetoprotparams.h"

struct DetRhs
{

    // CONSTRUCTOR
    DetRhs(TrailApoptosisParams *TrailApop_Params, GeneToProtParams *GtP_Params);


    VecDoub k; // kinetic rates
    VecDoub deg; // degradation rates
    VecDoub synth; // synthesis rates

    TrailApoptosisParams* TrailApop_Parameters;
    GeneToProtParams* GtP_Parameters;

    // derivs & jacobian
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx);
    void jacobian(const Doub x, VecDoub_I &y, VecDoub_O &dfdx,MatDoub_O &dfdy);


    // debug

    void displayAllHalfLives();


};

#endif // DETRHS_H
