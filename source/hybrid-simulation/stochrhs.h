#ifndef STOCHRHS_H
#define STOCHRHS_H

#include "trailapoptosisparams.h"
#include "genetoprotparams.h"

struct StochRhs
{

    // CONSTRUCTOR
    StochRhs(TrailApoptosisParams *TrailApop_Params, GeneToProtParams *GtP_Params);


    VecDoub k; // kinetic rates

    VecDoub deg; // degradation rates
    VecDoub synth; // synthesis rates

    TrailApoptosisParams* TrailApop_Parameters;
    GeneToProtParams* GtP_Parameters;


    int tindex;
    MatDoub *mrnatable;
    VecDoub *timetable;
    int nevents;



    void setMrnaTables(MatDoub *mrtab, VecDoub* ttab, int neventss);
    void findGoodTindex(const Doub x);
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx);
    void jacobian(const Doub x, VecDoub_I &y, VecDoub_O &dfdx,MatDoub_O &dfdy);


    // DEBUG

    void displayAllHalfLives();
    void displaySynthesisRates();
    void displayRhsPeqs();

};

#endif // STOCHRHS_H
