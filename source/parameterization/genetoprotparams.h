#ifndef GENETOPROTPARAMS_H
#define GENETOPROTPARAMS_H

#include "nr3.h"

struct GeneToProtParams
{

    // GeneToProt parameters

    VecDoub b,c,H,rm,K,rp;
    VecDoub Peq;

    static const Int NumGenes = 17;


    // constructor

    GeneToProtParams () ;
    GeneToProtParams(string path) ; //obsolete, kept if needed for legacy

    // method to change Mcl1 or Flip model
    void changeMcl1Model ( string mcl1_model_type ) ;
    void changeFlipModel ( string flip_model_type ) ;

    // analytical methods

    Doub giveProtMean(Int i);
    Doub giveProtVar(Int i);
    Doub giveProtCV(Int i);
    Doub giveProtMixingTime(Int i);
    Doub giveGeneMean(Int i);
    Doub giveMrnaMean(Int i);
    Doub giveMrnaVar(Int i);
    Doub giveGeneR(Int i);
    Doub giveMrnaProtCov(Int i);
    Doub giveGeneProtCov(Int i);
    Doub giveCp(Int i, Doub t);
    Doub giveCm(Int i, Doub t);
    Doub giveCg(Int i, Doub t);
    Doub giveProtAutocorr(Int i,Doub t);


    // display methods

    void showAllProtMeans();
    void showAllProtCVs();
    void showAllAutocs(Doub t);
    void showAllMrnaMeans();

};

#endif // GENETOPROTPARAMS_H
