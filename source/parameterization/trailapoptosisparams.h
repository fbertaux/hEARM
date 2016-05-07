#ifndef TRAILAPOPTOSISPARAMS_H
#define TRAILAPOPTOSISPARAMS_H

#include "nr3.h"



struct TrailApoptosisParams
{

    // description of species and reactions

    map <string,int> nativeNamesToIndex ;
    map <string,int> nativeNamesToNativeIndex ;
    map <string,int> activeNamesToIndex ;

    int num_species ;
    int num_reactions ;

    int* indexToNativeIndex ;
    int* nativeIndexToIndex ;


    // parameterization

    VecDoub KineticRates ; // per s
    VecDoub DegradationRates ; // per s


    // constructor
    TrailApoptosisParams () ; // default constructor

    // methods to changes paramters
    void setDefaultKineticRates () ;
    void removeFeedbackLoop () ;
    void setTrailHalfLife () ;
    void setAllActiveFormsHalfLife ( double hl_in_hours ) ; // except L and Pore*
    void setMcl1FormsHalfLife ( double hl_in_hours ) ;
    void setFlipFormsHalfLife ( double hl_in_hours ) ;

    // methods to know get infos about the parameters
    bool isNativeFromIndex ( int index ) ;
    int giveIndexFromNativeIndex ( int natIndex ) ;
    int giveNativeIndexFromIndex ( int index ) ;
    string giveNameFromIndex ( int index ) ;

};

#endif // TRAILAPOPTOSISPARAMS_H
