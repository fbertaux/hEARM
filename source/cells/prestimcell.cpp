#include "prestimcell.h"


//! DEFAULT CONSTRUCTOR
PreStimCell::PreStimCell(GeneToProtParams* gtp_pars,GeneMrnaSim* gm_sim)
    : GeneToProtParameters (gtp_pars) , GeneMrna_Sim (gm_sim) , CellAge(0.)
{
    // init states vector
    GeneMrnaState = VecDoub(3*GeneToProtParameters->NumGenes,0.);
    ProtState = VecDoub(GeneToProtParameters->NumGenes,0.);
    for (Int i=0;i<GeneToProtParameters->NumGenes;i++)
    {
        GeneMrnaState[3*i] = 1; GeneMrnaState[3*i+1] = 0; // gene is on for all genes
        GeneMrnaState[3*i+2] = (Int) GeneToProtParameters->giveMrnaMean(i); // Mrna to mean level
        ProtState[i] = GeneToProtParameters->giveProtMean(i); // Prot to mean level
    }
}


//! COPY CONSTRUCTOR
PreStimCell::PreStimCell(PreStimCell &toCopy) // copy constructor
    : GeneToProtParameters (toCopy.GeneToProtParameters) , GeneMrna_Sim (toCopy.GeneMrna_Sim) , CellAge(0.)
{
    // init states vector
    GeneMrnaState = VecDoub(3*GeneToProtParameters->NumGenes,0.);
    ProtState = VecDoub(GeneToProtParameters->NumGenes,0.);
    for (Int i=0;i<GeneToProtParameters->NumGenes;i++)
    {
        GeneMrnaState[3*i] = toCopy.GeneMrnaState[3*i];
        GeneMrnaState[3*i+1] = toCopy.GeneMrnaState[3*i+1];
        GeneMrnaState[3*i+2] = toCopy.GeneMrnaState[3*i+2];
        ProtState[i] = toCopy.ProtState[i];
    }
}




// simulation methods

void
PreStimCell::simulate(Doub duration)
{
    // prepare mrna simulator for this cell
    GeneMrna_Sim->prepareForSteps(GeneMrnaState);
    Doub dt,nt;
    Doub t = 0;
    while (t<duration)
    {
        nt = GeneMrna_Sim->step(GeneMrnaState,duration);
        dt = nt - t;
        t = GeneMrna_Sim->t;
        for (Int i=0;i<GeneToProtParameters->NumGenes;i++)
        {
            ProtState[i] = GeneToProtParameters->K[i]*GeneMrnaState[3*i+2]/GeneToProtParameters->rp[i]+(ProtState[i]-GeneToProtParameters->K[i]*GeneMrnaState[3*i+2]/GeneToProtParameters->rp[i])*exp(-GeneToProtParameters->rp[i]*dt);
        }
    }
    CellAge+=duration;
}

