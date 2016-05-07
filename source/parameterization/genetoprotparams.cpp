#include "genetoprotparams.h"

GeneToProtParams::GeneToProtParams ()
{
    b = VecDoub(17);
    c = VecDoub(17);
    H = VecDoub(17);
    rm = VecDoub(17);
    K = VecDoub(17);
    rp = VecDoub(17);
    Peq = VecDoub(17);

    int i=0 ;

    b[i]=0.0027923;	c[i]=0.00010774	;	rm[i]=2.1393e-05;	H[i]=0.0097889	;	rp[i]=7.1311e-06	;	K[i]=0.00041948    ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=1.734e-05;	c[i]=1.156e-05	;	rm[i]=0.00019254;	H[i]=0.008183	;	rp[i]=0.00048135	;	K[i]=0.05663       ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0027923;	c[i]=0.00010774	;	rm[i]=2.1393e-05;	H[i]=0.0097889	;	rp[i]=7.1311e-06	;	K[i]=0.0041948     ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0027923;	c[i]=0.00010774	;	rm[i]=2.1393e-05;	H[i]=0.0097889	;	rp[i]=7.1311e-06	;	K[i]=0.00041948    ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0028198;	c[i]=8.0229e-05	;	rm[i]=2.1393e-05;	H[i]=0.013146	;	rp[i]=7.1311e-06	;	K[i]=0.0041948     ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0027923;	c[i]=0.00010774	;	rm[i]=2.1393e-05;	H[i]=0.0097889	;	rp[i]=7.1311e-06	;	K[i]=0.0041948     ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0028237;	c[i]=7.6314e-05	;	rm[i]=2.1393e-05;	H[i]=0.01382	;	rp[i]=7.1311e-06	;	K[i]=0.041948      ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0027923;	c[i]=0.00010774	;	rm[i]=2.1393e-05;	H[i]=0.0097889	;	rp[i]=7.1311e-06	;	K[i]=0.41948       ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0028237;	c[i]=7.6314e-05	;	rm[i]=2.1393e-05;	H[i]=0.01382	;	rp[i]=7.1311e-06	;	K[i]=0.025169      ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=1.734e-05;	c[i]=1.156e-05	;	rm[i]=0.00019254;	H[i]=0.008183	;	rp[i]=0.00048135	;	K[i]=0.5663        ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0028117;	c[i]=8.8285e-05	;	rm[i]=2.1393e-05;	H[i]=0.011946	;	rp[i]=7.1311e-06	;	K[i]=0.033558      ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0028273;	c[i]=7.2692e-05	;	rm[i]=2.1393e-05;	H[i]=0.014509	;	rp[i]=7.1311e-06	;	K[i]=0.012584      ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0027923;	c[i]=0.00010774	;	rm[i]=2.1393e-05;	H[i]=0.0097889	;	rp[i]=7.1311e-06	;	K[i]=0.20974       ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0027923;	c[i]=0.00010774	;	rm[i]=2.1393e-05;	H[i]=0.0097889	;	rp[i]=7.1311e-06	;	K[i]=0.20974       ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0027923;	c[i]=0.00010774	;	rm[i]=2.1393e-05;	H[i]=0.0097889	;	rp[i]=7.1311e-06	;	K[i]=0.041948      ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0027923;	c[i]=0.00010774	;	rm[i]=2.1393e-05;	H[i]=0.0097889	;	rp[i]=7.1311e-06	;	K[i]=0.041948      ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;
    b[i]=0.0027923;	c[i]=0.00010774	;	rm[i]=2.1393e-05;	H[i]=0.0097889	;	rp[i]=7.1311e-06	;	K[i]=0.041948      ;   Peq[i]=giveMrnaMean(i)*K[i]/rp[i];    i++	;

}




void
GeneToProtParams::changeMcl1Model (string mcl1_model_type)
{
    int mcl1_index = 9 ;
    qDebug () << "Changing Mcl1, before : \n" << "b=" << b[mcl1_index] << "c=" << c[mcl1_index] << "rm=" << rm[mcl1_index] << "H=" << H[mcl1_index] << "rp=" << rp[mcl1_index] << "K=" << K[mcl1_index]  ;
    if ( mcl1_model_type == "NonFitted" )
    {
        //        0.0027912	0.00010875	9.627e-05	0.043643	0.00038508	0.45304
        b[mcl1_index]= 0.0027912 ;
        c[mcl1_index]= 0.00010875 ;
        H[mcl1_index]= 0.043643 ;
        rm[mcl1_index]= 9.627e-05 ;
        K[mcl1_index]= 0.45304 ;
        rp[mcl1_index]= 0.00038508 ;
    }
    else if (mcl1_model_type == "Fitted" )
    {
        //        1.734e-05	1.156e-05	0.00019254	0.008183	0.00048135	0.5663
        b[mcl1_index]=1.734e-05 ;
        c[mcl1_index]=1.156e-05 ;
        H[mcl1_index]=0.008183 ;
        rm[mcl1_index]=0.00019254 ;
        K[mcl1_index]=0.5663 ;
        rp[mcl1_index]=0.00048135 ;
    }
    else
    {
        qDebug () << "Unknown Mcl1 model type asked. I quit." ;
        exit (10) ;
    }
    qDebug () << "after : \n b=" << b[mcl1_index] << "c=" << c[mcl1_index] << "rm=" << rm[mcl1_index] << "H=" << H[mcl1_index] << "rp=" << rp[mcl1_index] << "K=" << K[mcl1_index]  ;
    // recompute Peq (but should not change normally)
    Peq[mcl1_index]=giveMrnaMean(mcl1_index)*K[mcl1_index]/rp[mcl1_index] ;
}

void
GeneToProtParams::changeFlipModel (string flip_model_type)
{
    int flip_index = 1 ;
    qDebug () << "Changing Flip, before : \n" << "b=" << b[flip_index] << "c=" << c[flip_index] << "rm=" << rm[flip_index] << "H=" << H[flip_index] << "rp=" << rp[flip_index] << "K=" << K[flip_index]  ;
    if ( flip_model_type == "NonFitted" )
    {
        //        0.0027912	0.00010875	9.627e-05	0.043643	0.00038508	0.045304
        b[flip_index]= 0.0027912 ;
        c[flip_index]= 0.00010875 ;
        H[flip_index]= 0.043643 ;
        rm[flip_index]= 9.627e-05 ;
        K[flip_index]= 0.045304 ;
        rp[flip_index]= 0.00038508 ;
    }
    else if (flip_model_type == "Fitted" )
    {
        //        1.734e-05	1.156e-05	0.00019254	0.008183	0.00048135	0.05663
        b[flip_index]=1.734e-05 ;
        c[flip_index]=1.156e-05 ;
        H[flip_index]=0.008183 ;
        rm[flip_index]=0.00019254 ;
        K[flip_index]=0.05663 ;
        rp[flip_index]=0.00048135 ;
    }
    else
    {
        qDebug () << "Unknown flip model type asked. I quit." ;
        exit (10) ;
    }
    qDebug () << "after : \n b=" << b[flip_index] << "c=" << c[flip_index] << "rm=" << rm[flip_index] << "H=" << H[flip_index] << "rp=" << rp[flip_index] << "K=" << K[flip_index]  ;
    // recompute Peq (but should not change normally)
    Peq[flip_index]=giveMrnaMean(flip_index)*K[flip_index]/rp[flip_index] ;
}

// Analytics Methods

Doub
GeneToProtParams::giveGeneR(Int i)
{
    return b[i]+c[i];
}

Doub
GeneToProtParams::giveGeneMean(Int i)
{
    return c[i]/giveGeneR(i);
}

Doub
GeneToProtParams::giveMrnaMean(Int i)
{
    return giveGeneMean(i)*H[i]/rm[i];
}

Doub
GeneToProtParams::giveProtMean(Int i)
{
    return giveMrnaMean(i)*K[i]/rp[i];
}

Doub
GeneToProtParams::giveProtVar(Int i)
{
    Doub fac1 = rp[i]*rm[i]*(rp[i]+rm[i]+giveGeneR(i))/(rm[i]+giveGeneR(i))/(rp[i]+giveGeneR(i))/(rp[i]+rm[i]);
    fac1 *= (1-giveGeneMean(i))/giveGeneMean(i)*giveProtMean(i)*giveProtMean(i);

    Doub fac2 = rp[i]/(rp[i]+rm[i])*giveProtMean(i)*giveProtMean(i)/giveMrnaMean(i);

    return fac1+fac2;
}

Doub
GeneToProtParams::giveProtCV(Int i)
{
    return sqrt(giveProtVar(i))/giveProtMean(i);
}

Doub
GeneToProtParams::giveMrnaProtCov(Int i)
{
    return giveProtVar(i)*giveMrnaMean(i)/giveProtMean(i);
}

Doub
GeneToProtParams::giveGeneProtCov(Int i)
{
    return giveProtMean(i)*(1-giveGeneMean(i))*rp[i]*rm[i]/(rp[i]+giveGeneR(i))/(rm[i]+giveGeneR(i));
}

Doub
GeneToProtParams::giveCp(Int i,Doub t)
{
    return exp(-rp[i]*t);
}

Doub
GeneToProtParams::giveCm(Int i,Doub t)
{
    return giveProtMean(i)/giveMrnaMean(i)*rp[i]/(rp[i]-rm[i])*(exp(-rm[i]*t)-exp(-rp[i]*t));
}

Doub
GeneToProtParams::giveCg(Int i,Doub t)
{
    return giveProtMean(i)/giveGeneMean(i)*rm[i]*rp[i]/(rp[i]-rm[i])/(rm[i]-giveGeneR(i))/(rp[i]-giveGeneR(i))*((rm[i]-giveGeneR(i))*exp(-rp[i]*t) + (rp[i]-rm[i])*exp(-giveGeneR(i)*t) + (giveGeneR(i)-rp[i])*exp(-rm[i]*t));
}

Doub
GeneToProtParams::giveProtAutocorr(Int i,Doub t)
{
    return (giveProtVar(i)*giveCp(i,t)+giveMrnaProtCov(i)*giveCm(i,t)+giveGeneProtCov(i)*giveCg(i,t))/giveProtVar(i);
}

// Display methods

void
GeneToProtParams::showAllProtMeans()
{
    for (int i=0;i<NumGenes;i++)
    {
        cout << "Mean prot " << i << " : " << giveProtMean(i) << endl;
    }
}

void
GeneToProtParams::showAllProtCVs()
{
    for (int i=0;i<NumGenes;i++)
    {
        cout << "CV prot " << i << " : " << giveProtCV(i) << endl;
    }
}








GeneToProtParams::GeneToProtParams(string path)
{
    ostringstream FilesPathOs;
    FilesPathOs << path << "parameters/GeneToProt_Params_ForC.txt";
    string FilesPath = FilesPathOs.str();
    cout << FilesPath << endl;
    ifstream myFlux(FilesPath.c_str() );

    if (!myFlux.is_open()) { cout << "Can't parse Gene to Prot parameters ! I quit." << endl; exit(1); }

    b = VecDoub(17);
    c = VecDoub(17);
    H = VecDoub(17);
    rm = VecDoub(17);
    K = VecDoub(17);
    rp = VecDoub(17);
    Peq = VecDoub(17);

    int indforc;
    Doub bb,cc,rmm,Hh,rpp,Kk;
    for (int i=0;i<NumGenes;i++)
    {
        myFlux >> indforc;
        myFlux >> bb;
        myFlux >> cc;
        myFlux >> rmm;
        myFlux >> Hh;
        myFlux >> rpp;
        myFlux >> Kk;

        b[i]=bb; c[i]=cc; rm[i]=rmm; H[i]=Hh;
        rp[i]=rpp; K[i]=Kk;

        Peq[i]=giveMrnaMean(i)*K[i]/rp[i];
    }
}
