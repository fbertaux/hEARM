#include "trailapoptosisparams.h"


TrailApoptosisParams::TrailApoptosisParams ()
{
    num_species = 58 ;
    num_reactions = 71 ;

    // native forms names to index mapping
    nativeNamesToIndex [ "R" ] = 1 ;
    nativeNamesToIndex [ "flip" ] = 3 ;
    nativeNamesToIndex [ "pC8" ] = 4 ;
    nativeNamesToIndex [ "Bar" ] = 6 ;
    nativeNamesToIndex [ "pC3" ] = 7 ;
    nativeNamesToIndex [ "pC6" ] = 9 ;
    nativeNamesToIndex [ "XIAP" ] = 11 ;
    nativeNamesToIndex [ "PARP" ] = 12 ;
    nativeNamesToIndex [ "Bid" ] = 14 ;
    nativeNamesToIndex [ "Mcl1" ] = 16 ;
    nativeNamesToIndex [ "Bax" ] = 17 ;
    nativeNamesToIndex [ "Bcl2" ] = 22 ;
    nativeNamesToIndex [ "Pore" ] = 23 ;
    nativeNamesToIndex [ "CyCm" ] = 25 ;
    nativeNamesToIndex [ "Smacm" ] = 28 ;
    nativeNamesToIndex [ "Apaf" ] = 31 ;
    nativeNamesToIndex [ "pC9" ] = 33 ;

    // native forms names to index mapping
    nativeNamesToNativeIndex [ "R" ] = 0 ;
    nativeNamesToNativeIndex [ "flip" ] = 1 ;
    nativeNamesToNativeIndex [ "pC8" ] = 2 ;
    nativeNamesToNativeIndex [ "Bar" ] = 3 ;
    nativeNamesToNativeIndex [ "pC3" ] = 4 ;
    nativeNamesToNativeIndex [ "pC6" ] = 5 ;
    nativeNamesToNativeIndex [ "XIAP" ] = 6 ;
    nativeNamesToNativeIndex [ "PARP" ] = 7 ;
    nativeNamesToNativeIndex [ "Bid" ] = 8 ;
    nativeNamesToNativeIndex [ "Mcl1" ] = 9 ;
    nativeNamesToNativeIndex [ "Bax" ] = 10 ;
    nativeNamesToNativeIndex [ "Bcl2" ] = 11 ;
    nativeNamesToNativeIndex [ "Pore" ] = 12 ;
    nativeNamesToNativeIndex [ "CyCm" ] = 13 ;
    nativeNamesToNativeIndex [ "Smacm" ] = 14 ;
    nativeNamesToNativeIndex [ "Apaf" ] = 15 ;
    nativeNamesToNativeIndex [ "pC9" ] = 16 ;

    // active forms names to index mapping
    activeNamesToIndex["R*"]=2;
    activeNamesToIndex["C8"]=5;
    activeNamesToIndex["C3"] = 8;
    activeNamesToIndex["C6"] = 10;
    activeNamesToIndex["cPARP"] = 13;
    activeNamesToIndex["tBid"] = 15;
    activeNamesToIndex["Bax*"] = 18;
    activeNamesToIndex["Bax*m"] = 19;
    activeNamesToIndex["Bax2"] = 20;
    activeNamesToIndex["Bax4"] = 21;
    activeNamesToIndex["Pore*"] = 24;
    activeNamesToIndex["CyCr"] = 26;
    activeNamesToIndex["CyC"] = 27;
    activeNamesToIndex["Smacr"] = 29;
    activeNamesToIndex["Smac"] = 30;
    activeNamesToIndex["Apaf*"] = 32;
    activeNamesToIndex["Apop"] = 34;
    activeNamesToIndex["C3_U"] = 35;
    activeNamesToIndex["L:R"] = 36;
    activeNamesToIndex["R*:flip"] = 37;
    activeNamesToIndex["R*:pC8"] = 38;
    activeNamesToIndex["C8:Bar"] = 39;
    activeNamesToIndex["C8:pC3"] = 40;
    activeNamesToIndex["C3:pC6"] = 41;
    activeNamesToIndex["C6:pC8"] = 42;
    activeNamesToIndex["C3:XIAP"] = 43;
    activeNamesToIndex["C3:PARP"] = 44;
    activeNamesToIndex["C8:Bid"] = 45;
    activeNamesToIndex["tBid:Mcl1"] = 46;
    activeNamesToIndex["tBid:Bax"] = 47;
    activeNamesToIndex["Bax*m:Bcl2"] = 48;
    activeNamesToIndex["Bax2:Bcl2"] = 49;
    activeNamesToIndex["Bax4:Bcl2"] = 50;
    activeNamesToIndex["Bax4:M"] = 51;
    activeNamesToIndex["M*:CyCm"] = 52;
    activeNamesToIndex["M*:Smacm"] = 53;
    activeNamesToIndex["CyC:Apaf"] = 54;
    activeNamesToIndex["Apop:pC3"] = 55;
    activeNamesToIndex["Apop:XIAP"] = 56;
    activeNamesToIndex["Smac:XIAP"] = 57;

//    cout << "species names to index maps filled." << endl ;

    // compute the mapping native index -> index
    nativeIndexToIndex = new int[17] ;
    for (int i=0 ; i<17 ; i++ )
    {
        nativeIndexToIndex [i] = giveIndexFromNativeIndex (i) ;
    }
//    cout << "the mapping between native species index and their general index is now stored." << endl ;

    // compute the mapping index -> native index
    indexToNativeIndex = new int[58] ;
    for (int i=0 ; i<58 ; i++ )
    {
        if ( isNativeFromIndex (i) )
        {
            indexToNativeIndex [i] = giveNativeIndexFromIndex (i) ;
        }
        else
        {
            indexToNativeIndex [i] = -1 ;
        }
    }
//    cout << "The inverse mapping has been also computed and stored." << endl ;

    // initialize kinetic rates ( with feedback loop )
    KineticRates = VecDoub (71,0.) ;
    setDefaultKineticRates () ;

    // remove feedback loop [ as in PCB paper ]
    removeFeedbackLoop () ;

    // set active half-lives to 5 [ except TRAIL and Pore*, as in PCB paper ]
    DegradationRates = VecDoub (58,0.) ;
    setAllActiveFormsHalfLife (5.) ;

    // set TRAIL and Pore* half-lives [ as in PCB paper ]
    DegradationRates [ activeNamesToIndex ["L"] ] = log (2.) / 9. / 3600. ;
    DegradationRates [ activeNamesToIndex ["Pore*"] ] = log (2.) / 1.9 / 3600. ;

    // set Mcl1 forms half-lives to 0.4 [ as in PCB paper ]
    setMcl1FormsHalfLife (0.4) ;

    // set Flip forms half-lives to 0.4 [ as in PCB paper ]
    setFlipFormsHalfLife (0.4) ;

}



void
TrailApoptosisParams::setDefaultKineticRates ()
{
    KineticRates[0]=4e-07; //R1_1
    KineticRates[1]=1e-06; //R1_2
    KineticRates[2]=0.01; //R1_3
    KineticRates[3]=1e-06; //R2_1
    KineticRates[4]=0.001; //R2_2
    KineticRates[5]=1e-07; //R3_1
    KineticRates[6]=0.001; //R3_2
    KineticRates[7]=1; //R3_3
    KineticRates[8]=1e-06; //R4_1
    KineticRates[9]=0.001; //R4_2
    KineticRates[10]=1e-07; //R5_1
    KineticRates[11]=0.001; //R5_2
    KineticRates[12]=1; //R5_3
    KineticRates[13]=1e-07; //R6_1
    KineticRates[14]=0.001; //R6_2
    KineticRates[15]=1; //R6_3
    KineticRates[16]=1e-07; //R7_1
    KineticRates[17]=0.001; //R7_2
    KineticRates[18]=1; //R7_3
    KineticRates[19]=2e-06; //R8_1
    KineticRates[20]=0.001; //R8_2
    KineticRates[21]=0.1; //R8_3
    KineticRates[22]=1e-06; //R9_1
    KineticRates[23]=0.001; //R9_2
    KineticRates[24]=20; //R9_3
    KineticRates[25]=1e-07; //R10_1
    KineticRates[26]=0.001; //R10_2
    KineticRates[27]=1; //R10_3
    KineticRates[28]=1e-06; //R11_1
    KineticRates[29]=0.001; //R11_2
    KineticRates[30]=1e-07; //R12_1
    KineticRates[31]=0.001; //R12_2
    KineticRates[32]=1; //R12_3
    KineticRates[33]=0.01; //R13_1
    KineticRates[34]=1; //R13_2
    KineticRates[35]=1.429e-05; //R14_1
    KineticRates[36]=0.001; //R14_2
    KineticRates[37]=1.429e-05; //R15_1
    KineticRates[38]=0.001; //R15_2
    KineticRates[39]=1.429e-05; //R16_1
    KineticRates[40]=0.001; //R16_2
    KineticRates[41]=1.429e-05; //R17_1
    KineticRates[42]=0.001; //R17_2
    KineticRates[43]=1.429e-05; //R18_1
    KineticRates[44]=0.001; //R18_2
    KineticRates[45]=1.429e-05; //R19_1
    KineticRates[46]=0.001; //R19_2
    KineticRates[47]=1; //R19_3
    KineticRates[48]=2.857e-05; //R20_1
    KineticRates[49]=0.001; //R20_2
    KineticRates[50]=10; //R20_3
    KineticRates[51]=2.857e-05; //R21_1
    KineticRates[52]=0.001; //R21_2
    KineticRates[53]=10; //R21_3
    KineticRates[54]=1; //R22_1
    KineticRates[55]=0.01; //R22_2
    KineticRates[56]=5e-07; //R23_1
    KineticRates[57]=0.001; //R23_2
    KineticRates[58]=1; //R23_3
    KineticRates[59]=5e-08; //R24_1
    KineticRates[60]=0.001; //R24_2
    KineticRates[61]=5e-09; //R25_1
    KineticRates[62]=0.001; //R25_2
    KineticRates[63]=1; //R25_3
    KineticRates[64]=1; //R26_1
    KineticRates[65]=0.01; //R26_2
    KineticRates[66]=2e-06; //R27_1
    KineticRates[67]=0.001; //R27_2
    KineticRates[68]=7e-06; //R28_1
    KineticRates[69]=0.001; //R28_2
    KineticRates[70]=0.01; //R29_1
}

void
TrailApoptosisParams::removeFeedbackLoop ()
{
    KineticRates[13]=0.;
    KineticRates[14]=0.;
    KineticRates[15]=0.;
}

void
TrailApoptosisParams::setAllActiveFormsHalfLife (double hl_in_hours)
{
    // iterate on active forms from the map
    typedef map<string, int>::iterator iter ;
    for ( iter iterator = activeNamesToIndex.begin() ; iterator != activeNamesToIndex.end() ; iterator++ )
    {
        if ( iterator->first != "L" && iterator->first != "Pore*" )
        {
            DegradationRates [ iterator->second ] = log(2.) / hl_in_hours / 3600. ;
        }
    }
}

void
TrailApoptosisParams::setMcl1FormsHalfLife (double hl_in_hours)
{
    DegradationRates [ nativeNamesToIndex ["Mcl1"] ] = log (2.) / hl_in_hours / 3600. ;
    DegradationRates [ activeNamesToIndex ["tBid:Mcl1"] ] = log (2.) / hl_in_hours / 3600. ;
}

void
TrailApoptosisParams::setFlipFormsHalfLife (double hl_in_hours)
{
    DegradationRates [ nativeNamesToIndex ["flip"] ] = log (2.) / hl_in_hours / 3600. ;
    DegradationRates [ activeNamesToIndex ["R*:flip"] ] = log (2.) / hl_in_hours / 3600. ;
}


bool
TrailApoptosisParams::isNativeFromIndex (int index)
{
    bool found = false ;
    // iterate on native forms from the map
    map<string, int>::iterator iter ;
    for ( iter = nativeNamesToIndex.begin() ; iter != nativeNamesToIndex.end() ; iter++ )
    {
        if ( iter->second == index ) found = true ;
    }
    return found ;
}

int
TrailApoptosisParams::giveIndexFromNativeIndex (int natIndex)
{
    map<string, int>::iterator iter = nativeNamesToNativeIndex.begin() ;
    while ( iter->second != natIndex ) iter++ ;
//    cout << "found that species " << iter->first.c_str() << " has the native index " << natIndex << endl ;
    return nativeNamesToIndex [ iter->first ] ;
}

int
TrailApoptosisParams::giveNativeIndexFromIndex (int index)
{
    map<string, int>::iterator iter = nativeNamesToIndex.begin() ;
    while ( iter->second != index ) iter++ ;
//    cout << "found that species " << iter->first.c_str() << " has the general index " << index << endl ;
    return nativeNamesToNativeIndex [ iter->first ] ;
}

string
TrailApoptosisParams::giveNameFromIndex (int index)
{
    map<string, int>::iterator iter = nativeNamesToIndex.begin() ;
    while ( iter->second != index && iter != nativeNamesToIndex.end() ) iter++ ;
    if ( iter != nativeNamesToIndex.end() )
    {
        return iter->first ;
    }
    else
    {
        map<string, int>::iterator iter = activeNamesToIndex.begin() ;
        while ( iter->second != index && iter != activeNamesToIndex.end() ) iter++ ;
        if ( iter != activeNamesToIndex.end() )
        {
            return iter->first ;
        }
        else { cout << "index not found? I quit." << endl ; exit(3) ; }
    }
}


