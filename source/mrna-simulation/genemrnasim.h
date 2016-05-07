#ifndef GENEMRNASIM_H
#define GENEMRNASIM_H


#include "nr3.h"
#include "sparse.h"
#include "ran.h"

void sparmatfill(NRvector<NRsparseCol> &sparmat, MatDoub &fullmat);

#include "genetoprotparams.h"

struct GeneMrnaSim
{

    Int mm;     // number of reactions
    Int nn;     // number of species
    VecDoub &b,&c,&H,&rm; // parameters

    VecDoub a;  // reaction rates
    MatDoub instate, outstate; // reactants matrix, change state matrix
    NRvector<NRsparseCol> outchg, depend; // change state sparse matrix, reaction dependancy sparse matrix
    VecInt pr;  // priority list for reaction
    Doub t; // time
    Doub asum; // sum of all reactions rates
    Ran ran; // random generator

    Bool readyForSteps;

    // CONSTRUCTOR
    GeneMrnaSim(GeneToProtParams *GtP_P,Int seed=1);
    GeneMrnaSim(Int mmm,Int nnn,VecDoub &bb,VecDoub &cc,VecDoub &HH,VecDoub &rmrm,Int seed=1);

    // PREPARATION
    void prepareForSteps(VecDoub &s); // compute asum and a[i]

    // SIMULATIONS
    Doub step(VecDoub &s , Doub targetTime);

    typedef Doub(GeneMrnaSim::*rateptr)(VecDoub &s);
    rateptr *dispatch;


    void describereactions();


    ~GeneMrnaSim() {delete [] dispatch;}

    Doub rate0(VecDoub &s) {return b[0]*s[0];}
    Doub rate1(VecDoub &s) {return c[0]*s[1];}
    Doub rate2(VecDoub &s) {return H[0]*s[0];}
    Doub rate3(VecDoub &s) {return rm[0]*s[2];}
    Doub rate4(VecDoub &s) {return b[1]*s[3];}
    Doub rate5(VecDoub &s) {return c[1]*s[4];}
    Doub rate6(VecDoub &s) {return H[1]*s[3];}
    Doub rate7(VecDoub &s) {return rm[1]*s[5];}
    Doub rate8(VecDoub &s) {return b[2]*s[6];}
    Doub rate9(VecDoub &s) {return c[2]*s[7];}
    Doub rate10(VecDoub &s) {return H[2]*s[6];}
    Doub rate11(VecDoub &s) {return rm[2]*s[8];}
    Doub rate12(VecDoub &s) {return b[3]*s[9];}
    Doub rate13(VecDoub &s) {return c[3]*s[10];}
    Doub rate14(VecDoub &s) {return H[3]*s[9];}
    Doub rate15(VecDoub &s) {return rm[3]*s[11];}
    Doub rate16(VecDoub &s) {return b[4]*s[12];}
    Doub rate17(VecDoub &s) {return c[4]*s[13];}
    Doub rate18(VecDoub &s) {return H[4]*s[12];}
    Doub rate19(VecDoub &s) {return rm[4]*s[14];}
    Doub rate20(VecDoub &s) {return b[5]*s[15];}
    Doub rate21(VecDoub &s) {return c[5]*s[16];}
    Doub rate22(VecDoub &s) {return H[5]*s[15];}
    Doub rate23(VecDoub &s) {return rm[5]*s[17];}
    Doub rate24(VecDoub &s) {return b[6]*s[18];}
    Doub rate25(VecDoub &s) {return c[6]*s[19];}
    Doub rate26(VecDoub &s) {return H[6]*s[18];}
    Doub rate27(VecDoub &s) {return rm[6]*s[20];}
    Doub rate28(VecDoub &s) {return b[7]*s[21];}
    Doub rate29(VecDoub &s) {return c[7]*s[22];}
    Doub rate30(VecDoub &s) {return H[7]*s[21];}
    Doub rate31(VecDoub &s) {return rm[7]*s[23];}
    Doub rate32(VecDoub &s) {return b[8]*s[24];}
    Doub rate33(VecDoub &s) {return c[8]*s[25];}
    Doub rate34(VecDoub &s) {return H[8]*s[24];}
    Doub rate35(VecDoub &s) {return rm[8]*s[26];}
    Doub rate36(VecDoub &s) {return b[9]*s[27];}
    Doub rate37(VecDoub &s) {return c[9]*s[28];}
    Doub rate38(VecDoub &s) {return H[9]*s[27];}
    Doub rate39(VecDoub &s) {return rm[9]*s[29];}
    Doub rate40(VecDoub &s) {return b[10]*s[30];}
    Doub rate41(VecDoub &s) {return c[10]*s[31];}
    Doub rate42(VecDoub &s) {return H[10]*s[30];}
    Doub rate43(VecDoub &s) {return rm[10]*s[32];}
    Doub rate44(VecDoub &s) {return b[11]*s[33];}
    Doub rate45(VecDoub &s) {return c[11]*s[34];}
    Doub rate46(VecDoub &s) {return H[11]*s[33];}
    Doub rate47(VecDoub &s) {return rm[11]*s[35];}
    Doub rate48(VecDoub &s) {return b[12]*s[36];}
    Doub rate49(VecDoub &s) {return c[12]*s[37];}
    Doub rate50(VecDoub &s) {return H[12]*s[36];}
    Doub rate51(VecDoub &s) {return rm[12]*s[38];}
    Doub rate52(VecDoub &s) {return b[13]*s[39];}
    Doub rate53(VecDoub &s) {return c[13]*s[40];}
    Doub rate54(VecDoub &s) {return H[13]*s[39];}
    Doub rate55(VecDoub &s) {return rm[13]*s[41];}
    Doub rate56(VecDoub &s) {return b[14]*s[42];}
    Doub rate57(VecDoub &s) {return c[14]*s[43];}
    Doub rate58(VecDoub &s) {return H[14]*s[42];}
    Doub rate59(VecDoub &s) {return rm[14]*s[44];}
    Doub rate60(VecDoub &s) {return b[15]*s[45];}
    Doub rate61(VecDoub &s) {return c[15]*s[46];}
    Doub rate62(VecDoub &s) {return H[15]*s[45];}
    Doub rate63(VecDoub &s) {return rm[15]*s[47];}
    Doub rate64(VecDoub &s) {return b[16]*s[48];}
    Doub rate65(VecDoub &s) {return c[16]*s[49];}
    Doub rate66(VecDoub &s) {return H[16]*s[48];}
    Doub rate67(VecDoub &s) {return rm[16]*s[50];}

};

#endif // GENEMRNASIM_H
