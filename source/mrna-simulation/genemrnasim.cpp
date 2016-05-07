#include "genemrnasim.h"

void sparmatfill(NRvector<NRsparseCol> &sparmat, MatDoub &fullmat)
{
    Int n,m,nz,nn=fullmat.nrows(),mm=fullmat.ncols();
    if (sparmat.size() != mm) throw("bad sizes");
    for (m=0;m<mm;m++)
    {
        for (nz=n=0;n<nn;n++) if (fullmat[n][m]) nz++;
        sparmat[m].resize(nn,nz);
        for (nz=n=0;n<nn;n++) if (fullmat[n][m])
        {
            sparmat[m].row_ind[nz] = n;
            sparmat[m].val[nz++] = fullmat[n][m];
        }
    }
}


// CONSTRUCTORS

GeneMrnaSim::GeneMrnaSim(GeneToProtParams *GtP_P, Int seed)
    : mm(4*GtP_P->NumGenes), nn(3*GtP_P->NumGenes),
    b(GtP_P->b), c(GtP_P->c), H(GtP_P->H), rm(GtP_P->rm),
    a(mm,0.), outchg(mm), depend(mm), pr(mm), t(0.),asum(0.),
    ran(seed), readyForSteps(false),
    dispatch(new rateptr[mm])
{
    Int i,j,k,d;
    describereactions();
    sparmatfill(outchg,outstate);
    MatDoub dep(mm,mm);
    for (i=0;i<mm;i++) for (j=0;j<mm;j++)
    {
        d = 0;
        for (k=0;k<nn;k++) d = d || (instate[k][i] && outstate[k][j]);
        dep[i][j] = d;
    }
    sparmatfill(depend,dep);
    for (i=0;i<mm;i++)
    {
        pr[i] = i;
    }
}

GeneMrnaSim::GeneMrnaSim(Int mmm, Int nnn, VecDoub &bb, VecDoub &cc, VecDoub &HH, VecDoub &rmrm,Int seed)
    : mm(mmm), nn(nnn), b(bb), c(cc), H(HH), rm(rmrm),
      a(mm,0.), outchg(mm), depend(mm), pr(mm), t(0.),asum(0.),
      ran(seed), readyForSteps(false),
      dispatch(new rateptr[mm])
{
    Int i,j,k,d;
    describereactions();
    sparmatfill(outchg,outstate);
    MatDoub dep(mm,mm);
    for (i=0;i<mm;i++) for (j=0;j<mm;j++)
    {
        d = 0;
        for (k=0;k<nn;k++) d = d || (instate[k][i] && outstate[k][j]);
        dep[i][j] = d;
    }
    sparmatfill(depend,dep);
    for (i=0;i<mm;i++)
    {
        pr[i] = i;
    }
}

// PREPARATION

void
GeneMrnaSim::prepareForSteps(VecDoub &s)
{
    t=0;
    Int i;
    asum = 0.;
    for (i=0;i<mm;i++)
    {
        pr[i] = i;
        a[i] = (this->*dispatch[i])(s);
        asum += a[i];
    }
    readyForSteps = true;
}


// SIMULATION

Doub
GeneMrnaSim::step(VecDoub &s, Doub targetTime)
{
    Int i,n,m,k=0;
    Doub tau,atarg,sum,anew;
    if (asum == 0.) {t *= 2.; return t;}
    tau = -log(ran.doub())/asum;
    if (t+tau>targetTime)
    {
        t=targetTime;
        return t;
    }
    atarg = ran.doub()*asum;
    sum = a[pr[0]];
    while (sum < atarg) sum += a[pr[++k]];
    m = pr[k];
//    cout << "reaction : " << m << endl;
//    cout << "type : " << m % 4 << endl;
//    cout << "mrna : " << m/4 << endl;
    if (k > 0) SWAP(pr[k],pr[k-1]);
    if (k == mm-1) asum = sum;
    n = outchg[m].nvals;
    for (i=0;i<n;i++)
    {
        k = outchg[m].row_ind[i];
        s[k] += outchg[m].val[i];
//        cout << "adding " << outchg[m].val[i] << " to s[ " << k << " ]" << endl;
    }
    n = depend[m].nvals;
    for (i=0;i<n;i++)
    {
        k = depend[m].row_ind[i];
        anew = (this->*dispatch[k])(s);
        asum += (anew - a[k]);
        a[k] = anew;
    }
    if (t*asum < 0.1)
        for (asum=0.,i=0;i<mm;i++) asum += a[i];
    return (t += tau);
}


// OTHERS

void
GeneMrnaSim::describereactions()
{

    instate = MatDoub(nn,mm,0.);
    instate[0][0]=1.;
    instate[1][1]=1.;
    instate[0][2]=1.;
    instate[2][3]=1.;
    instate[3][4]=1.;
    instate[4][5]=1.;
    instate[3][6]=1.;
    instate[5][7]=1.;
    instate[6][8]=1.;
    instate[7][9]=1.;
    instate[6][10]=1.;
    instate[8][11]=1.;
    instate[9][12]=1.;
    instate[10][13]=1.;
    instate[9][14]=1.;
    instate[11][15]=1.;
    instate[12][16]=1.;
    instate[13][17]=1.;
    instate[12][18]=1.;
    instate[14][19]=1.;
    instate[15][20]=1.;
    instate[16][21]=1.;
    instate[15][22]=1.;
    instate[17][23]=1.;
    instate[18][24]=1.;
    instate[19][25]=1.;
    instate[18][26]=1.;
    instate[20][27]=1.;
    instate[21][28]=1.;
    instate[22][29]=1.;
    instate[21][30]=1.;
    instate[23][31]=1.;
    instate[24][32]=1.;
    instate[25][33]=1.;
    instate[24][34]=1.;
    instate[26][35]=1.;
    instate[27][36]=1.;
    instate[28][37]=1.;
    instate[27][38]=1.;
    instate[29][39]=1.;
    instate[30][40]=1.;
    instate[31][41]=1.;
    instate[30][42]=1.;
    instate[32][43]=1.;
    instate[33][44]=1.;
    instate[34][45]=1.;
    instate[33][46]=1.;
    instate[35][47]=1.;
    instate[36][48]=1.;
    instate[37][49]=1.;
    instate[36][50]=1.;
    instate[38][51]=1.;
    instate[39][52]=1.;
    instate[40][53]=1.;
    instate[39][54]=1.;
    instate[41][55]=1.;
    instate[42][56]=1.;
    instate[43][57]=1.;
    instate[42][58]=1.;
    instate[44][59]=1.;
    instate[45][60]=1.;
    instate[46][61]=1.;
    instate[45][62]=1.;
    instate[47][63]=1.;
    instate[48][64]=1.;
    instate[49][65]=1.;
    instate[48][66]=1.;
    instate[50][67]=1.;

    outstate = MatDoub(nn,mm,0.);
    outstate[0][0]=-1.;
    outstate[1][0]=1.;
    outstate[0][1]=1.;
    outstate[1][1]=-1.;
    outstate[2][2]=1.;
    outstate[2][3]=-1.;
    outstate[3][4]=-1.;
    outstate[4][4]=1.;
    outstate[3][5]=1.;
    outstate[4][5]=-1.;
    outstate[5][6]=1.;
    outstate[5][7]=-1.;
    outstate[6][8]=-1.;
    outstate[7][8]=1.;
    outstate[6][9]=1.;
    outstate[7][9]=-1.;
    outstate[8][10]=1.;
    outstate[8][11]=-1.;
    outstate[9][12]=-1.;
    outstate[10][12]=1.;
    outstate[9][13]=1.;
    outstate[10][13]=-1.;
    outstate[11][14]=1.;
    outstate[11][15]=-1.;
    outstate[12][16]=-1.;
    outstate[13][16]=1.;
    outstate[12][17]=1.;
    outstate[13][17]=-1.;
    outstate[14][18]=1.;
    outstate[14][19]=-1.;
    outstate[15][20]=-1.;
    outstate[16][20]=1.;
    outstate[15][21]=1.;
    outstate[16][21]=-1.;
    outstate[17][22]=1.;
    outstate[17][23]=-1.;
    outstate[18][24]=-1.;
    outstate[19][24]=1.;
    outstate[18][25]=1.;
    outstate[19][25]=-1.;
    outstate[20][26]=1.;
    outstate[20][27]=-1.;
    outstate[21][28]=-1.;
    outstate[22][28]=1.;
    outstate[21][29]=1.;
    outstate[22][29]=-1.;
    outstate[23][30]=1.;
    outstate[23][31]=-1.;
    outstate[24][32]=-1.;
    outstate[25][32]=1.;
    outstate[24][33]=1.;
    outstate[25][33]=-1.;
    outstate[26][34]=1.;
    outstate[26][35]=-1.;
    outstate[27][36]=-1.;
    outstate[28][36]=1.;
    outstate[27][37]=1.;
    outstate[28][37]=-1.;
    outstate[29][38]=1.;
    outstate[29][39]=-1.;
    outstate[30][40]=-1.;
    outstate[31][40]=1.;
    outstate[30][41]=1.;
    outstate[31][41]=-1.;
    outstate[32][42]=1.;
    outstate[32][43]=-1.;
    outstate[33][44]=-1.;
    outstate[34][44]=1.;
    outstate[33][45]=1.;
    outstate[34][45]=-1.;
    outstate[35][46]=1.;
    outstate[35][47]=-1.;
    outstate[36][48]=-1.;
    outstate[37][48]=1.;
    outstate[36][49]=1.;
    outstate[37][49]=-1.;
    outstate[38][50]=1.;
    outstate[38][51]=-1.;
    outstate[39][52]=-1.;
    outstate[40][52]=1.;
    outstate[39][53]=1.;
    outstate[40][53]=-1.;
    outstate[41][54]=1.;
    outstate[41][55]=-1.;
    outstate[42][56]=-1.;
    outstate[43][56]=1.;
    outstate[42][57]=1.;
    outstate[43][57]=-1.;
    outstate[44][58]=1.;
    outstate[44][59]=-1.;
    outstate[45][60]=-1.;
    outstate[46][60]=1.;
    outstate[45][61]=1.;
    outstate[46][61]=-1.;
    outstate[47][62]=1.;
    outstate[47][63]=-1.;
    outstate[48][64]=-1.;
    outstate[49][64]=1.;
    outstate[48][65]=1.;
    outstate[49][65]=-1.;
    outstate[50][66]=1.;
    outstate[50][67]=-1.;
    dispatch[0] = &GeneMrnaSim::rate0;
    dispatch[1] = &GeneMrnaSim::rate1;
    dispatch[2] = &GeneMrnaSim::rate2;
    dispatch[3] = &GeneMrnaSim::rate3;
    dispatch[4] = &GeneMrnaSim::rate4;
    dispatch[5] = &GeneMrnaSim::rate5;
    dispatch[6] = &GeneMrnaSim::rate6;
    dispatch[7] = &GeneMrnaSim::rate7;
    dispatch[8] = &GeneMrnaSim::rate8;
    dispatch[9] = &GeneMrnaSim::rate9;
    dispatch[10] = &GeneMrnaSim::rate10;
    dispatch[11] = &GeneMrnaSim::rate11;
    dispatch[12] = &GeneMrnaSim::rate12;
    dispatch[13] = &GeneMrnaSim::rate13;
    dispatch[14] = &GeneMrnaSim::rate14;
    dispatch[15] = &GeneMrnaSim::rate15;
    dispatch[16] = &GeneMrnaSim::rate16;
    dispatch[17] = &GeneMrnaSim::rate17;
    dispatch[18] = &GeneMrnaSim::rate18;
    dispatch[19] = &GeneMrnaSim::rate19;
    dispatch[20] = &GeneMrnaSim::rate20;
    dispatch[21] = &GeneMrnaSim::rate21;
    dispatch[22] = &GeneMrnaSim::rate22;
    dispatch[23] = &GeneMrnaSim::rate23;
    dispatch[24] = &GeneMrnaSim::rate24;
    dispatch[25] = &GeneMrnaSim::rate25;
    dispatch[26] = &GeneMrnaSim::rate26;
    dispatch[27] = &GeneMrnaSim::rate27;
    dispatch[28] = &GeneMrnaSim::rate28;
    dispatch[29] = &GeneMrnaSim::rate29;
    dispatch[30] = &GeneMrnaSim::rate30;
    dispatch[31] = &GeneMrnaSim::rate31;
    dispatch[32] = &GeneMrnaSim::rate32;
    dispatch[33] = &GeneMrnaSim::rate33;
    dispatch[34] = &GeneMrnaSim::rate34;
    dispatch[35] = &GeneMrnaSim::rate35;
    dispatch[36] = &GeneMrnaSim::rate36;
    dispatch[37] = &GeneMrnaSim::rate37;
    dispatch[38] = &GeneMrnaSim::rate38;
    dispatch[39] = &GeneMrnaSim::rate39;
    dispatch[40] = &GeneMrnaSim::rate40;
    dispatch[41] = &GeneMrnaSim::rate41;
    dispatch[42] = &GeneMrnaSim::rate42;
    dispatch[43] = &GeneMrnaSim::rate43;
    dispatch[44] = &GeneMrnaSim::rate44;
    dispatch[45] = &GeneMrnaSim::rate45;
    dispatch[46] = &GeneMrnaSim::rate46;
    dispatch[47] = &GeneMrnaSim::rate47;
    dispatch[48] = &GeneMrnaSim::rate48;
    dispatch[49] = &GeneMrnaSim::rate49;
    dispatch[50] = &GeneMrnaSim::rate50;
    dispatch[51] = &GeneMrnaSim::rate51;
    dispatch[52] = &GeneMrnaSim::rate52;
    dispatch[53] = &GeneMrnaSim::rate53;
    dispatch[54] = &GeneMrnaSim::rate54;
    dispatch[55] = &GeneMrnaSim::rate55;
    dispatch[56] = &GeneMrnaSim::rate56;
    dispatch[57] = &GeneMrnaSim::rate57;
    dispatch[58] = &GeneMrnaSim::rate58;
    dispatch[59] = &GeneMrnaSim::rate59;
    dispatch[60] = &GeneMrnaSim::rate60;
    dispatch[61] = &GeneMrnaSim::rate61;
    dispatch[62] = &GeneMrnaSim::rate62;
    dispatch[63] = &GeneMrnaSim::rate63;
    dispatch[64] = &GeneMrnaSim::rate64;
    dispatch[65] = &GeneMrnaSim::rate65;
    dispatch[66] = &GeneMrnaSim::rate66;
    dispatch[67] = &GeneMrnaSim::rate67;
}

