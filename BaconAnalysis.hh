static double sStart = 900.;
static double sEnd = 1020.;
static double tStart = 1020;
static double tTripletCut = 2300.;
static double tEnd = 10000;
static double tailStart = 3000.0;
static double tailStop = 10000.0;
static int nDigi = 10000;

static double singletFraction = 0.14;
//static double tres = 7.86;
//static double tTriplet = 1600.0; //2100.0;
//static double tSinglet = 5.0;
//static double tXe = 20.0;
//static double kxe = 8.8E-5;
//static double kplusZero = 1.3E-4;
//static double xTrigger = 1000;
//double k1Zero = kxe * 131. / 40.;

//
enum
{
    SUMVARS = 12
};

enum
{
    ERUN,
    ESET,
    EBASE,
    EBASEEND,
    EACCEPT,
    ESINGLET,
    EDOUBLET,
    ETRIPLET,
    NTOTAL,
    NGOOD,
    OVERSHOOT,
    EMINBASE
};
Float_t sumVars[SUMVARS];
TString sumNames[SUMVARS];
TString setSumNames() {
    sumNames[ERUN] = TString("run");
    sumNames[ESET] = TString("set");
    sumNames[EBASE] = TString("base");
    sumNames[EBASEEND] = TString("baseend");
    sumNames[EACCEPT] = TString("accept");
    sumNames[ESINGLET] = TString("singlet");
    sumNames[EDOUBLET] = TString("dublet");
    sumNames[ETRIPLET] = TString("triplet");
    sumNames[NTOTAL] = TString("ntotal");
    sumNames[NGOOD] = TString("ngood");
    sumNames[OVERSHOOT] = TString("over");
    sumNames[EMINBASE] = TString("minbase");

    TString sumString;
    for (int is = 0; is < SUMVARS; ++is)
    {
        if (is < SUMVARS - 1)
            sumString += sumNames[is] + TString(":");
        else
            sumString += sumNames[is];
    }
    return sumString;
}
enum
{
    EVVARS = 11
};

enum
{
    EVEVENT,
    EVRUN,
    EVSET,
    EVFLAG,
    EVSUM,
    EVSINGLET,
    EVTRIPLET,
    EVLATE,
    EVLATETIME,
    EVWFSINGLET,
    EVWFMIN
};
float evVars[EVVARS];
TString evNames[EVVARS];
TString setEvNames()
{
    evNames[0] = TString("ev");
    evNames[1] = TString("run");
    evNames[2] = TString("set");
    evNames[3] = TString("flag");
    evNames[4] = TString("sum");
    evNames[5] = TString("singlet");
    evNames[6] = TString("triplet");
    evNames[7] = TString("late");
    evNames[8] = TString("latetime");
    evNames[9] = TString("wfsinglet");
    evNames[10] = TString("wfmin");
    TString evString;
    for (int is = 0; is < EVVARS; ++is)
    {
        if (is < EVVARS - 1)
            evString += evNames[is] + TString(":");
        else
            evString += evNames[is];
    }
    return evString;
}

enum
{
    PREVARS = 9
};
enum
{
    PREEVENT,
    PRERUN,
    PRESET,
    PREFLAG,
    PRESUM,
    PRESINGLET,
    PRETRIPLET,
    PRELATE,
    PRELATETIME
};

float preVars[PREVARS];
TString preNames[PREVARS];
TString setPreNames()
{
    preNames[0] = TString("ev");
    preNames[1] = TString("run");
    preNames[2] = TString("set");
    preNames[3] = TString("flag");
    preNames[4] = TString("sum");
    preNames[5] = TString("singlet");
    preNames[6] = TString("triplet");
    preNames[7] = TString("late");
    preNames[8] = TString("latetime");
    TString preString;
    for (int is = 0; is < PREVARS; ++is)
    {
        if (is < PREVARS - 1)
            preString += preNames[is] + TString(":");
        else
            preString += preNames[is];
    }
    return preString;
}
