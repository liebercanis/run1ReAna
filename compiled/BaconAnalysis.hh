enum
{
    SUMVARS = 11
};

enum
{
    ERUN,
    ESET,
    EBASE,
    EBASEEND,
    EACCEPT,
    ETOTAL,
    ESINGLET,
    EDOUBLET,
    ETRIPLET,
    NGOOD,
    OVERSHOOT
};

enum
{
    EVVARS = 10
};

enum
{
    EVRUN,
    EVSET,
    EVFLAG,
    EVSUM,
    EVSINGLET,
    EVTRIPLET,
    EVLATE,
    EVLATETIME,
    EVWFSINGLET,
    EVWFMIN,
};
static double singletStart = 0.99;
static double singletEnd = 1.1;
static double afterLow = 1.30;
static double afterHigh = 1.60; //1.60;
static double afterLow2 = 1.80; //1.84;
static double afterHigh2 = 2.2;
//
