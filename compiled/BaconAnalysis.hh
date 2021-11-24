
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
    ETOTAL,
    ESINGLET,
    EDOUBLET,
    ETRIPLET,
    NGOOD,
    OVERSHOOT,
    EMINBASE
};



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

static double singletStart = 0.99;
static double singletEnd = 1.02;
static double afterLow = 1.30;
static double afterHigh = 1.60; //1.60;
static double afterLow2 = 1.88; //1.84;
static double afterHigh2 = 2.1;
//
