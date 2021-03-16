#include "TPmtSimMatchStats.hxx"
ClassImp(TPmtSimMatchStats)

  TPmtSimMatchStats::TPmtSimMatchStats(): TNamed("TPmtSimMatchStats","TPmtSimMatchStats")
{
  clear();
}

//TPmtSimMatchStats::~TPmtSimMatchStats(){}

void TPmtSimMatchStats::clear()
{
  entries=0;
  tSim=0;
  tHit=0;
  tMatch=0;
  tNot=0;
  tMiss=0;
  eff=0;
  effError=0;
  over=0;
  overError=0;
}

void TPmtSimMatchStats::fill(unsigned nsim, unsigned nhit, unsigned nmatch, unsigned nnot, unsigned nmiss)
{
  ++entries;
  tSim += nsim;
  tHit += nhit;
  tMatch += nmatch;
  tNot += nnot;
  tMiss += nmiss;
  eff   = double(tMatch)/double(tSim);
  over = double(tNot)/double(entries);
}

void TPmtSimMatchStats::print() 
{
  printf(" \t sim match stats: events %i sim %i hit %i match %i miss %i not %i efficiency/hit %f  over hits/event %f \n",
      entries,tSim,tHit,tMatch,tMiss,tNot,eff,over);
}
