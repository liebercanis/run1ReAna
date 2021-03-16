#include "TPmtEvent.hxx"
ClassImp(TPmtEvent)

TPmtEvent::TPmtEvent(): TNamed("TPmtEvent","TPmtEvent")
{
  clear();
}

//TPmtEvent::~TPmtEvent(){}

void TPmtEvent::clear()
{
  event=0;
  time.clear();	 
  volt1.clear();
  volt2.clear();
 }

