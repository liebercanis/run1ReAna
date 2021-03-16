// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TBaconEvent_Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "TBaconEvent.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TBaconEvent(void *p = 0);
   static void *newArray_TBaconEvent(Long_t size, void *p);
   static void delete_TBaconEvent(void *p);
   static void deleteArray_TBaconEvent(void *p);
   static void destruct_TBaconEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TBaconEvent*)
   {
      ::TBaconEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TBaconEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TBaconEvent", ::TBaconEvent::Class_Version(), "TBaconEvent.hxx", 16,
                  typeid(::TBaconEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TBaconEvent::Dictionary, isa_proxy, 4,
                  sizeof(::TBaconEvent) );
      instance.SetNew(&new_TBaconEvent);
      instance.SetNewArray(&newArray_TBaconEvent);
      instance.SetDelete(&delete_TBaconEvent);
      instance.SetDeleteArray(&deleteArray_TBaconEvent);
      instance.SetDestructor(&destruct_TBaconEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TBaconEvent*)
   {
      return GenerateInitInstanceLocal((::TBaconEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TBaconEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TBaconEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TBaconEvent::Class_Name()
{
   return "TBaconEvent";
}

//______________________________________________________________________________
const char *TBaconEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBaconEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TBaconEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBaconEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TBaconEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBaconEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TBaconEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBaconEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TBaconEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class TBaconEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TBaconEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(TBaconEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TBaconEvent(void *p) {
      return  p ? new(p) ::TBaconEvent : new ::TBaconEvent;
   }
   static void *newArray_TBaconEvent(Long_t nElements, void *p) {
      return p ? new(p) ::TBaconEvent[nElements] : new ::TBaconEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_TBaconEvent(void *p) {
      delete ((::TBaconEvent*)p);
   }
   static void deleteArray_TBaconEvent(void *p) {
      delete [] ((::TBaconEvent*)p);
   }
   static void destruct_TBaconEvent(void *p) {
      typedef ::TBaconEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TBaconEvent

namespace ROOT {
   static TClass *vectorlETPulsegR_Dictionary();
   static void vectorlETPulsegR_TClassManip(TClass*);
   static void *new_vectorlETPulsegR(void *p = 0);
   static void *newArray_vectorlETPulsegR(Long_t size, void *p);
   static void delete_vectorlETPulsegR(void *p);
   static void deleteArray_vectorlETPulsegR(void *p);
   static void destruct_vectorlETPulsegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TPulse>*)
   {
      vector<TPulse> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TPulse>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TPulse>", -2, "vector", 210,
                  typeid(vector<TPulse>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETPulsegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TPulse>) );
      instance.SetNew(&new_vectorlETPulsegR);
      instance.SetNewArray(&newArray_vectorlETPulsegR);
      instance.SetDelete(&delete_vectorlETPulsegR);
      instance.SetDeleteArray(&deleteArray_vectorlETPulsegR);
      instance.SetDestructor(&destruct_vectorlETPulsegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TPulse> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TPulse>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETPulsegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TPulse>*)0x0)->GetClass();
      vectorlETPulsegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETPulsegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETPulsegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TPulse> : new vector<TPulse>;
   }
   static void *newArray_vectorlETPulsegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TPulse>[nElements] : new vector<TPulse>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETPulsegR(void *p) {
      delete ((vector<TPulse>*)p);
   }
   static void deleteArray_vectorlETPulsegR(void *p) {
      delete [] ((vector<TPulse>*)p);
   }
   static void destruct_vectorlETPulsegR(void *p) {
      typedef vector<TPulse> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TPulse>

namespace {
  void TriggerDictionaryInitialization_TBaconEvent_Dict_Impl() {
    static const char* headers[] = {
"TBaconEvent.hxx",
0
    };
    static const char* includePaths[] = {
"/usr/local/root/include",
"/.",
"/home/admin/root-6.14.06/include",
"/home/gold/bacon/pmtLocal/obj/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TBaconEvent_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TPulse.hxx")))  __attribute__((annotate("$clingAutoload$TBaconEvent.hxx")))  TPulse;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$TBaconEvent.hxx")))  TBaconEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TBaconEvent_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TBaconEvent.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TBaconEvent", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TBaconEvent_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TBaconEvent_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TBaconEvent_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TBaconEvent_Dict() {
  TriggerDictionaryInitialization_TBaconEvent_Dict_Impl();
}
