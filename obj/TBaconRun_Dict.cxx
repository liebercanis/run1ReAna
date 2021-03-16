// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TBaconRun_Dict

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
#include "TBaconRun.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TBaconRun(void *p = 0);
   static void *newArray_TBaconRun(Long_t size, void *p);
   static void delete_TBaconRun(void *p);
   static void deleteArray_TBaconRun(void *p);
   static void destruct_TBaconRun(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TBaconRun*)
   {
      ::TBaconRun *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TBaconRun >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TBaconRun", ::TBaconRun::Class_Version(), "TBaconRun.hxx", 14,
                  typeid(::TBaconRun), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TBaconRun::Dictionary, isa_proxy, 4,
                  sizeof(::TBaconRun) );
      instance.SetNew(&new_TBaconRun);
      instance.SetNewArray(&newArray_TBaconRun);
      instance.SetDelete(&delete_TBaconRun);
      instance.SetDeleteArray(&deleteArray_TBaconRun);
      instance.SetDestructor(&destruct_TBaconRun);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TBaconRun*)
   {
      return GenerateInitInstanceLocal((::TBaconRun*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TBaconRun*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TBaconRun::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TBaconRun::Class_Name()
{
   return "TBaconRun";
}

//______________________________________________________________________________
const char *TBaconRun::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBaconRun*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TBaconRun::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBaconRun*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TBaconRun::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBaconRun*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TBaconRun::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBaconRun*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TBaconRun::Streamer(TBuffer &R__b)
{
   // Stream an object of class TBaconRun.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TBaconRun::Class(),this);
   } else {
      R__b.WriteClassBuffer(TBaconRun::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TBaconRun(void *p) {
      return  p ? new(p) ::TBaconRun : new ::TBaconRun;
   }
   static void *newArray_TBaconRun(Long_t nElements, void *p) {
      return p ? new(p) ::TBaconRun[nElements] : new ::TBaconRun[nElements];
   }
   // Wrapper around operator delete
   static void delete_TBaconRun(void *p) {
      delete ((::TBaconRun*)p);
   }
   static void deleteArray_TBaconRun(void *p) {
      delete [] ((::TBaconRun*)p);
   }
   static void destruct_TBaconRun(void *p) {
      typedef ::TBaconRun current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TBaconRun

namespace {
  void TriggerDictionaryInitialization_TBaconRun_Dict_Impl() {
    static const char* headers[] = {
"TBaconRun.hxx",
0
    };
    static const char* includePaths[] = {
"/usr/local/root/include",
"/.",
"/usr/local/root_v6.12.06/include",
"/Users/gold/baconLocal/git/bacon/run1ReAna/obj/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TBaconRun_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TBaconRun.hxx")))  TBaconRun;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TBaconRun_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TBaconRun.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TBaconRun", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TBaconRun_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TBaconRun_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TBaconRun_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TBaconRun_Dict() {
  TriggerDictionaryInitialization_TBaconRun_Dict_Impl();
}
