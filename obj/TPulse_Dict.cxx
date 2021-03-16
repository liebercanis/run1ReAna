// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TPulse_Dict

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
#include "TPulse.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TPulse(void *p = 0);
   static void *newArray_TPulse(Long_t size, void *p);
   static void delete_TPulse(void *p);
   static void deleteArray_TPulse(void *p);
   static void destruct_TPulse(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TPulse*)
   {
      ::TPulse *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TPulse >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TPulse", ::TPulse::Class_Version(), "TPulse.hxx", 15,
                  typeid(::TPulse), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TPulse::Dictionary, isa_proxy, 4,
                  sizeof(::TPulse) );
      instance.SetNew(&new_TPulse);
      instance.SetNewArray(&newArray_TPulse);
      instance.SetDelete(&delete_TPulse);
      instance.SetDeleteArray(&deleteArray_TPulse);
      instance.SetDestructor(&destruct_TPulse);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TPulse*)
   {
      return GenerateInitInstanceLocal((::TPulse*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TPulse*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TPulse::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TPulse::Class_Name()
{
   return "TPulse";
}

//______________________________________________________________________________
const char *TPulse::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPulse*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TPulse::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPulse*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TPulse::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPulse*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TPulse::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPulse*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TPulse::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPulse.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TPulse::Class(),this);
   } else {
      R__b.WriteClassBuffer(TPulse::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TPulse(void *p) {
      return  p ? new(p) ::TPulse : new ::TPulse;
   }
   static void *newArray_TPulse(Long_t nElements, void *p) {
      return p ? new(p) ::TPulse[nElements] : new ::TPulse[nElements];
   }
   // Wrapper around operator delete
   static void delete_TPulse(void *p) {
      delete ((::TPulse*)p);
   }
   static void deleteArray_TPulse(void *p) {
      delete [] ((::TPulse*)p);
   }
   static void destruct_TPulse(void *p) {
      typedef ::TPulse current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TPulse

namespace {
  void TriggerDictionaryInitialization_TPulse_Dict_Impl() {
    static const char* headers[] = {
"TPulse.hxx",
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
#line 1 "TPulse_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TPulse.hxx")))  TPulse;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TPulse_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TPulse.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TPulse", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TPulse_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TPulse_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TPulse_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TPulse_Dict() {
  TriggerDictionaryInitialization_TPulse_Dict_Impl();
}
