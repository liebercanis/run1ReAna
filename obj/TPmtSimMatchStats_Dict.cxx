// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TPmtSimMatchStats_Dict

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
#include "TPmtSimMatchStats.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TPmtSimMatchStats(void *p = 0);
   static void *newArray_TPmtSimMatchStats(Long_t size, void *p);
   static void delete_TPmtSimMatchStats(void *p);
   static void deleteArray_TPmtSimMatchStats(void *p);
   static void destruct_TPmtSimMatchStats(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TPmtSimMatchStats*)
   {
      ::TPmtSimMatchStats *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TPmtSimMatchStats >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TPmtSimMatchStats", ::TPmtSimMatchStats::Class_Version(), "TPmtSimMatchStats.hxx", 15,
                  typeid(::TPmtSimMatchStats), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TPmtSimMatchStats::Dictionary, isa_proxy, 4,
                  sizeof(::TPmtSimMatchStats) );
      instance.SetNew(&new_TPmtSimMatchStats);
      instance.SetNewArray(&newArray_TPmtSimMatchStats);
      instance.SetDelete(&delete_TPmtSimMatchStats);
      instance.SetDeleteArray(&deleteArray_TPmtSimMatchStats);
      instance.SetDestructor(&destruct_TPmtSimMatchStats);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TPmtSimMatchStats*)
   {
      return GenerateInitInstanceLocal((::TPmtSimMatchStats*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TPmtSimMatchStats*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TPmtSimMatchStats::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TPmtSimMatchStats::Class_Name()
{
   return "TPmtSimMatchStats";
}

//______________________________________________________________________________
const char *TPmtSimMatchStats::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPmtSimMatchStats*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TPmtSimMatchStats::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPmtSimMatchStats*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TPmtSimMatchStats::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPmtSimMatchStats*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TPmtSimMatchStats::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPmtSimMatchStats*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TPmtSimMatchStats::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPmtSimMatchStats.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TPmtSimMatchStats::Class(),this);
   } else {
      R__b.WriteClassBuffer(TPmtSimMatchStats::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TPmtSimMatchStats(void *p) {
      return  p ? new(p) ::TPmtSimMatchStats : new ::TPmtSimMatchStats;
   }
   static void *newArray_TPmtSimMatchStats(Long_t nElements, void *p) {
      return p ? new(p) ::TPmtSimMatchStats[nElements] : new ::TPmtSimMatchStats[nElements];
   }
   // Wrapper around operator delete
   static void delete_TPmtSimMatchStats(void *p) {
      delete ((::TPmtSimMatchStats*)p);
   }
   static void deleteArray_TPmtSimMatchStats(void *p) {
      delete [] ((::TPmtSimMatchStats*)p);
   }
   static void destruct_TPmtSimMatchStats(void *p) {
      typedef ::TPmtSimMatchStats current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TPmtSimMatchStats

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 210,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace {
  void TriggerDictionaryInitialization_TPmtSimMatchStats_Dict_Impl() {
    static const char* headers[] = {
"TPmtSimMatchStats.hxx",
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
#line 1 "TPmtSimMatchStats_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$TPmtSimMatchStats.hxx")))  TPmtSimMatchStats;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TPmtSimMatchStats_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TPmtSimMatchStats.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TPmtSimMatchStats", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TPmtSimMatchStats_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TPmtSimMatchStats_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TPmtSimMatchStats_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TPmtSimMatchStats_Dict() {
  TriggerDictionaryInitialization_TPmtSimMatchStats_Dict_Impl();
}
