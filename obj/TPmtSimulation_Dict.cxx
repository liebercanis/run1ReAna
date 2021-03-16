// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TPmtSimulation_Dict

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
#include "TPmtSimulation.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TPmtSimulation(void *p = 0);
   static void *newArray_TPmtSimulation(Long_t size, void *p);
   static void delete_TPmtSimulation(void *p);
   static void deleteArray_TPmtSimulation(void *p);
   static void destruct_TPmtSimulation(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TPmtSimulation*)
   {
      ::TPmtSimulation *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TPmtSimulation >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TPmtSimulation", ::TPmtSimulation::Class_Version(), "TPmtSimulation.hxx", 15,
                  typeid(::TPmtSimulation), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TPmtSimulation::Dictionary, isa_proxy, 4,
                  sizeof(::TPmtSimulation) );
      instance.SetNew(&new_TPmtSimulation);
      instance.SetNewArray(&newArray_TPmtSimulation);
      instance.SetDelete(&delete_TPmtSimulation);
      instance.SetDeleteArray(&deleteArray_TPmtSimulation);
      instance.SetDestructor(&destruct_TPmtSimulation);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TPmtSimulation*)
   {
      return GenerateInitInstanceLocal((::TPmtSimulation*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TPmtSimulation*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TPmtSimulation::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TPmtSimulation::Class_Name()
{
   return "TPmtSimulation";
}

//______________________________________________________________________________
const char *TPmtSimulation::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPmtSimulation*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TPmtSimulation::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPmtSimulation*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TPmtSimulation::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPmtSimulation*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TPmtSimulation::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPmtSimulation*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TPmtSimulation::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPmtSimulation.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TPmtSimulation::Class(),this);
   } else {
      R__b.WriteClassBuffer(TPmtSimulation::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TPmtSimulation(void *p) {
      return  p ? new(p) ::TPmtSimulation : new ::TPmtSimulation;
   }
   static void *newArray_TPmtSimulation(Long_t nElements, void *p) {
      return p ? new(p) ::TPmtSimulation[nElements] : new ::TPmtSimulation[nElements];
   }
   // Wrapper around operator delete
   static void delete_TPmtSimulation(void *p) {
      delete ((::TPmtSimulation*)p);
   }
   static void deleteArray_TPmtSimulation(void *p) {
      delete [] ((::TPmtSimulation*)p);
   }
   static void destruct_TPmtSimulation(void *p) {
      typedef ::TPmtSimulation current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TPmtSimulation

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 210,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace {
  void TriggerDictionaryInitialization_TPmtSimulation_Dict_Impl() {
    static const char* headers[] = {
"TPmtSimulation.hxx",
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
#line 1 "TPmtSimulation_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$TPmtSimulation.hxx")))  TPmtSimulation;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TPmtSimulation_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TPmtSimulation.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TPmtSimulation", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TPmtSimulation_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TPmtSimulation_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TPmtSimulation_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TPmtSimulation_Dict() {
  TriggerDictionaryInitialization_TPmtSimulation_Dict_Impl();
}
