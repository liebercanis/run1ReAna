// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TPmtRun_Dict

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
#include "TPmtRun.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TPmtRun(void *p = 0);
   static void *newArray_TPmtRun(Long_t size, void *p);
   static void delete_TPmtRun(void *p);
   static void deleteArray_TPmtRun(void *p);
   static void destruct_TPmtRun(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TPmtRun*)
   {
      ::TPmtRun *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TPmtRun >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TPmtRun", ::TPmtRun::Class_Version(), "TPmtRun.hxx", 16,
                  typeid(::TPmtRun), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TPmtRun::Dictionary, isa_proxy, 4,
                  sizeof(::TPmtRun) );
      instance.SetNew(&new_TPmtRun);
      instance.SetNewArray(&newArray_TPmtRun);
      instance.SetDelete(&delete_TPmtRun);
      instance.SetDeleteArray(&deleteArray_TPmtRun);
      instance.SetDestructor(&destruct_TPmtRun);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TPmtRun*)
   {
      return GenerateInitInstanceLocal((::TPmtRun*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TPmtRun*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TPmtRun::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TPmtRun::Class_Name()
{
   return "TPmtRun";
}

//______________________________________________________________________________
const char *TPmtRun::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPmtRun*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TPmtRun::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPmtRun*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TPmtRun::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPmtRun*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TPmtRun::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPmtRun*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TPmtRun::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPmtRun.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TPmtRun::Class(),this);
   } else {
      R__b.WriteClassBuffer(TPmtRun::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TPmtRun(void *p) {
      return  p ? new(p) ::TPmtRun : new ::TPmtRun;
   }
   static void *newArray_TPmtRun(Long_t nElements, void *p) {
      return p ? new(p) ::TPmtRun[nElements] : new ::TPmtRun[nElements];
   }
   // Wrapper around operator delete
   static void delete_TPmtRun(void *p) {
      delete ((::TPmtRun*)p);
   }
   static void deleteArray_TPmtRun(void *p) {
      delete [] ((::TPmtRun*)p);
   }
   static void destruct_TPmtRun(void *p) {
      typedef ::TPmtRun current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TPmtRun

namespace ROOT {
   static TClass *vectorlETPmtEventgR_Dictionary();
   static void vectorlETPmtEventgR_TClassManip(TClass*);
   static void *new_vectorlETPmtEventgR(void *p = 0);
   static void *newArray_vectorlETPmtEventgR(Long_t size, void *p);
   static void delete_vectorlETPmtEventgR(void *p);
   static void deleteArray_vectorlETPmtEventgR(void *p);
   static void destruct_vectorlETPmtEventgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TPmtEvent>*)
   {
      vector<TPmtEvent> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TPmtEvent>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TPmtEvent>", -2, "vector", 210,
                  typeid(vector<TPmtEvent>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETPmtEventgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TPmtEvent>) );
      instance.SetNew(&new_vectorlETPmtEventgR);
      instance.SetNewArray(&newArray_vectorlETPmtEventgR);
      instance.SetDelete(&delete_vectorlETPmtEventgR);
      instance.SetDeleteArray(&deleteArray_vectorlETPmtEventgR);
      instance.SetDestructor(&destruct_vectorlETPmtEventgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TPmtEvent> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TPmtEvent>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETPmtEventgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TPmtEvent>*)0x0)->GetClass();
      vectorlETPmtEventgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETPmtEventgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETPmtEventgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TPmtEvent> : new vector<TPmtEvent>;
   }
   static void *newArray_vectorlETPmtEventgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TPmtEvent>[nElements] : new vector<TPmtEvent>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETPmtEventgR(void *p) {
      delete ((vector<TPmtEvent>*)p);
   }
   static void deleteArray_vectorlETPmtEventgR(void *p) {
      delete [] ((vector<TPmtEvent>*)p);
   }
   static void destruct_vectorlETPmtEventgR(void *p) {
      typedef vector<TPmtEvent> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TPmtEvent>

namespace {
  void TriggerDictionaryInitialization_TPmtRun_Dict_Impl() {
    static const char* headers[] = {
"TPmtRun.hxx",
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
#line 1 "TPmtRun_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TPmtEvent.hxx")))  __attribute__((annotate("$clingAutoload$TPmtRun.hxx")))  TPmtEvent;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$TPmtRun.hxx")))  TPmtRun;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TPmtRun_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TPmtRun.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TPmtRun", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TPmtRun_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TPmtRun_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TPmtRun_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TPmtRun_Dict() {
  TriggerDictionaryInitialization_TPmtRun_Dict_Impl();
}
