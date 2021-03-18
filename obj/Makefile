# makefile for run1ReAna 
# MGold, Jan 2016
# not yet ported to OSX!
########################################################
LIBNAME :=run1ReAnaRoot
ARCH  = $(shell root-config --arch)
include $(ROOTSYS)/etc/Makefile.arch


SYS := $(shell uname)
LOCAL        := $(shell pwd)
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
FREE 				 := $(shell freetype-config --libs)

CXXFLAGS     :=  $(ROOTCFLAGS) -I$(LOCAL) 
LIBS         :=  $(ROOTGLIBS) $(CINTFLAGS)


ifeq ($(ARCH),macosx64)
				LIBRARY := lib$(LIBNAME).$(DllSuf)
				LIBLINK := $(subst .$(DllSuf),.so,$(LIBRARY))
				LIBS += $(FREE) 
else
				LIBRARY := lib$(LIBNAME).so
				LIBLINK :=
				CXXFLAGS += -fPIC
endif 



#------------------------------------------------------------------------------
SKIP_OBJ = 
SKIP_DICT = 
OBJS := $(filter-out $(SKIP_OBJ),$(patsubst %.hxx,%.o,$(wildcard *.hxx)))
DICTS := $(filter-out $(SKIP_DICT), $(patsubst %.hxx,%_Dict.cxx,$(wildcard *.hxx)))
#------------------------------------------------------------------------------


all: lib 

lib: $(DICTS) $(OBJS)
	@echo "lib step $(LD) $(SOFLAGS) $(LIBS) $(OBJS) $(DICTS) -o $(LIBRARY) "
	$(LD) $(SOFLAGS) $(CXXFLAGS) $(LIBS) $(OBJS) $(DICTS) -o $(LIBRARY) 
	@echo "$@ done"


CINTFLAGS   := -I$(ROOTSYS)/include -I/. 
#-DDebug rootcling?
%_Dict.cxx: %.hxx %_LinkDef.h
	@echo "rootcint -f $@ -c $(CINTCFLAGS) $^"
	@rootcint -f $@ -c $(CINTFLAGS) $^

	
%.o:  %.cxx %.cxx
	@echo "compiling $@ ... $(CXX) $(CXXFLAGS) -c $<"
	$(CXX) $(CXXFLAGS) -c $<


.PHONY : clean show

show: 
	@echo sys =$(SYS)
	@echo arch =$(ARCH)
	@echo OBJS = $(OBJS)
	@echo DICTS = $(DICTS)
	@echo LIBS = $(LIBS)
	@echo LDFLAGS= $(LDFLAGS)
	@echo ROOTCFLAGS= $(ROOTCFLAGS)
	@echo CINTFLAGS= $(CINTFLAGS)
	@echo CXXFLAGS= $(CXXFLAGS)
	@echo EXTRA_CFLAGS= $(EXTRA_CFLAGS)
	@echo DllSuf= $(DllSuf)
	@echo OutPutOpt= $(OutPutOpt)
	@echo CAPTROOTCONFIG=$(CAPTROOTCONFIG)	
	@echo LIBRARY= $(LIBRARY)
	@echo LIBLINK $(LIBLINK)
	@echo LIBNAME=$(LIBNAME)	
	@echo CXX=$(CXX)

clean:
	@echo "removing $(OBJS) *_Dict.* $(LIBRARY) *.pcm"
	@rm -f $(OBJS) *_Dict.* $(LIBRARY) *.pcm 
