#
#
#

CXXFLAGS += -I/usr/local/include/boost-1_34_1
CXXFLAGS += -Wall

#TBINS = 
SIMPLEBINS = $(filter-out $(COMPLEXBIN),$(BINS) $(TBINS))

skip_files := 

LINK_SHARED_MODULES = yes
ifdef USESHLIBS
  SHAREDLIB = libSamantha_Pho.so
else
  LIB = libSamantha_Pho.a
endif
 

LIBCCFILES = $(filter-out $(skip_files), $(wildcard  *.cc))


include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk





