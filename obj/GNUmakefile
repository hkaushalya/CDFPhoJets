#
#
#

CXXFLAGS += -I/usr/local/include/boost-1_34_1

#TBINS = 
SIMPLEBINS = $(filter-out $(COMPLEXBIN),$(BINS) $(TBINS))

skip_files := 

LINK_SHARED_MODULES = yes
ifdef USESHLIBS
  SHAREDLIB = libSamantha_obj.so
else
  LIB = libSamantha_obj.a
endif
 

LIBCCFILES = $(filter-out $(skip_files), $(wildcard  *.cc))


include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
