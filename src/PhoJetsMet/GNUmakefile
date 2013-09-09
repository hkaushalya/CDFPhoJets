#
#
#


#TBINS = 
SIMPLEBINS = $(filter-out $(COMPLEXBIN),$(BINS) $(TBINS))

skip_files := 

LINK_SHARED_MODULES = yes
ifdef USESHLIBS
  SHAREDLIB = libSamantha_PhoJetsMet.so
else
  LIB = libSamantha_PhoJetsMet.a
endif
 

LIBCCFILES = $(filter-out $(skip_files), $(wildcard  *.cc))


include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk





