#
#
#

skip_files := 

LINK_SHARED_MODULES = yes

ifdef USESHLIBS
  SHAREDLIB = libSamantha_MetModel.so
else
  LIB = libSamantha_MetModel.a
endif


LIBCXXFILES = $(filter-out $(skip_files), $(wildcard $(sharedlib_o_dir)*_dict.cxx))


include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_rootcint.mk




