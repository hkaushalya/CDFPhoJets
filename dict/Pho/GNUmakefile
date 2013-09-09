#
#
#
echo:
	@ In dict/Pho 

skip_files := 

LINK_SHARED_MODULES = yes

ifdef USESHLIBS
  SHAREDLIB = libSamantha_Pho.so
else
  LIB = libSamantha_Pho.a
endif


LIBCXXFILES = $(filter-out $(skip_files), $(wildcard $(sharedlib_o_dir)*_dict.cxx))


include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_rootcint.mk




