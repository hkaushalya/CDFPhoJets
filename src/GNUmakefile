#############################################################
#
# pass to subdirectories
#
#############################################################

sharedlib: codegen
bin: sharedlib

.PHONY: pho

override SRT_QUAL := debug

SUBDIRS = Pho

pho:
	$(MAKE) -C dict/Pho $(OVERRIDES) codegen 
	$(MAKE) -C dict/Pho $(OVERRIDES) sharedlib
	$(MAKE) -C Pho $(OVERRIDES) sharedlib
	
include SoftRelTools/standard.mk

