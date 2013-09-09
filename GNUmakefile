#############################################################
#
# pass to subdirectories
#
#############################################################

sharedlib: codegen
bin: sharedlib

#.PHONY: obj utils Pho MetModel 
.PHONY: obj utils pho metmodel 

override SRT_QUAL := debug

SUBDIRS = obj src utils 


obj:
	@echo compiling obj
	$(MAKE) -C dict/obj $(OVERRIDES) codegen 
	$(MAKE) -C dict/obj $(OVERRIDES) sharedlib
	$(MAKE) -C obj $(OVERRIDES) sharedlib
utils:
	@echo compiling utils
	$(MAKE) -C dict/utils $(OVERRIDES) codegen 
	$(MAKE) -C dict/utils $(OVERRIDES) sharedlib
	$(MAKE) -C utils/ $(OVERRIDES) sharedlib
pho:
	@echo compiling Pho
	$(MAKE) -C dict/Pho $(OVERRIDES) codegen 
	$(MAKE) -C dict/Pho $(OVERRIDES) sharedlib
	$(MAKE) -C src/Pho $(OVERRIDES) sharedlib
metmodel:
	@echo compiling MetModel
	$(MAKE) -C dict/MetModel $(OVERRIDES) codegen
	$(MAKE) -C dict/MetModel $(OVERRIDES) sharedlib
	$(MAKE) -C src/MetModel $(OVERRIDES) sharedlib


include SoftRelTools/standard.mk

