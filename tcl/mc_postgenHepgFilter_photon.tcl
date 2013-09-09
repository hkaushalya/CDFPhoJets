# mc_postgenHepgFilter_photon_jl1.tcl
#
#This sets up the general hepg filter

# First you need to set up HistogramManager because CP
# insists on throwing an exception if there's no Histogram manager.
# I.e. unless you have HepRootManager or alike in your path BEFORE
# any of the GenTrig stuff your exe will crash at run time.

if [file exists HepHist.dat]  { exec rm HepHist.dat }
mod enable HepRootManager
mod talk HepRootManager
   histfile set HepHist.dat
   createHistoFile set false
exit

# Second you set up the general hepg filter:

#mod enable GenTrigBFilter
#mod talk GenTrigBFilter
#  GenTrigBFilter
#    AbsPDG set 1
#    AbsPDG set 1
#    CodePDG set 22
#    PtMin set 10
#    EtaMin set -3.0
#    EtaMax set 3.0
#    show
#  exit
#exit

mod enable PartFilter
mod talk PartFilter
  Particle1
    idhepCode1 set 22
    Et1  set 22.0
    Eta1 set 1.1
    show
  exit
exit




