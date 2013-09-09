# mc_Pythia_photon_jl1.tcl
#
module enable Pythia
talk Pythia
   PythiaMenu
# set process to none, override with msub below
     msel set 0
     cmEnergy set 1960
     inListLevel set 0
     evListLevel set 1
# put the Z onshell
     commonMenu
# diphotons
       set_msub -index=14  -value=1
       set_msub -index=29  -value=1
       set_msub -index=115 -value=1
#      set q**2 scale to -t_hat
       set_mstp -index=32 -value=5
#      set minimum pt of hard process in GeV, approx=min photon pt
       set_ckin -index=3 -value=8
#      set min and max rapidities for lower- and higher-eta pho and jet
       set_ckin -index=13 -value=-2.0
       set_ckin -index=14 -value=2.0
       set_ckin -index=15 -value=-2.0
       set_ckin -index=16 -value=2.0
     exit
#     decayFileMode set read
#     decayFile set $env(PYTHIA_DIR)/z_ee.dcy
   exit
exit
#-----------------------------------------------------------------------
# Pythia tunings
#-----------------------------------------------------------------------
source $env(PROJECT_DIR)/mcProduction/tcl/mc_Pythia_WZPt_tune.tcl
source $env(PROJECT_DIR)/mcProduction/tcl/mc_Pythia_pdf_CTEQ5L.tcl
source $env(PROJECT_DIR)/mcProduction/tcl/mc_Pythia_underlying_event_A.tcl
