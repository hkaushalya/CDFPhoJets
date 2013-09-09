# mc_Pythia_photon_sam.tcl : this has most setting similar to mc_Pythia_photon_20cen.tcl
# this is for Photon+jets production. to study the ISR/FSR effects
# at HEPG level
########################################################################################3

# deafult =0, more 1, less -1
set ISR       [ getenv ISR 0 ]
set FSR       [ getenv FSR 0 ]
set Q2        [ getenv Q2  0 ]
set minPtHat  [ getenv minPtHat 20.0 ]

puts "mc_Pythia_photon_sam: min ptHat  = $minPtHat"
puts "mc_Pythia_photon_sam:legend: 0==default, 1==more/up, -1==less/down"
puts "mc_Pythia_photon_sam: Q2         = $Q2"
puts "mc_Pythia_photon_sam: ISR        = $ISR"
puts "mc_Pythia_photon_sam: FSR        = $FSR"

module enable Pythia
talk Pythia
   PythiaMenu
# set process to none, override with msub below
     msel set 0
     cmEnergy set 1960
     inListLevel set 0
     evListLevel set 1

     commonMenu

#      photon+jets
       set_msub -index=14  -value=1
       set_msub -index=29  -value=1
       set_msub -index=115 -value=1
#      set q**2 scale to -t_hat
       set_mstp -index=32 -value=5

#      modify q**2 scale by this factor to derive q**2 uncertainty
	if { $Q2 == 0 } {
		puts "mc_Pythia_photon_sam: Q2 is NOT scaled. Using PYTHIA default."
	}

	if { $Q2 == 1 } {
		puts "mc_Pythia_photon_sam: Using scaled up Q2"
		 set_parp -index=34 -value=2.0
	}
	if { $Q2 == -1 } {
		puts "mc_Pythia_photon_sam: Using scaling down Q2"
		 set_parp -index=34 -value=0.5
	}
		 
		 
#      set minimum pt of hard process in GeV, approx=min photon pt
#      All the matrix elements in this group, MSEL=1,2 , are for massless quarks
#      (although  final state quarks are of course put on the mass shell). As a consequence,
#      cross sections are divergent for P(perp)->0 (tranverse momentum),
#		 and some kind of regularization is required.
#      This is done by setting the P(perp)-min value in CKIN(3)
#		 To generate hard scatterings with 50 GeV<Pt<100 GeV, for instance, use 
#      CKIN(3)=50D0 
#      CKIN(4)=100D0 

       set_ckin -index=3 -value=$minPtHat

#      set min and max rapidities for lower- and higher-eta pho and jet
       set_ckin -index=13 -value=-1.1
       set_ckin -index=14 -value=1.1
       set_ckin -index=15 -value=-3.3
       set_ckin -index=16 -value=3.3
 
#----------------------------------------
#this must be set for ISR/FSR -> MSTP(3)=1 
#----------------------------------------
 		set_mstp -index=3 -value=1 

#-------- controls the ISR --------------
#from Joint Physics group recommnedations
#-----------------------------------------
#DEFAULT
	if { $ISR == 0 } {
       puts "mc_Pythia_photon_sam: Using DEFAULT ISR values." 
		 set_parp -index=61 -value=0.146
		 set_parp -index=64 -value=1.0
	}
#MORE ISR
	if { $ISR == 1 } {
       puts "mc_Pythia_photon_sam: Using MORE ISR values." 
		 set_parp -index=61 -value=0.292
		 set_parp -index=64 -value=0.5
	}
#LESS ISR
	if { $ISR == -1 } {
       puts "mc_Pythia_photon_sam: Using LESS ISR values." 
		 set_parp -index=61 -value=0.073
		 set_parp -index=64 -value=2.0
	}		 
#-------- controls the FSR --------------
#from Joint Physics group recommnedations
#-----------------------------------------
#DEFAULT
	if { $FSR == 0 } {
       puts "mc_Pythia_photon_sam: Using DEFAULT FSR values." 
		 set_parp -index=72 -value=0.146
		 set_parp -index=71 -value=4.0
	}
#MORE FSR
	if { $FSR == 1 } {
       puts "mc_Pythia_photon_sam: Using MORE FSR values." 
		 set_parp -index=72 -value=0.292
		 set_parp -index=71 -value=8.0
	}
#LESS FSR
	if { $FSR == -1 } {
       puts "mc_Pythia_photon_sam: Using LESS FSR values." 
		set_parp -index=72 -value=0.073
		set_parp -index=71 -value=2.0
	}

		 
     exit
	  
   exit
exit
#-----------------------------------------------------------------------
# Pythia tunes: pdf_CTEQ5L, undelying_event_A
#-----------------------------------------------------------------------
source $env(PROJECT_DIR)/mcProduction/tcl/mc_Pythia_pdf_CTEQ5L.tcl
source $env(PROJECT_DIR)/mcProduction/tcl/mc_Pythia_underlying_event_A.tcl
