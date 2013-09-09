if [ file exists $env(SRT_PRIVATE_CONTEXT)/TclUtils ] { 
  source $env(SRT_PRIVATE_CONTEXT)/TclUtils/scripts/getenv.tcl
} else { 
  source $env(SRT_PUBLIC_CONTEXT)/TclUtils/scripts/getenv.tcl
}


set TheCalibInputFile  [ getenv TheCalibInputFile    "" ]
set lum_log            [ getenv lum_log              "" ]
set NEVENTS            [ getenv NEVENTS              "" ]

#-----------------------------------------------------------------
# ACTIVATE calibration paths
#-----------------------------------------------------------------

set BEAMLINES  1
set TOF        1
set CALORTIME  1
set MINIPLUG   0
set CPR2       1
set CALORMB    1
set PES        1
set VAL        1

#-----------------------------------------------------------------

# singlestream
set singleStream 0

# TwoStream
set twoStream 0

# photon split
set photonSplit 0

# bad runs j split
set jSplit 0

# bad runs h split
set hSplit 0

set rawDataStream rest
set mcProd 0
set writeMetaData 0
set procName PROD
creator set $procName


#  --------------------------------------------------------------------
#  Some initializations.

source  Production/calib/calibexe_init.tcl
source  Production/calib/setup_managers_calib.tcl


#  --------------------------------------------------------------------
# OUTPUTS
# For beamlines
module enable HepRootManager
talk HepRootManager
  histfile set calibhist_$TheCalibInputFile.root   
exit
# For mb-calibration:
talk MbsCalibrationModule
  mbsCaltextFile set mbscal_$TheCalibInputFile.text
exit

#  For ToF (files not active waiting for fixes)

  set TOFTIMEWALK ttw_$TheCalibInputFile.root               
  set TOFLIGHTSPEED tls_$TheCalibInputFile.root                

#  For CalorTime
  set TIMING_NTUPLE hadtimingntuple_$TheCalibInputFile
##  For Miniplug
#  set MP_STNFILE mpstntuple_$TheCalibInputFile.root
#for CPR2
set  CPR2_STNFILE cpr2stntuple_$TheCalibInputFile.root
#  For PES
set PESNTUPLE pesgain_$TheCalibInputFile.root

set MB_STRIPPING mbstrip_$TheCalibInputFile.root

#  --------------------------------------------------------------------
#  Activate calibration ntuplizers

# o PUFF necessary info.
source Production/calib/setup_path_calibexe_dc_puff.tcl


if {$VAL} then {
    set     USING_PRODUCTIONEXE_TCL 0
    source  Production/setup_output.tcl
    source  Production/calib/setup_path_calibexe_dc_mbstrip.tcl
}

# o CALIBRATION: BeamLines
if {$BEAMLINES} then {
    source  Production/calib/setup_path_calibexe_dc_beamlines.tcl  
}

# o CALIBRATION: TOF (IT HANGS AT THE END)
if {$TOF} then {
    source  Production/calib/setup_path_calibexe_dc_tof.tcl    
}

# o CALIBRATION: CALOR TIMMING
if {$CALORTIME} then {
    source  Production/calib/setup_path_calibexe_dc_hadcalortime.tcl  
}

# o CALIBRATION: CPR2
#making stntuples here (be carefull if you want too!-Stntuple modules cannot be cloned)
if {$CPR2} then {
    source  Production/calib/setup_path_calibexe_dc_cpr2.tcl   
}
#not using miniplug anymore!
## o CALIBRATION: MINIPLUG
#if {$MINIPLUG} then {
#    source  Production/calib/setup_path_calibexe_dc_miniplug.tcl   
#}

# o CALIBRATION: MB CALORIMETER
if {$CALORMB} then {
    source  Production/calib/setup_path_calibexe_dc_mbcalo.tcl   
}


# o CALIBRATION: PES GAINS
if {$PES} then {
    source  Production/calib/setup_path_calibexe_dc_pes_gain.tcl  
}



#  --------------------------------------------------------------------

module input DHInput
alias include input
talk DHInput
   dropList add TofPulsesColl
   dropList add TofHitBarColl 
   dropList add TofMatchesColl
   dropList add TofTrackPulseColl
   dropList add TofTrackView
   dropList add TofT0Coll
#   cache set DCACHE
   include file ./$TheCalibInputFile
   luminosityLog set $lum_log 
#   include file root://fcdfdata053.fnal.gov//cdf/scratch/cdfopr/validation_data/6.1.0pre6_calib_maxopt/gphysr_calibexe/gg02f205.0088prod 
   statusFile set 1
exit

set excludeFile  [ getenv ExcludeEventsTcl  nofile ]
if [ file exists $excludeFile ] {
          source $env(ExcludeEventsTcl)
}


if { $NEVENTS <= 0 } then { cont } else { cont -nev $NEVENTS }

show timer
show path
exit
