# this is to make HEPG only stntuple format job output.
#------------------------------------------------------------------------------
# save only hepg in the output file
#------------------------------------------------------------------------------

set MC_JOB 1
set DO_HAD_LEV_JETS    [ getenv DO_HAD_LEV_JETS   0               ]
set DO_PAR_LEV_JETS    [ getenv DO_PAR_LEV_JETS   0               ]
set HERWIG             [ getenv HERWIG            0               ]

source $WORK_DIR/JetCluModule.tcl

useRCP set f
module disable  PuffModule
module enable  SimInitManager 


path disable AllPath
path create  StntuplePath    ManagerSequence 

module disable  JetCluModule-cone0.4   JetCluModule-cone0.7   JetCluModule-cone1.0   

if { $DO_PAR_LEV_JETS != 0 } {
  module enable  JetCluModule-par-cone0.4 JetCluModule-par-cone0.7 \
                 JetCluModule-par-cone1.0
  path   append StntuplePath JetCluModule-par-cone0.4 \
                             JetCluModule-par-cone0.7 \
                             JetCluModule-par-cone1.0
}

if { $DO_HAD_LEV_JETS != 0 } {
  module enable  JetCluModule-had-cone0.4 JetCluModule-had-cone0.7 \
                 JetCluModule-had-cone1.0 
  path   append StntuplePath JetCluModule-had-cone0.4  \
                             JetCluModule-had-cone0.7  \
                             JetCluModule-had-cone1.0
}

module enable  InitStntuple StntupleMaker FillStntuple
path append StntuplePath     InitStntuple              \
			     StntupleMaker             \
			     FillStntuple

path enable StntuplePath

module talk StntupleMaker
processName   set $PROCESS_NAME

if { $DO_HAD_LEV_JETS != 0 || $DO_PAR_LEV_JETS != 0 } {
# this will be "JetBlock"
  jetCollName set PROD@JetCluModule-par-cone0.4
}

if { $DO_PAR_LEV_JETS != 0 } {
  jetCollName add PROD@JetCluModule-par-cone0.4 \
                  PROD@JetCluModule-par-cone0.7 \
                  PROD@JetCluModule-par-cone1.0
}

if { $DO_HAD_LEV_JETS != 0 } {
  jetCollName add PROD@JetCluModule-had-cone0.4 \
                  PROD@JetCluModule-had-cone0.7 \
                  PROD@JetCluModule-had-cone1.0 
}

#-----------------------------------------------------------------------
# define output STNTUPLE file
#-----------------------------------------------------------------------
  histfile $OUTPUT_STNTUPLE
#-----------------------------------------------------------------------
#        turn OFF splitting (ROOT 3.01/06)
#-----------------------------------------------------------------------
  splitMode set 1
#
  makeCalData          set 0
  makeCcrData          set 0
  makeCesData          set 0
  makePesCorrectedData set 0
  makeClcData          set 0
  makeCprData          set 0
  makeCp2Data          set 0
  makeEmtData          set 0
  makeHatData          set 0
  makeEmTiming         set 0
  makeElectrons        set 0
  makePhoenixElectrons set 0
  makeDcasData         set 0
if { $DO_HAD_LEV_JETS != 0 || $DO_PAR_LEV_JETS != 0 } {
  makeJets             set 1
} else {
  makeJets             set 0
}
  makeMet              set 0
  makeMuons            set 0
  makeCosmic           set 0
  makePhotons          set 0
  makeClusters         set 0
  makeTaus             set 0
  makeTofMatches       set 0
  makeTracks           set 0
  makeTrigger          set 0
  makeTrigSim          set 0
  makeVertices         set 0
  makeZVertices        set 0
  makeSecVtxTag        set 0
  makeJetProb          set 0
  makePi0s             set 0
  makeFwdDetData       set 0
  makeConversions      set 0
  makeL3Summary        set 0
  makeGenp          set 1
  makeObsp          set 0
  makeCotData       set 0
  makePesData       set 0
  makeTopTElectrons set 0
  makeTopLElectrons set 0
  makeTopTJets      set 0
  makeTopLJets      set 0
  makeTopMet        set 0
  makeL3Muons       set 0
  makeL3Taus        set 0
  makeTopTMuons     set 0
  makeTopLMuons     set 0
  makeTopSummary    set 0
  makeSvxData       set 0
  makeSiStrips      set 0
  makeTrackLinks    set 0
  makeTofData       set 0
  makeTopSummary    set 0
  makeXft           set 0
  storeXftHits      set 0
  storeXftPixels    set 0
  storeXftTrack     set 0
  makeSvt           set 0
  makeSiGeantIsect  set 0
  makeSiIsect       set 0
  makeCmuData       set 0
  makeCmpData       set 0
  makeCmxData       set 0
  makeBmuData       set 0
  makeBsuData       set 0
  makeTsuData       set 0
exit

module disable ConfigManager



