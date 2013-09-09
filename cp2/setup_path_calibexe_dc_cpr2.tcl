# calibexe for CPR2 calibrations
# outputs an Stntuple
# vars to be set: CPR2_STNFILE

# Define the CPR2 Stntuple path and prereq:
module enable ManagerSequence
#module enable CP2QModule
module enable InitStntuple
module enable StntupleMaker
module enable FillStntuple

path create CP2path  \
                        ManagerSequence         \
                        InitStntuple            \
                        StntupleMaker           \
                        FillStntuple
			

path enable CP2path


module talk StntupleMaker

  histfile $CPR2_STNFILE
  processName   set PROD 
  splitMode            set 1
#        
  makeCalData          set 0
  makeCcrData          set 0
  makeCesData          set 0
  makePesCorrectedData set 0
  makeClcData          set 0
  makeCprData          set 0
  makeEmtData          set 0
  makeHatData          set 0
  makeElectrons        set 0
  makeDcasData         set 0
  makeJets             set 0
  makeMet              set 0
  makeMuons            set 0
  makeCosmic           set 0
  makePhotons          set 0
  makeClusters         set 0
  makeTaus             set 0
  makeTracks           set 0
  makeTrigger          set 0
  makeTrigSim          set 0
  makeVertices         set 0
  makeZVertices        set 0
  makeSecVtxTag        set 0
  makeJetProb          set 0
  makePi0s             set 0
  makeFwdDetData       set 0
  makeL3Summary        set 1
  makeSvt              set 0
  makePhoenixElectrons set 0
  makeCp2Data 	       set 1 
  show
exit


