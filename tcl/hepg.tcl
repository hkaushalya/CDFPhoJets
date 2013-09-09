# this is to make HEPG only stntuple format job output.
#this sources the stnmaker-hepg.tcl
#-----------------------------------------------------------------------
#  example TCL file to process 100 events
#-----------------------------------------------------------------------
set     WORK_DIR          $env(PWD)
source  $WORK_DIR/TclUtils/scripts/getenv.tcl
#
creator set               PROD
set     PROCESS_NAME      PROD
set     OUTPUT_STNTUPLE   ./stn.$env(CAF_JID)_$env(CAF_SECTION).root
puts "OUTPUT_STNTUPLE = $OUTPUT_STNTUPLE"
set     DO_EMBEDDING      ""
set     L3_SOURCE         ""; # TL3D
#set     env(ALL_JETS)     ""; # 1= do additional jet collections and tags
#
module input DHInput
module talk  DHInput
#
#  change the line below to include the input file you need
#
include file ./*.root
#
show
#
exit

source $WORK_DIR/stnmaker-hepg.tcl


ev begin
#ev begin -nev 100
show
exit

exit

