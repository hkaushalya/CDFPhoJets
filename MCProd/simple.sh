#! /bin/sh

JOBNO=$1
DATASET_ID=$2
MYOUTPUT_LOCATION=$3

if [ -d $3 ]
then
	echo "Output location $MYOUTPUT_LOCATION found"
else
	echo "Output location $MYOUTPUT_LOCATION is not found! exiting!"
fi	

BOOK=cdfpexo
JOB_OUTPUT_DIR=. #current dir
USER_GENMODE=$4
FARM=cdfgrid

export ISR=$5
export FSR=$6
export Q2=$7
export minPtHat=$8
export minPhoPt=$9
echo "========= simple.sh : GOT VALUES"
echo "JOBNO          = $JOBNO"
echo "USER_GENMODE   = $USER_GENMODE"
echo "BOOK           = $BOOK"
echo "DATASET_ID     = $DATASET_ID"
echo "JOB_OUTPUT_DIR = $JOB_OUTPUT_DIR"
echo "FARM           = $FARM"
echo "ISR            = $ISR"
echo "FSR            = $FSR"
echo "Q2             = $Q2"
echo "minPtHat       = $minPtHat"
echo "FilterminPhoEt = $minPhoPt"
echo "================= env before cdfsoft setup==========="
env
echo "============================"


#MCprod seem to set this when run in CAF
#source ~cdfsoft/cdf2.shrc
#setup cdfsoft2 development
./mcProduction/scripts/run1segment -J $1 -b $BOOK -d $DATASET_ID -o $JOB_OUTPUT_DIR -f $FARM -x USER_GENMODE=$USER_GENMODE
echo "===========ntupling ========"
source ~cdfsoft/cdf2.shrc
setup cdfsoft2 6.1.4
export DO_HAD_LEV_JETS=1
export DO_PAR_LEV_JETS=1
env | grep "LEV_JETS"
./stnmaker_prod.exe hepg.tcl

echo "=================files after stn making ======================= " 
ls -ltrh
echo "=================env before rcp =============================== "
env
echo "=============================================================== "

#for file in `ls | grep stn | grep $CAF_JID`
#do
#	echo "copying $file using rcp"
	#fcp -c $KRB5BIN_DIR/rcp -N -r $file $MYOUTPUT_LOCATION 
#	/usr/krb5/bin/rcp -N -r $file $MYOUTPUT_LOCATION 
#done
#echo "=================env after  rcp =============================== "
#env
#echo "=============================================================== "


rm -rf stnmaker_prod.exe
rm -rf CalTrigger GNUmakefile Production RootObjs SimulationMods TclUtils
rm -rf TriggerMods XFTSim XTRPSim bin caf_remove.lst cafsubmitMCProd cdfopr cesData
rm -rf commands dab dbt etc farmsonly fcp include mcProduction oracle rcpdb setupmcprodtestjob shlib svtsim
rm -rf foo* bar tcl scripts cdfpexo book MCTests cint *.dat
rm -rf *concatenate *concat_tcl *.output
rm -rf *.tmp
rm -rf gen*
#rm -rf *.root*
#rm -rf pexo* # not good. removes the log files too!
echo "======== end simple.sh"
