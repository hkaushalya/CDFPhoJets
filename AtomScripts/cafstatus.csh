#! /bin/tcsh
##########################################
# Created this to monitor the CAF jobs
# given farm, JID, and sections this will
# indicate the events to be processed and 
# events processed.
##########################################

if ($# =~ 0) then
	echo "cafstatus FARM JID N-SECTIONS"
	exit 1
endif

if ($# !~ 3) then
	echo "Require Farm, JID, and Number of job sections."
	exit 1
endif
	
set CAF=$1
set JID=$2
set SECTIONS=$3
set loop=1
set loopto=0
#set tempfile="${CAF}_${JID}.caftmp"
#if (-f $tempfile) then
#	set name=`date`
	#$tempfile = `echo ${name}_${tempfile}.tmp`
#	$tempfile = `echo ${name}_${CAF}_${JID}.${tempfile}`
#endif
	#echo $tempfile
@ loopto = $SECTIONS + $loop
echo "Looking in $CAF for JID=$JID for $SECTIONS job sections."
source ~cdfsoft/cdf2.cshrc
setup cdfsoft2 development
#echo "JID-Sec :: Processed / Requested"
#while ( $loop !~ $loopto )
#		echo -n "${JID}-${loop} :: "
	 	#`source ~cdfsoft/cdf2.cshrc && setup cdfsoft2 6.1.4 && CafMon --farm $CAF cat $2 $loop job_${loop}.out`
#	 	CafMon --farm $CAF cat $2 $loop job_${loop}.out >&! $tempfile
#		set processed=`grep "Processed" $tempfile | tail -1 | awk '{print $8}'`
#		set need=`grep -m1 "TStnRun2InputModule::BeginJob: chained" $tempfile | awk '{print $5}'`
#		echo "${JID}-${loop} :: $processed / $need" & 
#		@ loop = $loop + 1
#end
#/bin/rm $tempfile

set section=1

while ( $section !~ $SECTIONS )
		set tempfile = `echo ${JID}_${section}.out`
		echo "processing section $section : output file = $tempfile"
	 	#CafMon --farm $CAF cat $2 $section job_${section}.out >&! $tempfile
	 	CafMon --farm $CAF cat $2 $section job_${section}.out | grep -i "Opened" $tempfile
		@ section = $section + 1
end
		
		
