#! /bin/sh
##########################################
# Created this to monitor the CAF jobs
# given farm, JID, and sections this will
# indicate the events to be processed and 
# events processed.
##########################################

if [ $# -eq 0 ] then
	echo "cafstatus FARM JID N-SECTIONS"
	exit 1
fi

if [ $# -ne 3 ] then
	echo "Require Farm, JID, and Number of job sections."
	exit 1
fi
	
CAF=$1
JID=$2
SECTIONS=$3
loop=1
loopto=0
tempfile="${CAF}_${JID}.caftmp"
if (-f $tempfile) then
	set name=`date`
	#$tempfile = `echo ${name}_${tempfile}.tmp`
	$tempfile = `echo ${name}_${CAF}_${JID}.${tempfile}`
fi
	#echo $tempfile
loopto=`expr $SECTIONS + $loop`
echo "Looking in $CAF for JID=$JID for $SECTIONS job sections."
source ~cdfsoft/cdf2.cshrc
setup cdfsoft2 development
echo "JID-Sec :: Processed / Requested"
while [ $loop -ne $loopto ]
	do
		echo -n "${JID}-${loop} :: "
	 	#`source ~cdfsoft/cdf2.cshrc && setup cdfsoft2 6.1.4 && CafMon --farm $CAF cat $2 $loop job_${loop}.out`
	 	CafMon --farm $CAF cat $2 $loop job_${loop}.out >&! $tempfile
		processed=`grep "Processed" $tempfile | tail -1 | awk '{print $8}'`
		need=`grep -m1 "TStnRun2InputModule::BeginJob: chained" $tempfile | awk '{print $5}'`
		echo "${JID}-${loop} :: $processed / $need" & 
		loop=`expr $loop + 1`
done
/bin/rm $tempfile
