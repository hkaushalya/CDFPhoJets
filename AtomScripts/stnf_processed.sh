#!/bin/sh
#created this to monitor the progress of stntuple making jobs
#running on ATOM cluster

if [ $# -lt 1 ]
then 
	echo "require file name. exiting!"
	exit -1
else
	for file in $*
	do
		#since I am using SAM dataset there is noway to determine the total number of files
		#that will be processed by a signle job(i have set SAM FILE LIMIT to ZERO!)
		#total=`cat $file | grep 'phys\"' | grep -v Open | grep -v Clos |wc -l`
		done=`cat $file | grep 'Closing Input File' |wc -l`
		#echo "$file => Files done/ File requested = $done / $total" 
		echo "$file => Files done = $done" 
	done
fi
