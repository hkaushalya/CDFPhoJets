#!/bin/sh
#This need the output of the 'getzeros'
#search for CP2 HV off periods and spits out the run::section-setion
#where the HV was off.
#

if [ $# -lt 1 ]
then 
	echo "require input file. exiting!"
	exit -1
else
	debug=0
	if [ $# -eq 2 ]
	then
		debug=$2
	fi
		
	echo "File >> $1"
	runs=`cat $1 | cut -d: -f1| sort | uniq`
	#runs=`cat $1 | grep -v [A-Za-z] | cut -f1| sort | uniq`
	echo -n "Runs HV off ::"
	echo $runs

	for run in $runs
	do
		if [ $debug -eq 1 ]
		then
			echo "Searching for run = $run" 
		fi
		
		sec=`cat $1 | grep -e "^${run}" | cut -d: -f3 | sort | uniq`
		#sec=`cat $1 | grep $run | cut -f2 | sort | uniq`
		narr=${#sec}  # this or many other forms I tried did not work! 
		#echo "${sec}" 
		minsec=99999
		maxsec=-1

		for section in $sec
		do
			if [ $debug -eq 1 ]
			then
				echo "$section"
			fi

			if [ $section -gt $maxsec ]
			then
				maxsec=${section}
			fi

			if [ $section -lt $minsec ]
			then
				minsec=${section}
			fi
		done
		
		if [ $debug -eq 0 ]
		then
			#echo "for run $run::${minsec}-${maxsec}"
			echo "$run::${minsec}-${maxsec}"
		else
			#the info in brackets can indicate if I am fooled by intermittent stuck bits.
			#in anycase better to look at dump to make sure if there is a period
			# or just bit stucked sections.
			#echo "for run $run::$minsec - $maxsec (had $narr elements in the searh array)"
			echo "for run $run::${minsec}-${maxsec}"
		fi

	done 
fi
