#!/bin/sh
if [ $# -lt 1 ]
then 
	echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
	echo ":: This is to get the status (how many events are processed) of a job ::"
	echo ":: from a log file. Usefull when running in ATOM cluster. I needed    ::"
	echo ":: this to be fast and efficient by ignoring the billions of Warnings ::"
	echo ":: spitted out from the JetFilterModule.                              ::"
	echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
	echo ":: Require file/s name/s. exiting!                                    ::"
	echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
	exit -1
else
		tail -s 2 -F $* | grep -v Warn | grep -v Hello 
fi
