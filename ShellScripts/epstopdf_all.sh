#!/bin/sh
echo "Converting all eps file in ${PWD} to pdf.. "
	#for file in `ls -1 | grep eps$`
	for file in `ls -1tr | grep eps$`  #preserve time order
	do
		echo "converting $file to pdf"
	 	`epstopdf $file`
	done 
echo "Conversion done." 
