#!/bin/sh
echo "Converting all eps file in ${PWD} to pdf.. "
	#for file in `ls -1 | grep eps$`
	for file in `ls -1tr | grep eps$`  #preserve time order
	do
		echo "converting $file to png"
		`convert -density 100 $file -flatten ${file%.*}.png`
	done 
echo "Conversion done." 
