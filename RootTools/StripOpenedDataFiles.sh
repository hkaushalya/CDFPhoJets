#! /bin/sh

#this will scan my job out files and make a
#list of files that were opened for processing
if [ $# -lt 1 ]
then
	echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
	echo ":: Wrote this to scan my CAF job out files and strip out the data     ::"
	echo ":: files processed by the job. This will create a file with those file::"
	echo ":: name, which can be fed into the FindProcFiles.sh to be compared to ::"
	echo ":: the reference list of files.                                       ::"
	echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
	echo ":: Pl. specify file/s to scan.                                        ::"
	echo ":: Argument 1. list of files to scan                                  ::"
	echo ":: Argument 2. output file name                                       ::"
	echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
	exit 1
fi

outfile=filesopened.txt
#$if [ $# -gt 1 ]
#$then
#	outfile=$2 # this is wrong. if I do this, it will take the second argument in the $*
# when someone gives a list of files using a wild card
#fi

for FILE in  $*
do
	grep Open $FILE | cut -d / -f13 >> $outfile
done
echo output written to $outfile
