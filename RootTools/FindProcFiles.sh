#! /bin/sh

#this will scan my job out files and make a
#list of data files that were processed in the job
#I think I made this when I had the duplicate file issue.

if [ $# -lt 2 ]
then
	echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
	echo ":: Wrote this to find the duplcate files processed in my CAF jobs.    ::"
	echo ":: Need two files. One with the list of reference file names. Other   ::"
	echo ":: with the list file names need to compared to the reference. This   ::"
	echo ":: will create a file 'notfound.files' with a list of file names that ::"
	echo ":: are in the reference list but not found in the one compared.       ::"
	echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
	echo ":: Pl. specify file/s to scan.                                        ::"
	echo ":: parameter 1. name of the reference file                            ::"
	echo ":: parameter 2. name of the file need to be compared                  ::"
	echo ":: eg. - FindProcFiles.sh ref_files.txt new_files.txt                 ::"
	echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
	exit 1
fi
refwords=`cat $2`
#echo $refwords
words=`cat $1`
#echo $words


outfile=notfound.files

totalfiles=0
foundfiles=0
notfoundfiles=0

for word1 in $refwords
do
	totalfiles=`expr $totalfiles + 1`
	found=0
	for word2 in $words
	do
		if [ $word1 = $word2 ]
		then
			found=1
			foundfiles=`expr $foundfiles + 1`
			#echo found $word1
			break
		fi
	done

	if [ $found -eq 0 ]
	then
		notfoundfiles=`expr $notfoundfiles + 1`
		echo $word1 >> $outfile
	fi
done
echo output written to $outfile
echo Files not found $notfoundfiles/$totalfiles
