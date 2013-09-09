#! /bin/tcsh
#this will grep for a pattern in specific file type
#this is pretty simple. I'll make is more fancy later
#Sep -7, 2009

#echo [`date`] "Starting grep"
echo "File type : *." $1
echo "Pattern   :  " $2
echo "Using default settings: ignore case/file names/ line numbers/colors"

#opt=
#line number = -n
#file name = -H
#files with match only = -l
#ignore case -i
#only matching part of a  line = -o
#find . -type f -name '*.cc' -exec grep -l -n --color "GetElectronIdWord" {} \;

#this is the default case
find . -type f -name "*.$1" -exec grep -H -n -i --color "$2" {} \;
