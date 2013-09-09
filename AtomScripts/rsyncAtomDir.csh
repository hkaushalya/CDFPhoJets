#! /bin/tcsh
if ($# < 3) then
	echo "=====================================================================|"
	echo "| argument 1 = absolute path of the source dir/file                  |"
	echo "| argument 2 = destination dir                                       |"                               
	echo "| argument 3 = file type  (eg stn.1429)                                           |"                               
	echo "| require absolute path of the ATOM dir and destintaion dir. exiting!|"
	echo "=====================================================================|"
	exit -1
endif
set nodes="Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium"
#set nodes="Gadolinium"
set SOURCE_DIR="$1"
set DEST_DIR="$2"
set FILE_TYPE="$3"
echo "SOURCE_DIR = $SOURCE_DIR"
echo "DEST_DIR   = $DEST_DIR"
echo "FILE_TYPE  = $FILE_TYPE"

foreach node ( $nodes )
	echo "//${node}"
	set filelist=`ssh $node "ls -1 $SOURCE_DIR/ | grep $FILE_TYPE"`
	foreach file ( $filelist )
		echo $file
		rsync -av ${node}:${SOURCE_DIR}${file}  ${DEST_DIR}
	end
end
