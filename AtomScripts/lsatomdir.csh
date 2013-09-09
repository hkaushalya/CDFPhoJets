#! /bin/tcsh
if ($# < 1) then
	echo "require absolute path of the dir. exiting!"
	exit -1
endif
#set nodes="Lanthanum Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium Actinium"
set nodes="Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium"
echo "Listing dir $1"

echo $1
foreach node ( $nodes )
	# I am writing things in a format so I can copy and paste
	# this into the cp2 gain check driver script.
	echo "//${node}"
	set filelist=`ssh $node "ls -1 $1/*.stn*"`
	foreach file ( $filelist )
	#echo $file
		echo 'inChain.Add("root://'${node}'.fnal.gov:5151/'${file}'");'
	end
end
