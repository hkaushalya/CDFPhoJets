#! /bin/tcsh
#copy the file startrootd.csh to ATOM node home
#and run this script to start a rootd in each node

set nodes="Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium"

foreach node ( $nodes )
	echo ":::::: $node"
	ssh $node "~/krbhelper.csh ~/startrootd.csh"
end
