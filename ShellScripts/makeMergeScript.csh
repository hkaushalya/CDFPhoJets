#! /bin/tcsh
if ($# < 1) then
	echo "require absolute path of the dir. exiting!"
	exit -1
endif

echo "{"
echo 'TChain ch("Stuple");'
	set filelist=`ls -1 $1 | grep root`
	foreach file ( $filelist )
		echo 'ch.Add("'${file}'");'
	end
echo 'std::cout << "Chained " << ch.GetListOfFiles()->GetEntries() << " files with " << ch.GetNtrees() << " trees and entries " << ch.GetEntries() << std::endl;'
echo 'ch.Merge("StupleV8_MCPythia_WW_NoJetVtxCuts.root");'
#echo 'ch.Merge("MergedStuple.root");'
echo "}"

