#! /bin/tcsh
set nodes="Lanthanum Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium Actinium"
set srcdir="/cdf/scratch/samantha/cosmicskim"
#set srcdir="/cdf/scratch/samantha/zee"
set dest="lanthanum"
set destdir="/cdf/scratch/samantha/cosmicmerge"
#set destdir="/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicJetEmMismatchedStntuples/"
set file="\*\.log"
set i=0
  foreach node ( $nodes )
  		echo "${node} :: ${i}"
		#set p=32
		#set q=1
    	#set id = `printf "%02i_%02i_%02i.log" $p $q $i`
  	 	#scp  $node\:${srcdir}/${id} $destdir 
		
  	 	scp  $node\:${srcdir}/Cosmics_phodata\*.root ${dest}\:${destdir}
		#mv $destdir/CosmicsAndJetEmMismatch.root $destdir/CosmicsAndJetEmMismatch_${i}.root
		
#  	scp  $node\:${srcdir}/Stuple_wenmc_2_$i.root $dest\:$destdir &
#		scp  $node\:${srcdir}/Stuple_wenmc_3_$i.root $dest\:$destdir &
#		scp  $node\:${srcdir}/Stuple_wenmc_4_$i.root $dest\:$destdir &
#		scp  $node\:${srcdir}/Stuple_wmnmc_1_$i.root $dest\:$destdir 
#		scp  $node\:${srcdir}/Stuple_wmnmc_2_$i.root $dest\:$destdir &
#		scp  $node\:${srcdir}/Stuple_wmnmc_3_$i.root $dest\:$destdir &
#		scp  $node\:${srcdir}/Stuple_wmnmc_4_$i.root $dest\:$destdir 
#		scp  $node\:${srcdir}/Stuple_wtnmc_1_$i.root $dest\:$destdir &
#		scp  $node\:${srcdir}/Stuple_zeemc_1_$i.root $dest\:$destdir &
#		scp  $node\:${srcdir}/Stuple_zeemc_2_$i.root $dest\:$destdir 
#		scp  $node\:${srcdir}/Stuple_zeemc_3_$i.root $dest\:$destdir &
#		scp  $node\:${srcdir}/Stuple_zmmmc_1_$i.root $dest\:$destdir &
#		scp  $node\:${srcdir}/Stuple_zmmmc_2_$i.root $dest\:$destdir 
#		scp  $node\:${srcdir}/Stuple_zmmmc_3_$i.root $dest\:$destdir &
#		scp  $node\:${srcdir}/Stuple_zmmmc_4_$i.root $dest\:$destdir &
#		scp  $node\:${srcdir}/Stuple_zttmc_1_$i.root $dest\:$destdir 
#		scp  $node\:${srcdir}/Stuple_zttmc_2_$i.root $dest\:$destdir 
		@ i = $i + 1
	end


