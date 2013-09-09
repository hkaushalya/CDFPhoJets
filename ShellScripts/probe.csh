#! /bin/tcsh

set nodes="Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium"


set newdir="/cdf/scratch/samantha/stuple"
set a=1
set sleep=300
while ( $a == 1 )
	clear
foreach node ( $nodes )
	echo ">>>>>>>>>>>>>>> $node <<<<<<<<<<<<<<<<<<<<<<"
	
	#if [ $1 -eq 1 ]
	#then 
    	#ssh $node "ps -al r -U $USER | grep root"
   	#ssh $node "ps -al -U $USER | grep root"

	#elif [ $1 -eq 2 ]
	#then
	#	set list=`ssh $node "ls /cdf/scratch/samantha/zee/dataset*"`
	# 	foreach file ( $list ) 
	# 		echo $file
   #		ssh $node "~/summary $file"
	# 	end

#	else 
		
		#set list=`ssh $node "ls /cdf/scratch/samantha/phodata/dataset*"`
		#set list=`ssh $node "ls /cdf/scratch/samantha/cosmicskim/dataset*"`
		#set list=`ssh $node "ls ${newdir}/dataset*"`
		#set list=`ssh $node "ls ${newdir}/dataset_10_02*"`
		#set list=`ssh $node "ls ${newdir}/dataset_32_03*"`
		set list=`ssh $node "ls ${newdir}/dataset_10_04*"`
		#set list=`ssh $node "ls /cdf/scratch/samantha/zee/dataset*"`
	 	foreach file ( $list ) 
	 		#echo $file
   		ssh $node "~/processed $file"
	 	end
		set list=`ssh $node "ls ${newdir}/dataset_10_05*"`
		#set list=`ssh $node "ls /cdf/scratch/samantha/zee/dataset*"`
	 	foreach file ( $list ) 
	 		#echo $file
   		ssh $node "~/processed $file"
	 	end


	#fi
end
	echo "sleeping for $sleep seconds"
	sleep $sleep
end
