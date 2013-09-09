#! /bin/tcsh

#set nodes="Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium nbay01 nbay02 nbay03 nbay04 nbay05"
#set nodes="nbay01 nbay02 nbay03 nbay04 nbay05"
set nodes="nbay04 nbay05"

foreach node ( $nodes )
	echo ">>>>>>>>>>>>>>> $node <<<<<<<<<<<<<<<<<<<<<<"
	
		#set list1=`ssh $node "ps aux | grep "root.exe" | grep samantha "`
		#set list1=`ssh $node "ps aux"`
		#echo ${list1}
		#set list=`ssh $node "ps aux | grep "root.exe" | grep -v "grep" | cut -d ' '  -f2"`
		#ssh $node "ps ax -u samantha | grep 'root.exe' | cut -c 1-5"
		#set list=`ssh $node "ps ax -u samantha | grep "root.exe" | grep -v "grep" | awk '{print ${1}}'"`
		set list=`ssh $node "ps ax -u samantha | grep "root.exe" | grep -v "grep" | cut -c 1-5"`
	 	foreach ps ( $list ) 
   		echo "renice job ${ps}"
   		ssh $node "renice 19 ${ps}"
	 	end
end
