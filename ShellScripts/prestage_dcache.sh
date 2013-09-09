#! /bin/bash

. ~cdfsoft/cdf2.shrc

prestagedelay=10
nprestage=0
maxprestage=10
blockdelay=20

process_file () {
    echo Prestaging file `date`
    dccp -P $1
    rv=$?
    let nprestage+=1
    if [ $nprestage -ge $maxprestage ];then
	nprestage=0
	# Check if the last file is staged or not
	if dc_check $1 2> /dev/null | grep -q FAILED ; then
	    echo Sleeping for $blockdelay seconds
	    sleep $blockdelay
	fi
    fi
    sleep $prestagedelay
    return $rv
}

run_sam_consumer () {

# Make sure that the project name is defined
: ${SAM_PROJECT:?}

local cid cpid filename status

echo Establishing consumer...
cid=$(sam establish consumer -s --appName=demo --appVersion=1)|| return 1
cid=${cid/CID /}
echo Consumer established with CID $cid

echo Establishing consumer process...
cpid=$(sam establish consumer process --cid=$cid ) || return 1
cpid=${cpid/CPID /}
echo Consumer process established with CPID $cpid

echo Requesting first file...
filename=$(sam get next file --cpid=$cpid --timeout=600) || return 1
filename=${filename/File /}

while [ "$filename" != "END OF STREAM" ]; do
    echo Got file $filename
    status=bad
    process_file $filename
    if [ $? -eq 0 ]; then
	status=ok
    fi
    sam release file --cpid=$cpid --status=$status --file=$filename
    echo Requesting next file...
    filename=$(sam get next file --cpid=$cpid --timeout=600) || return 1
    filename=${filename/File /}
done
echo Received end of file stream

}

# first arg is dataset definition
run_sam_project () {

local defname extraprojectargs

echo Starting SAM project ${SAM_PROJECT:?}...

if [ "$SAM_FILECUT" ]; then
    extraprojectargs="--filecut=$SAM_FILECUT"
fi

sam start project --group=test --defname=${SAM_DEFNAME:?} $extraprojectargs || exit 1

echo Project started

run_sam_consumer
local rval=$?

sam dump project
sam stop project --force 

return $rval
}

setup sam
setup dcap

SAM_DEFNAME=$1
SAM_FILECUT=$2

# Make project name if none defined
if [ -z "$SAM_PROJECT"]; then
    SAM_PROJECT=`whoami`_`date +%H%M%S_%m%d%Y`
fi
export SAM_PROJECT

run_sam_project
