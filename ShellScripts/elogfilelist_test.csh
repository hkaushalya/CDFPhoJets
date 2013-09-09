#! /bin/tcsh

set newdir=$1
set filetype=$2
set files=`ls -1 *.${filetype}`
foreach file ( $files)
  if ( -f $file ) then
    set id = `printf "-f %s" $file`
	 echo "${id}"
  endif
end
