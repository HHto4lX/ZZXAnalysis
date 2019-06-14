#!/bin/tcsh

set nonomatch
set list = ( */*.root */*.corrupted */*.recovered */*.gz */*.txt */core* */jobid  */LSFJOB*/ */error/*.err */log/*.log */output/*.out)

foreach f ( ${list} )
    if ( -e $f ) then
	rm -r $f
    endif
end

