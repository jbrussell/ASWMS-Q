#!/bin/csh
# script to run plot_wk and choose a specific 
# to be output
# BR represents the mode you want to output

set TABLE='pa5.t0to50'
set BR=$1

plot_wk <<!
table $TABLE.table_hdr
search
1 0. 51.
99 0 0
branch $TABLE.table_hdr.branch
search
1 0. 51.
99 0 0
output $TABLE.$BR b
$BR

quit
!
echo "Done!"
