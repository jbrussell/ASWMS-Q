#!/bin/csh
# script to run plot_wk and choose a specific 
# to be output
# B1 represents the first mode you want
# B2 represents the last mode you want

set TABLE='pa5.t0to50'
set B1=$1
set B2=$2

plot_wk <<!
table $TABLE.table_hdr
search
1 0. 51.
99 0 0
branch $TABLE.table_hdr.branch
search
1 0. 51.
99 0 0
outbran
$TABLE.$1_$2
$B1
$B2
Branches $B1 to $B2
quit
!
echo "Done!"
