#!/usr/bin/bash

# Path to the data
# Note that the existence of a two sub-repositories named respectively
# XMM-Newton and Swift-UVOT with the tables of single positions per 
# catalogue in it named respectively SUSS6_ra_dec_per_src.fits 
# and uvot_ra_dec_per_src.fits is necessary. Also a subrepository where 
# the output files will be saved named Match_OM within Swift-UVOT is required.
datapath="../../../Data"

echo "Performing match calculation UVOT-OM catalogue skyerr"
stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/XMM-Newton/SUSS6_ra_dec_per_src.fits \
in2=$datapath/Swift-UVOT/uvot_ra_dec_per_src.fits matcher=skyerr \
values1="RA DEC 3*POSERR" values2="RA DEC 3*POSERR" \
ocmd="addcol angDist skyDistanceDegrees(RA_OM,DEC_OM,RA_UVOT,DEC_UVOT)*3600" \
join=1and2 find=best ofmt=fits params=1 suffix1='_OM' suffix2='_UVOT' \
out=$datapath/Swift-UVOT/Match_OM/match_uvot_om progress=none

# The random matches are evaluated moving sucessively the right ascensions
# and declinations by 1' and 2'. The result of the counting is saved in the
# the file count_random_match.txt.
echo "Computing random match values"
echo "Random match count" > $datapath/Swift-UVOT/Match_OM/count_random_match.txt
shift_list='+0.016667 -0.016667 +0.03333 -0.03333'
for shift in $shift_list; do 
    echo "Shift: "$shift
    
    # Shift of the right ascension
    stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/XMM-Newton/SUSS6_ra_dec_per_src.fits \
    in2=$datapath/Swift-UVOT/uvot_ra_dec_per_src.fits matcher=skyerr \
    values1="RA"$shift" DEC 3*POSERR" values2="RA DEC 3*POSERR" \
    join=1and2 find=best params=1 omode=count progress=none >> $datapath/Swift-UVOT/Match_OM/count_random_match.txt

    # Shift of the declination
    stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/XMM-Newton/SUSS6_ra_dec_per_src.fits \
    in2=$datapath/Swift-UVOT/uvot_ra_dec_per_src.fits matcher=skyerr \
    values1="RA DEC"$shift" 3*POSERR" values2="RA DEC 3*POSERR" \
    join=1and2 find=best params=1 omode=count progress=none >> $datapath/Swift-UVOT/Match_OM/count_random_match.txt
done






