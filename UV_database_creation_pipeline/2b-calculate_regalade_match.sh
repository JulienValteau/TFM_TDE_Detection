#!/usr/bin/bash

# Path to the data
# Note that the existence of a two sub-repositories named respectively
# XMM-Newton and Swift-UVOT with the tables of single positions per 
# catalogue in it named respectively SUSS6_ra_dec_per_src.fits 
# and uvot_ra_dec_per_src.fits is necessary and a subrepository Regalade
# with the extracted columns from the complete Regalade catalog using TOPCAT
# named regalade_ra_dec_plus.fits is necessary. Also a subrepository where 
# the output files will be saved named Match_OM within Swift-UVOT is required
datapath="../../../Data"
echo "Performing match calculation OM-Regalade catalogue skyerr"
stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/XMM-Newton/SUSS6_ra_dec_per_src.fits \
in2=$datapath/Regalade/regalade_ra_dec_plus.fits matcher=skyerr \
values1="RA DEC 3*POSERR" values2="ra dec 0" \
ocmd="addcol angDist skyDistanceDegrees(RA_1,DEC_1,ra_2,dec_2)*3600" \
join=1and2 find=best ofmt=fits params=1 out=$datapath/XMM-Newton/Match_regalade/match_om_regalade

echo "Computing random match values"
stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/XMM-Newton/SUSS6_ra_dec_per_src.fits \
in2=$datapath/Regalade/regalade_ra_dec.fits matcher=skyerr \
values1="RA+0.016667 DEC 3*POSERR" values2="ra dec 0" \
join=1and2 find=best params=1 omode=count > $datapath/XMM-Newton/Match_regalade/count_random_match.txt

datapath="../../../Data"
echo "Performing match calculation UVOT-Regalade catalogue skyerr"
stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/Swift-UVOT/uvot_ra_dec_per_src.fits \
in2=$datapath/Regalade/regalade_ra_dec_plus.fits matcher=skyerr \
values1="RA DEC 3*POSERR" values2="ra dec 0" \
ocmd="addcol angDist skyDistanceDegrees(RA_1,DEC_1,ra_2,dec_2)*3600" \
join=1and2 find=best ofmt=fits params=1 out=$datapath/Swift-UVOT/Match_regalade/match_uvot_regalade

echo "Computing random match values"
stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/Swift-UVOT/uvot_ra_dec_per_src.fits \
in2=$datapath/Regalade/regalade_ra_dec.fits matcher=skyerr \
values1="RA+0.016667 DEC 3*POSERR" values2="ra dec 0" \
join=1and2 find=best params=1 omode=count > $datapath/Swift-UVOT/Match_regalade/count_random_match.txt
