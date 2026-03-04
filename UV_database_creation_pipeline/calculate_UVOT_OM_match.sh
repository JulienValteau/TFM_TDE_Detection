#!/usr/bin/bash

datapath="../../../Data"

echo "Performing match calculation UVOT-OM catalogue skyerr"
stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/XMM-Newton/SUSS6_ra_dec_per_src.fits \
in2=$datapath/Swift-UVOT/uvot_ra_dec_per_src.fits matcher=skyerr \
values1="RA DEC 3*POSERR" values2="RA DEC 3*POSERR" \
ocmd="addcol angDist skyDistanceDegrees(RA_1,DEC_1,RA_2,DEC_2)*3600" \
join=1and2 find=best ofmt=fits params=1 out=$datapath/Swift-UVOT/Match_OM/match_uvot_om

stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/XMM-Newton/SUSS6_ra_dec_per_src.fits \
in2=$datapath/Swift-UVOT/uvot_ra_dec_per_src.fits matcher=skyerr \
values1="RA+0.016667 DEC 3*POSERR" values2="RA DEC 3*POSERR" \
join=1and2 find=best params=1 omode=count > $datapath/Swift-UVOT/Match_OM/count_random_match.txt




