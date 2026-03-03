#!/usr/bin/bash

datapath="../../../Data"
echo "Performing match calculation OM-Regalade catalogue skyerr"
stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/XMM-Newton/SUSS6_ra_dec_per_src.fits \
in2=$datapath/Regalade/regalade_ra_dec.fits matcher=skyerr \
values1="RA DEC 3*POSERR" values2="ra dec 0" \
ocmd="addcol angDist skyDistanceDegrees(RA,DEC,ra,dec)*3600" \
join=1and2 find=best ofmt=fits params=1 out=$datapath/XMM-Newton/Match_regalade/match_om_regalade

stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/XMM-Newton/SUSS6_ra_dec_per_src.fits \
in2=$datapath/Regalade/regalade_ra_dec.fits matcher=skyerr \
values1="RA+1/60 DEC 3*POSERR" values2="ra dec 0" \
join=1and2 find=best params=1 omode=count > $datapath/XMM-Newton/Match_regalade/count_random_match.txt

datapath="../../../Data"
echo "Performing match calculation UVOT-Regalade catalogue skyerr"
stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/Swift-UVOT/uvot_ra_dec_per_src.fits \
in2=$datapath/Regalade/regalade_ra_dec.fits matcher=skyerr \
values1="RA DEC 3*POSERR" values2="ra dec 0" \
ocmd="addcol angDist skyDistanceDegrees(RA,DEC,ra,dec)*3600" \
join=1and2 find=best ofmt=fits params=1 out=$datapath/Swift-UVOT/Match_regalade/match_om_regalade

stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/Swift-UVOT/uvot_ra_dec_per_src.fits \
in2=$datapath/Regalade/regalade_ra_dec.fits matcher=skyerr \
values1="RA+1/60 DEC 3*POSERR" values2="ra dec 0" \
join=1and2 find=best params=1 omode=count > $datapath/Swift-UVOT/Match_regalade/count_random_match.txt
