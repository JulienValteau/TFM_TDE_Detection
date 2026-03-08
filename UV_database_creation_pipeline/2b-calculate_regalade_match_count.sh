#!/usr/bin/bash

datapath="../../../Data"
echo "Random match count" > $datapath/XMM-Newton/Match_regalade/count_random_match.txt
echo "Computing random match values"

echo "Dividing the table in various part"
python Utilities/divide_time_table.py "Regalade" "regalade_ra_dec.fits"

echo "Match of each time table with the main catalogue"
j=0
for f in $( ls $datapath/Regalade/temp_part*.fits); do
    stilts tmatch2 ifmt1=fits ifmt2=fits in1=$datapath/XMM-Newton/SUSS6_ra_dec_per_src_rand.fits \
    in2=$f matcher=skyerr \
    values1="RA DEC 3*POSERR" values2="ra dec 0" \
    join=1and2 find=best params=1 omode=count >> $datapath/XMM-Newton/Match_regalade/count_random_match.txt
	let j=j+1
done

rm $datapath/Regalade/temp_*.fits
